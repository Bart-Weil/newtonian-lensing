#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lensing.h"
#include "linalg.h"
#include "tracing.h"

#include <unistd.h>

const int red[3] = {192, 57, 43}; 
const int green[3] = {22, 160, 133}; 
const int blue[3] = {133, 193, 233}; 
const int yellow[3] = {243, 156, 18}; 
const int grey[3] = {200, 200, 200}; 
const int black[3] = {0, 0, 0}; 
const float line_width = 0.5;
const float line_space = 3;
float test_camera_loc[3] = {0, 0, 0};
float test_camera_orient[3] = {0, M_PI/2, 0};

void set_rgb_four_colour_plane(int *dest, float* coord) {
	// custom colour fetching goes here
	int y = ((int) fabs(*(coord+1)))/line_space;
	int z = ((int) fabs(*(coord+2)))/line_space;
	
	if (y%2==0 && z%2==0) {
		memcpy(dest, red, CHANNEL*sizeof(int));
	} else if (y%2==1 && z%2==0) {
		memcpy(dest, green, CHANNEL*sizeof(int));
	} else if (y%2==0 && z%2==1) {
		memcpy(dest, blue, CHANNEL*sizeof(int));
	} else {
		memcpy(dest, yellow, CHANNEL*sizeof(int));
	}
}

void set_rgb_grid_plane(int *dest, float* coord) {
	// custom colour fetching goes here
	float y = fabs(fmod(*(coord+1), (double) line_space));
	float z = fabs(fmod(*(coord+2), (double) line_space));
	
	if ((y<(line_width/2)) || (y>=(line_space-line_width/2))) {
		memcpy(dest, grey, CHANNEL*sizeof(int));
	} else if ((z<(line_width/2)) || (z>=(line_space-line_width/2))) {
		memcpy(dest, grey, CHANNEL*sizeof(int));
	} else {
		memcpy(dest, black, CHANNEL*sizeof(int));
	}
}

bool check_collision_grid_plane(float *v_1, float *v_2) {
	float x_1 = *v_1;
	float x_2 = *v_2;

	return ((x_1 < 100) && (x_2 >= 100));
}

void get_collision_grid_plane(float *dest, float *v_1, float *v_2) {
	float x_1 = *v_1;
	float y_1 = *(v_1+1);
	float z_1 = *(v_1+2);

	float x_2 = *v_2;
	float y_2 = *(v_2+1);
	float z_2 = *(v_2+2);

	float scale = (100-x_1)/(x_2-x_1);
	*dest = 100;
	*(dest+1) = y_1 + scale*(y_2-y_1);
	*(dest+2) = z_1 + scale*(z_2-z_1);
}

int *get_screen_loc(camera *target_camera, int x, int y) {
	return (target_camera->screen + (target_camera->pixel_x*CHANNEL)*y + CHANNEL*x);
}

void get_3D_rotation_matrix(float *dest, float *src) {
	// src = {yaw, pitch, roll} (yaw rotates radians in xy plane)
	// applies roll, then pitch, then yaw
	float sin_yaw = sin(*src);
	float cos_yaw = cos(*src);

	float sin_pitch = sin(*(src+1));
	float cos_pitch = cos(*(src+1));

	float sin_roll = sin(*(src+2));
	float cos_roll = cos(*(src+2));

	set_entry_NxK_float_matrix(dest, 0, 0, cos_yaw*cos_pitch, 
		D, D);
	set_entry_NxK_float_matrix(dest, 0, 1, cos_yaw*sin_pitch*sin_roll - sin_yaw*cos_roll, 
		D, D);
	set_entry_NxK_float_matrix(dest, 0, 2, cos_yaw*sin_pitch*cos_roll + sin_yaw*sin_roll,
		D, D);

	set_entry_NxK_float_matrix(dest, 1, 0, sin_yaw*cos_pitch,
		D, D);
	set_entry_NxK_float_matrix(dest, 1, 1, sin_yaw*sin_pitch*sin_roll + cos_yaw*cos_roll,
		D, D);
	set_entry_NxK_float_matrix(dest, 1, 2, sin_yaw*sin_pitch*cos_roll - cos_yaw*sin_roll,
		D, D);

	set_entry_NxK_float_matrix(dest, 2, 0, -sin_pitch,
		D, D);
	set_entry_NxK_float_matrix(dest, 2, 1, cos_pitch*sin_roll,
		D, D);
	set_entry_NxK_float_matrix(dest, 2, 2, cos_pitch*cos_roll,
		D, D);
}

ray *init_rays(camera *target_camera) {
	// by default camera orients in the z direction with up pointing in x direction;
	ray *rays = malloc(sizeof(ray) * target_camera->pixel_x * target_camera->pixel_y);
	float *screen_x_axis_pre = alloc_nD_float_array(D);
	float *screen_y_axis_pre = alloc_nD_float_array(D);
	float *screen_x_axis = alloc_nD_float_array(D);
	float *screen_y_axis = alloc_nD_float_array(D);
	*screen_x_axis_pre = 1;
	*(screen_y_axis_pre+1) = 1;

	float *camera_rotate = alloc_NxK_float_matrix(D, D);
	get_3D_rotation_matrix(camera_rotate, target_camera->orientation);
	multiply_NxK_float_matrix(screen_x_axis, camera_rotate, screen_x_axis_pre, D, D, D, 1);
	multiply_NxK_float_matrix(screen_y_axis, camera_rotate, screen_y_axis_pre, D, D, D, 1);

	float *screen_coord_pre = alloc_nD_float_array(D);
	float *screen_coord = alloc_nD_float_array(D);
	*screen_coord_pre = -((float) target_camera->screen_x)/2;
	*(screen_coord_pre+1) = -((float) target_camera->screen_y)/2;
	*(screen_coord_pre+2) = 1;

	multiply_NxK_float_matrix(screen_coord, camera_rotate, screen_coord_pre, D, D, D, 1);

	// stores intersection of current_ray with screen
	float *screen_pos = alloc_nD_float_array(D);
	// stores x and y components of internsection
	float *screen_pos_x = alloc_nD_float_array(D);
	float *screen_pos_y = alloc_nD_float_array(D);
	// coefficients for screen_x_axis, screen_y_axis to be added to screen_coord
	float x_ratio = (target_camera->screen_x)/((float) target_camera->pixel_x);
	float y_ratio = (target_camera->screen_y)/((float) target_camera->pixel_y);
	float x_axis_coeff, y_axis_coeff;

	ray *current_ray = rays;

	for (int y=0; y<target_camera->pixel_y; ++y) {
		y_axis_coeff = y_ratio * (y+0.5);
		scalar_multi_nD_float_array(screen_pos_y, y_axis_coeff, screen_y_axis, D);
		for (int x=0; x<target_camera->pixel_x; ++x) {
			x_axis_coeff = x_ratio * (x+0.5);
			scalar_multi_nD_float_array(screen_pos_x, x_axis_coeff, screen_x_axis, D);

			// stores velocity of current_ray
			float *v = alloc_nD_float_array(D);
			current_ray->coord = alloc_nD_float_array(D);

			add_nD_float_array(screen_pos, screen_coord, screen_pos_x, D);
			add_nD_float_array(screen_pos, screen_pos, screen_pos_y, D);
			add_nD_float_array(screen_pos, screen_pos, target_camera->location, D);
			unit_displacement_nD_float_array(v, screen_pos, target_camera->location, D);
			scalar_multi_nD_float_array(v, C, v, D);
			memcpy((current_ray)->coord, target_camera->location, D*sizeof(float));
			(current_ray++)->v = v;
		}
	}
	free(screen_x_axis_pre);
	free(screen_y_axis_pre);
	free(screen_coord_pre);

	free(screen_x_axis);
	free(screen_y_axis);
	free(camera_rotate);
	free(screen_coord);
	free(screen_pos);
	free(screen_pos_x);
	free(screen_pos_y);

	return rays;
}

int run_ray(ray *target_ray, camera *target_camera, int x, int y,
	render_object *render_objects, int n_objects) {

	float *old_coord = alloc_nD_float_array(D);
	float *collision = alloc_nD_float_array(D);
	int *rgb = get_screen_loc(target_camera, x, y);
	for (int step= 0;step<MAX_STEP;++step) {
		memcpy(old_coord, target_ray->coord, D*sizeof(float));
		if (update_ray(target_ray)==-1) {
			free(old_coord);
			free(collision);

			return -1;
		}
		for (int i=0; i<n_objects; ++i) {
			if ((*((render_objects+i)->check_collision))
				(old_coord, target_ray->coord)) {

				(*((render_objects+i)->get_collision))
					(collision, old_coord, target_ray->coord);
				
				(*((render_objects+i)->set_rgb))(rgb, collision);
				free(old_coord);
				free(collision);

				return 0;
			}
		}
	}
	free(old_coord);
	free(collision);

	return -1;
}

void render_camera_view(camera *target_camera, char *name) {
	FILE *f = fopen(name, "wb");
	fprintf(f, "P6\n%i %i 255\n", target_camera->pixel_x, target_camera->pixel_y);
	int *rgb;
	for (int y=0; y<target_camera->pixel_y; ++y) {
		for (int x=0; x<target_camera->pixel_x; ++x) {
			rgb = get_screen_loc(target_camera, x, y);
			fputc(*rgb, f);
			fputc(*(rgb+1), f);
			fputc(*(rgb+2), f);
		}
	}

	fclose(f);
}

void free_rays(ray *rays, int n_rays) {
	for (int i=0; i<n_rays; ++i) {
		free((rays+i)->v);
		free((rays+i)->coord);
	}
	free(rays);
}

int main() {
	camera *camera_1 = malloc(sizeof(camera));
	assert(camera_1!=NULL);

	camera_1->location = test_camera_loc;
	camera_1->orientation = test_camera_orient;
	camera_1->pixel_x = PIXEL_X;
	camera_1->pixel_y = PIXEL_Y;
	camera_1->screen_x = SCREEN_X;
	camera_1->screen_y = SCREEN_Y;

	int *test_screen = malloc(CHANNEL*(camera_1->pixel_x)*(camera_1->pixel_y)*sizeof(int));
	assert(test_screen!=NULL);

	camera_1->screen = test_screen;

	render_object *four_colour_plane = malloc(sizeof(render_object));

	// infinite plane parallel to yz plane with grey lines every 250m in
	// the x and z directions
	four_colour_plane->set_rgb = set_rgb_grid_plane;
	four_colour_plane->check_collision = check_collision_grid_plane;
	four_colour_plane->get_collision = get_collision_grid_plane;

	ray *rays = init_rays(camera_1);

	#pragma omp parallel for
	for (int i=0; i<(camera_1->pixel_x)*(camera_1->pixel_y); ++i) {
		printf("%f\n", ((float) i)/((camera_1->pixel_x)*(camera_1->pixel_y)));
		run_ray((rays+i), camera_1, i%(camera_1->pixel_x), i/(camera_1->pixel_x), four_colour_plane, 1);
	}

	render_camera_view(camera_1, "out.ppm");

	free(test_screen);
	free(four_colour_plane);
	free_rays(rays, (camera_1->pixel_x)*(camera_1->pixel_y));

	free(camera_1);

	return 0;
}