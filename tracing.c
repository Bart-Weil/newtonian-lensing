#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linalg.h"
#include "object.h"
#include "ray.h"
#include "tracing.h"

#include <unistd.h>

#define CHANNEL 3
#define D 3

#define MAX_STEP 100000

#define FOV 120.0f

#define BLUR 0.5f
#define RPP 5

#define UNIFORM(x) (float)rand()/(float)(RAND_MAX/x);

static int *get_screen_loc(camera_t *tgt_camera, int x, int y);
static void init_rays(ray_t *dest, camera_t *tgt_camera);
static status_t run_ray(ray_t *tgt_ray, scene_t *tgt_scene, int x, int y);

static int *get_screen_loc(camera_t *tgt_camera, int x, int y) {
	return (tgt_camera->screen + (tgt_camera->pixel_x*CHANNEL)*y + CHANNEL*x);
}

static void init_rays(ray_t *rays, camera_t *tgt_camera) {
	// by default camera orients in the z direction with up pointing in x direction;
	vec3 screen_x_axis_pre = {1.0, 0.0, 0.0};
	vec3 screen_y_axis_pre = {0.0, 1.0, 0.0};
	vec3 screen_x_axis = {0.0};
	vec3 screen_y_axis = {0.0};

	mat3 camera_rotate = {{0.0}};
	get_3D_rotation_matrix(camera_rotate, tgt_camera->rot);
	multi_mat3_vec3(screen_x_axis, camera_rotate, screen_x_axis_pre);
	multi_mat3_vec3(screen_y_axis, camera_rotate, screen_y_axis_pre);

	vec3 screen_coord;
	vec3 screen_coord_pre;

	float screen_x = tgt_camera->pixel_x /
		pow(10, ((int) log10(tgt_camera->pixel_x)));
	float screen_y = tgt_camera->pixel_y /
		pow(10, ((int) log10(tgt_camera->pixel_y)));

	*screen_coord_pre = -((float) screen_x)/2;
	*(screen_coord_pre+1) = -((float) screen_y)/2;
	*(screen_coord_pre+2) = tan(0.5*M_PI-(0.5*tgt_camera->fov*M_PI/180.0f));

	multi_mat3_vec3(screen_coord, camera_rotate, screen_coord_pre);

	// stores intersection of current_ray with screen
	vec3 screen_pos = {0.0};
	// stores x and y components of internsection
	vec3 screen_pos_x = {0.0};
	vec3 screen_pos_y = {0.0};
	// coefficients for screen_x_axis, screen_y_axis to be added to screen_coord
	float x_ratio = (screen_x)/((float) tgt_camera->pixel_x);
	float y_ratio = (screen_y)/((float) tgt_camera->pixel_y);
	float x_axis_coeff, y_axis_coeff;

	// stores components of blurred ray velocity
	vec3 v_init = {0.0};
	float x_fact = 0.0;
	float y_fact = 0.0;
	vec3 blur_x = {0.0};
	vec3 blur_y = {0.0};

	ray_t *current_ray = rays;

	for (int y=0; y<tgt_camera->pixel_y; ++y) {
		y_axis_coeff = y_ratio * (y+0.5);
		scalar_multi_vec3(screen_pos_y, y_axis_coeff, screen_y_axis);

		for (int x=0; x<tgt_camera->pixel_x; ++x) {
			x_axis_coeff = x_ratio * (x+0.5);
			scalar_multi_vec3(screen_pos_x, x_axis_coeff, screen_x_axis);

			add_vec3(screen_pos, screen_coord, screen_pos_x);
			add_vec3(screen_pos, screen_pos, screen_pos_y);
			add_vec3(screen_pos, screen_pos, tgt_camera->pos);
			unit_displacement_vec3(v_init, screen_pos, tgt_camera->pos);

			for (int r=0; r<tgt_camera->rpp; ++r) {
				memcpy(current_ray->pos, tgt_camera->pos, D*sizeof(float));

				memcpy(current_ray->v, v_init, D*sizeof(float));
				// Do not alter course of first ray (important if rpp == 1)
				if (r != 0) {
					x_fact = 0.5*tgt_camera->blur - UNIFORM(tgt_camera->blur);
					y_fact = 0.5*tgt_camera->blur - UNIFORM(tgt_camera->blur);
				}
				scalar_multi_vec3(blur_x, x_fact, screen_x_axis);
				scalar_multi_vec3(blur_y, y_fact, screen_y_axis);
				add_vec3(current_ray->v, current_ray->v, blur_x);
				add_vec3(current_ray->v, current_ray->v, blur_y);
				current_ray++;
			}
		}
	}
}

static status_t run_ray(ray_t *tgt_ray, scene_t *tgt_scene, int x, int y) {
	int *rgb = get_screen_loc(tgt_scene->scene_camera, x, y);

	vec3 old_pos = {0.0};

	status_t status;

	for (int step=0;step<MAX_STEP;++step) {
		memcpy(old_pos, tgt_ray->pos, D*sizeof(float));
		status = update_ray(tgt_ray);

		if (status != RUNNING) { return status; }

		for (int i=0; i<tgt_scene->n_objects; ++i) {
			object_t tgt_object = (tgt_scene->objects)[i];

			status = tgt_object.obj_renderer(rgb, tgt_object, old_pos, tgt_ray->pos);

			if (status != RUNNING) { return status; }
		}
	}

	return TIMEOUT;
}

void render(scene_t *tgt_scene) {
	camera_t *scene_cam = tgt_scene->scene_camera;

	int count = 0;

	ray_t *rays = malloc(sizeof(ray_t)*(scene_cam->pixel_x)*(scene_cam->pixel_y)*(scene_cam->rpp));
	assert(rays!=NULL);

	init_rays(rays, tgt_scene->scene_camera);

	int num_rays = (scene_cam->pixel_x)*(scene_cam->pixel_y)*(scene_cam->rpp);

	#pragma omp parallel for
	for (int i=0; i<num_rays; ++i) {
		int pixel_i = i/(scene_cam->rpp);
		run_ray((rays+i), tgt_scene, pixel_i%(scene_cam->pixel_x), pixel_i/(scene_cam->pixel_x));
		count++;
//		printf("%.2f\n", ((float) count)/num_rays);
	}

	free(rays);

	// Take average for blur
	for (int i=0; i<((scene_cam->pixel_x)*(scene_cam->pixel_y)*CHANNEL); ++i) {
		*(scene_cam->screen+i) = *(scene_cam->screen+i)/(scene_cam->rpp);
	}
}

void save_render(scene_t *tgt_scene, char *name) {
	camera_t *scene_cam = tgt_scene->scene_camera;
	FILE *f = fopen(name, "wb");
	fprintf(f, "P6\n%i %i 255\n", scene_cam->pixel_x, scene_cam->pixel_y);
	int *rgb;
	for (int y=0; y<scene_cam->pixel_y; ++y) {
		for (int x=0; x<scene_cam->pixel_x; ++x) {
			rgb = get_screen_loc(scene_cam, x, y);
			fputc(*rgb, f);
			fputc(*(rgb+1), f);
			fputc(*(rgb+2), f);
		}
	}

	fclose(f);
}

void init_scene(scene_t *tgt_scene) {
	camera_t *tgt_camera = tgt_scene->scene_camera;
	int pixel_x = tgt_camera->pixel_x;
	int pixel_y = tgt_camera->pixel_y;
	tgt_camera->screen = malloc(CHANNEL*pixel_x*pixel_y*sizeof(int));
	assert(tgt_scene->scene_camera->screen!=NULL);
	object_t *tgt_object;

	for (int i=0; i<tgt_scene->n_objects; ++i) {
		tgt_object = tgt_scene->objects+i;
		tgt_object->obj_init(tgt_object);
	}
}

void term_scene(scene_t *tgt_scene) {
	free(tgt_scene->scene_camera->screen);

	object_t *tgt_object;
	for (int i=0; i<tgt_scene->n_objects; ++i) {
		tgt_object = tgt_scene->objects+i;
		tgt_object->obj_term(tgt_object);
	}
}
