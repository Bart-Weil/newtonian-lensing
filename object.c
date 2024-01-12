#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "linalg.h"
#include "object.h"
#include "tracing.h"

#define CHANNEL 3

#define EPS 1e-10

static void set_rgb(int *rgb, int *col);

int red[CHANNEL] = {192, 57, 43}; 
int green[CHANNEL] = {22, 160, 133}; 
int blue[CHANNEL] = {133, 193, 233}; 
int yellow[CHANNEL] = {243, 156, 18}; 
int grey[CHANNEL] = {200, 200, 200}; 
int black[CHANNEL] = {0, 0, 0}; 

void plane_circle_r_20_texture(int *rgb, vec3 p, object_t plane) {
	vec3 d;
	scalar_multi_vec3(d, -1.0, plane.pos);
	add_vec3(d, p, d);

	printf("%f\n", sq_magnitude_vec3(d));

	if (sq_magnitude_vec3(d) < 400) {
		set_rgb(rgb, grey);
	} else {
		set_rgb(rgb, red);
	}
}

status_t plane_renderer(int *rgb, object_t plane, vec3 pos_1, vec3 pos_2) {
	// plane by default lies in the xy plane with a normal on the +z axis

	vec3 l;
	scalar_multi_vec3(l, -1.0, pos_1);
	add_vec3(l, pos_2, l);

	vec3 n_pre = {0.0, 0.0, 1.0};
	vec3 n;
	multi_mat3_vec3(n, plane.rot_mat, n_pre);

	float denom = dotp_vec3(n, l);

	if (fabs(denom) < EPS) { return RUNNING; }

	vec3 d_num;
	scalar_multi_vec3(d_num, -1.0, pos_1);
	add_vec3(d_num, plane.pos, d_num);

	float num = dotp_vec3(d_num, n);

	if (fabs(num) < EPS) {
		(*(plane.obj_texture))(rgb, pos_1, plane);

		return SEEN;
	}

	float d = num/denom;

	if (d < 0.0 || d > 1.0) { return RUNNING; }

	vec3 p;

	scalar_multi_vec3(l, d, l);
	add_vec3(p, l, pos_1);

	(*(plane.obj_texture))(rgb, p, plane);

	return SEEN;
}

void no_texture_plane_init(object_t *plane) {
	get_3D_rotation_matrix(plane->rot_mat, plane->rot);
	plane->tex = NULL;
}

void no_texture_plane_term(object_t *plane) { return; }

void set_rgb(int *rgb, int *col) {
	*rgb     = *col;
	*(rgb+1) = *(col+1);
	*(rgb+2) = *(col+2);
}

