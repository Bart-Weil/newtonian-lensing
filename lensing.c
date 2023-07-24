#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "lensing.h"
#include "linalg.h"

float gravitational_acceleration(ray* target_ray, body* target_body) {
	float r = euclidian_dist_nD_float_array(target_body->coord, target_ray->coord, D);

	return G*(target_body->mass)/((r*r));
}

int update_position_single(ray* target_ray, body* target_body) {
	// check if ray is inside schwarzchild radius
	if (euclidian_dist_nD_float_array(target_ray->coord, target_body->coord, D)
		<= target_body->schwarzchild_radius) {
		return -1;
	}
	float accel_magnitude = 2 * gravitational_acceleration(target_ray, target_body);
	// store old velocity for SUVAT
	float *v_old = alloc_nD_float_array(D);
	memcpy(v_old, target_ray->v, D*sizeof(float));
	// calculate change in velocity
	float *accel_t = alloc_nD_float_array(D);
	unit_displacement_nD_float_array(accel_t,
		target_body->coord, target_ray->coord, D);
	scalar_multi_nD_float_array(accel_t, STEP*accel_magnitude, accel_t, D);
	add_nD_float_array(target_ray->v, target_ray->v, accel_t, D);
	// scale velocity to C
	float v_scale = magnitude_nD_float_array(target_ray->v, D);
	scalar_multi_nD_float_array(target_ray->v, C/v_scale, target_ray->v, D);
	// calculate change in position
	float *d_coord = alloc_nD_float_array(D);
	add_nD_float_array(d_coord, v_old, target_ray->v, D);
	scalar_multi_nD_float_array(d_coord, ((double) STEP)/2, d_coord, D);
	add_nD_float_array(target_ray->coord, target_ray->coord, d_coord, D);

	free(accel_t);
	free(v_old);
	free(d_coord);

	return 0;
}

int update_position_continuum(ray* target_ray, continuum target_continuum, int n_bodies) {
	for (int i=0; i<n_bodies; ++i) {
		if (update_position_single(target_ray, *(target_continuum+i)) == -1) {
			return -1;
		}
	}
	return 0;
}

float coord_1[3] = {70, 0, 0};
float coord_2[3] = {45, 0, 0};
float coord_3[3] = {65, 0, 0};

body body_1 = {
	.mass = 4e27l,
	.schwarzchild_radius = k_schwarzchild*4e27,
	.coord = coord_1
};

body body_2 = {
	.mass = 2e27l,
	.schwarzchild_radius = k_schwarzchild*2e27,
	.coord = coord_2
};
body body_3 = {
	.mass = 3e27l,
	.schwarzchild_radius = k_schwarzchild*3e27,
	.coord = coord_3
};

body *test_continuum[N_BODIES] = {&body_1, &body_2};

int update_ray(ray* target_ray) {

	return update_position_continuum(target_ray, test_continuum, N_BODIES);
}