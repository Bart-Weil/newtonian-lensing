#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "linalg.h"

void print_vec3(vec3 v) {
	printf("[");
	for (int i=0; i<(3-1); ++i) {
		printf("%f, ", v[i]);
	}
	printf("%f]\n", v[2]);
}

void add_vec3(vec3 r, vec3 v1, vec3 v2) {
	for (int i=0; i<3; ++i) {
		r[i] = v1[i] + v2[i];
	}
}

void scalar_multi_vec3(vec3 r, float c, vec3 v) {
	for (int i=0; i<3; ++i) {
		r[i] = v[i] * c;
	}
}

float dotp_vec3(vec3 v1, vec3 v2) {
	float r = 0.0;
	for (int i=0; i<3; ++i) {
		r += v1[i] * v2[i];
	}

	return r;
}

float sq_magnitude_vec3(vec3 v) { return dotp_vec3(v, v); }

float sq_euclidian_dist_vec3(vec3 v1, vec3 v2) {
	float dist = 0;
	float diff;
	for (int i=0; i<3; ++i) {
		diff = v1[i] - v2[i];
		dist += diff * diff;
	}

	return dist;
}

void unit_displacement_vec3(vec3 r, vec3 v1, vec3 v2) {
	assert(r!=v1 && r!=v2);
	float inv_dist = 1.0/sqrt(sq_euclidian_dist_vec3(v1, v2));
	scalar_multi_vec3(r, -1, v2);
	add_vec3(r, v1, r);
	scalar_multi_vec3(r, inv_dist, r);
}

void print_mat3(mat3 m) {
	for (int i=0; i<3; ++i) {
		print_vec3(&(m[i][0]));
	}
}

static float multi_mat3_entry(mat3 m1, mat3 m2, int i, int j) {
	float sum = 0;
	for (int k=0; k<3; ++k) {
		sum += m1[i][k] * m2[k][j];
	}

	return sum;
}

void multi_mat3_mat3(mat3 r, mat3 m1, mat3 m2) {
	assert(r!=m1 && r!=m2);
	for (int i=0; i<3; ++i) {
		for (int j=0; j<3; ++j) {
			r[i][j] = multi_mat3_entry(m1, m2, i, j);
		}
	}
}

void multi_mat3_vec3(vec3 r, mat3 m, vec3 v) {
	assert(r!=v);
	for (int i=0; i<3; ++i) {
		r[i] = dotp_vec3(v, m[i]);
	}
}

void get_3D_rotation_matrix(mat3 r, vec3 rotation) {
	// src = {yaw, pitch, roll} (yaw rotates radians in xy plane)
	// applies roll, then pitch, then yaw
	float sin_yaw = sin(rotation[0]);
	float cos_yaw = cos(rotation[0]);

	float sin_pitch = sin(rotation[1]);
	float cos_pitch = cos(rotation[1]);

	float sin_roll = sin(rotation[2]);
	float cos_roll = cos(rotation[2]);

	r[0][0] = cos_yaw*cos_pitch;
	r[0][1] = cos_yaw*sin_pitch*sin_roll - sin_yaw*cos_roll;
	r[0][2] = cos_yaw*sin_pitch*cos_roll + sin_yaw*sin_roll;

	r[1][0] = sin_yaw*cos_pitch;
	r[1][1] = sin_yaw*sin_pitch*sin_roll + cos_yaw*cos_roll;
	r[1][2] = sin_yaw*sin_pitch*cos_roll - cos_yaw*sin_roll;

	r[2][0] = -sin_pitch;
	r[2][1] = cos_pitch*sin_roll;
	r[2][2] = cos_pitch*cos_roll;
}
