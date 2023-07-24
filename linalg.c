#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linalg.h"

float *alloc_nD_float_array(int d) {
	float *out_array = calloc(d, sizeof(float));
	assert(out_array!=NULL);

	return out_array;
}

void print_nD_float_array(float *v, int d) {
	printf("(");
	for (int i=0; i<(d-1); ++i) {
		printf("%f, ", *(v+i));
	}
	printf("%f)\n", *(v+d-1));
}

void add_nD_float_array(float *dest, float *v_1, float *v_2, int d) {
	for (int i=0; i<d; ++i) {
		*(dest+i) = *(v_1+i) + *(v_2+i);
	}
}

void scalar_multi_nD_float_array(float *dest, float c, float *v, int d) {
	for (int i=0; i<d; ++i) {
		*(dest+i) = *(v+i) * c;
	}
}

float magnitude_nD_float_array(float *v, int d) {
	float magnitude = 0;
	for (int i=0; i<d; ++i) {
		float entry = *(v+i);
		magnitude += entry * entry;
	}
	return sqrt(magnitude);
}

float euclidian_dist_nD_float_array(float *v_1, float *v_2, int d) {
	float dist = 0;
	float diff;
	for (int i=0; i<d; ++i) {
		diff = *(v_1+i) - *(v_2+i);
		dist += diff * diff;
	}

	return (float) sqrt(dist);
}

void unit_displacement_nD_float_array(float *dest, float *v_1, float *v_2, int d) {
	float inv_dist = 1.0/euclidian_dist_nD_float_array(v_1, v_2, d);
	float *negate_v_2 = alloc_nD_float_array(d);
	scalar_multi_nD_float_array(negate_v_2, -1, v_2, d);
	add_nD_float_array(dest, v_1, negate_v_2, d);
	scalar_multi_nD_float_array(dest, inv_dist, dest, d);
	free(negate_v_2);
}

float *alloc_NxK_float_matrix(int N, int K) {
	float *out_array = calloc(N*K, sizeof(float));
	assert(out_array!=NULL);

	return out_array;
}

void print_NxK_float_matrix(float *m, int N, int K) {
	for (int i=0; i<N*K; i+=K) {
		print_nD_float_array((m+i), K);
	}
}

void set_entry_NxK_float_matrix(float *dest, int i, int j, float x, int N, int K) {
	*(dest + i*K + j) = x;
}

float multiply_NxK_float_matrix_entry(float *m_1, float *m_2, int i, int j,
	int N_1, int K_1, int N_2, int K_2) {

	float sum = 0;
	int i_1, i_2;
	for (int k=0; k<K_1; ++k) {
		i_1 = K_1*i + k;
		i_2 = j + K_2*k;
		sum += *(m_1 + i_1) * *(m_2 + i_2);
	}

	return sum;
}

void multiply_NxK_float_matrix(float *dest, float *m_1, float *m_2, 
	int N_1, int K_1, int N_2, int K_2) {
	assert((dest!=m_1&&dest!=m_2));
	assert(K_1==N_2);
	for (int i=0; i<N_1; ++i) {
		for (int j=0; j<K_2; ++j) {
			*(dest + K_2*i + j) = 
				multiply_NxK_float_matrix_entry(m_1, m_2, i, j, 
					N_1, K_1, N_2, K_2);
		}
	}
}
