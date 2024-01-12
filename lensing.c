#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ray.h"
#include "linalg.h"

#define C 299792458l
#define D 3
#define G 0.000000000066743

#define STEP 0.00000000001l

#define N_BODIES 1

#define S_RADIUS(m) m*2*G/(C*C)

typedef struct body{
	float mass;
	float s_radius;
	vec3 pos;
} body_t;

body_t continuum[N_BODIES];

/* Pass by value to compute on the stack without copying from memory */
static float grav_accel(ray_t tgt_ray, body_t tgt_body) {
	float r_2 = sq_euclidian_dist_vec3(tgt_body.pos, tgt_ray.pos);

	return 2 * G*(tgt_body.mass)/r_2;
}

status_t update_pos_continuum(ray_t* tgt_ray) {
	vec3 accel;
	vec3 v_old;
	vec3 d_pos;

	memcpy(v_old, tgt_ray->v, 3*sizeof(float));

	body_t tgt_body;
	for (int i=0; i<N_BODIES; ++i) {
		tgt_body = *(continuum+i);
		/* check if ray is inside schwarzchild radius */
		if (sq_euclidian_dist_vec3(tgt_ray->pos, tgt_body.pos)
			<= (tgt_body.s_radius)*(tgt_body.s_radius)) {
			return ERR;
		}

		/* compute resultant acceleration */
		float accel_magnitude = grav_accel(*tgt_ray, tgt_body);
		unit_displacement_vec3(accel, tgt_body.pos, tgt_ray->pos);
		scalar_multi_vec3(accel, STEP*accel_magnitude, accel);

		add_vec3(tgt_ray->v, tgt_ray->v, accel);
	}
	/* scale velocity to C */
	float v_scale = C/sqrt(sq_magnitude_vec3(tgt_ray->v));
	scalar_multi_vec3(tgt_ray->v, v_scale, tgt_ray->v);
	/* calculate change in position */
	add_vec3(d_pos, v_old, tgt_ray->v);
	scalar_multi_vec3(d_pos, ((double) STEP)/2, d_pos);
	add_vec3(tgt_ray->pos, tgt_ray->pos, d_pos);

	return RUNNING;
}

int init_tracing() {
	continuum[0] = (body_t) {.mass = 1e26, .s_radius = S_RADIUS(1e26), 
		.pos = {0,0,60}};
	return 0;
}

status_t update_ray(ray_t* tgt_ray) {
	return update_pos_continuum(tgt_ray);
}

int terminate_tracing() { return 0; }
