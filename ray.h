#ifndef RAY_H
#define RAY_H

#include "linalg.h"

typedef enum {
	RUNNING,
	SEEN,
	ERR,
	TIMEOUT,
} status_t;

/* ray struct required by tracing.c */
typedef struct {
	vec3 pos;
	vec3 v;
} ray_t;

/* optional initialisation for tracing */
extern int init_tracing();
extern status_t update_ray(ray_t *tgt_ray);
/* optional termination for tracing */
extern int terminate_tracing();

#endif
