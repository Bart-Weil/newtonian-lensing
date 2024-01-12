#ifndef OBJECT_H
#define OBJECT_H

#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"
#include "ray.h"

struct object_t;

typedef status_t (*renderer)(int *rgb, struct object_t object, vec3 pos_1, vec3 pos_2);
typedef void (*texture)(int *rgb, vec3 p, struct object_t object);
typedef void (*init)(struct object_t *object);
typedef void (*term)(struct object_t *object);

typedef struct object_t {
	vec3 pos;
	vec3 rot;
	texture obj_texture;
	renderer obj_renderer;
	init obj_init;
	term obj_term;
	mat3 rot_mat;
	int *tex;
} object_t;

extern void plane_circle_r_20_texture(int *rgb, vec3 p, object_t plane);
extern status_t plane_renderer(int *rgb, object_t plane, vec3 pos_1, vec3 pos_2);
extern void no_texture_plane_init(object_t *plane);
extern void no_texture_plane_term(object_t *plane);

#endif
