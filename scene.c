#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "linalg.h"
#include "object.h"
#include "tracing.h"

int main() {

	camera_t test_camera = {
		.pos = {0.0, 0.0, 0.0},
		.rot = {0.0, 0.0, 0.0},
		.pixel_x = 100,
		.pixel_y = 100,
		.fov = 120.0,
		.rpp = 1
	};

	camera_t *test_camera_p = &test_camera;

	object_t test_plane = {
		.pos = {-40.0, 0.0, 80.0},
		.rot = {0.0, 0.0, 0.0},
		.obj_texture = &plane_circle_r_20_texture,
		.obj_renderer = &plane_renderer,
		.obj_init = &no_texture_plane_init,
		.obj_term = &no_texture_plane_term
	};

	scene_t test_scene = {
		test_scene.scene_camera = test_camera_p,
		test_scene.objects = &test_plane,
		test_scene.n_objects = 1
	};

	scene_t *test_scene_p = &test_scene;

	init_scene(test_scene_p);
	init_tracing();

	render(test_scene_p);
	save_render(test_scene_p, "out.ppm");


	term_scene(test_scene_p);
	terminate_tracing();

	return 0;
}
