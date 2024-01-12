#include "linalg.h"
#include "object.h"

typedef struct {
	vec3 pos;
	vec3 rot;
	int pixel_x;
	int pixel_y;
	float fov;
	float blur;
	int rpp;
	int *screen;
} camera_t;

typedef struct {
	camera_t *scene_camera;
	object_t *objects;
	size_t n_objects;
} scene_t;

extern void render(scene_t *tgt_scene);

extern void save_render(scene_t *tgt_scene, char *name);

extern void init_scene(scene_t *tgt_scene);

extern void term_scene(scene_t *tgt_scene);
