#define PIXEL_X 1920
#define PIXEL_Y 1080

#define SCREEN_X 1.92
#define SCREEN_Y 1.08

#define CHANNEL 3
#define C 299792458l
#define D 3

#define MAX_STEP 100000
#define N_OBJECTS 1

typedef void (*rgb_setter)(int *, float *);
typedef bool (*collision_checker)(float *, float *);
typedef void (*collision_getter)(float *, float *, float *);

typedef struct {
	rgb_setter set_rgb;
	collision_checker check_collision;
	collision_getter get_collision;
} render_object;

typedef struct {
	int *screen;
	float *location;
	float *orientation;
	int pixel_x;
	int pixel_y;
	float screen_x;
	float screen_y;
} camera;
