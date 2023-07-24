#define C 299792458l
#define D 3
#define G 0.000000000066743

#define STEP 0.000000001l

#define N_BODIES 2

#define k_schwarzchild 2*G/(C*C)

typedef struct {
	double mass;
	float schwarzchild_radius;
	float *coord;
} body;

typedef struct {
	float *coord;
	float *v;
} ray;

typedef body** continuum;

extern ray *init_ray(float *coord, float *v);
extern int update_ray(ray *target_ray);
