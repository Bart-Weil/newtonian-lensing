typedef float vec3[3];
typedef float mat3[3][3];

extern void print_vec3(vec3 v);

extern void add_vec3(vec3 r, vec3 v1, vec3 v2);

extern void scalar_multi_vec3(vec3 r, float c, vec3 v);

extern float dotp_vec3(vec3 v1, vec3 v2);

extern float sq_magnitude_vec3(vec3 v);

extern float sq_euclidian_dist_vec3(vec3 v1, vec3 v2);

extern void unit_displacement_vec3(vec3 r, vec3 v1, vec3 v2);

extern void print_mat3(mat3 m);

extern void multi_mat3_mat3(mat3 r, mat3 m1, mat3 m2);

extern void multi_mat3_vec3(vec3 r, mat3 m, vec3 v);

extern void get_3D_rotation_matrix(mat3 r, vec3 rotation);
