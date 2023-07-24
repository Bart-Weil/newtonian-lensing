extern float *alloc_nD_float_array(int d);

extern void print_nD_float_array(float *v, int d);

extern void add_nD_float_array(float *dest, float *v_1, float *v_2, int d);

extern void scalar_multi_nD_float_array(float *dest, float c, float *v, int d);

extern float magnitude_nD_float_array(float *v, int d);

extern float euclidian_dist_nD_float_array(float *v_1, float *v_2, int d);

extern void unit_displacement_nD_float_array(float *dest, float *v_1, float *v_2, int d);

extern float *alloc_NxK_float_matrix(int N, int K);

extern void print_NxK_float_matrix(float *m, int N, int K);

extern void set_entry_NxK_float_matrix(float *dest, int i, int j, float x, int N, int K);

extern void multiply_NxK_float_matrix(float *dest, float *m_1, float *m_2, int N_1, int K_1, int N_2, int K_2);
