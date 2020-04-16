#include "functions.c"

void get_info(int *N, int *angles, double *thetamax, int *projections);
void get_sinogram(FILE **fp_to_sinogram, int angles, int projections, double **sinogram);
void show_sinogram(double **sinogram, int angles, int projections); //for testing purposes
void get_density(double *linear_density);
void get_comparison(double *comparison_index, int *compare_frequency);
void create_point_positions(double **point_positions_matrix, int linear_elements, double linear_density, int N);
void write_points(double **point_positions_matrix, int linear_elements); //for testing purposes
void get_thread_number(int *thr_number);
void initialize_image(int N, double **I);
void *make_projection(void *t);
void normalize_projection(double **U, int N, int projections, double surface_density);
void print_projection_matrices(double **U, int N, int projections); //for testing purposes
void contribution(int projections, int N, int current_angle, double **Area, double **U, double **sinogram, double **I, double surface_density);
int converged(double *I, double **comparison_I, int N, double comparison_index);
int ask_for_repetition(void);
void normalize_image_bmp(int N, double **I, double *A);
int create_bmp(int w, int h, double **red);
void normalize_image(int N, double **I);
void output_image(int N, double **I);