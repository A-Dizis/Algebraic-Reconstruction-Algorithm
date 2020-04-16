#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <string.h>
#include "functions.h"

int main(void)
{
	int N;				   //number of pixel analysis(NxN)
	int angles;			   //number of available angles
	double thetamax;	   //the max angle in which projections are taken(in degrees)
	int projections;	   //number of available projections
	FILE *fp_to_sinogram;  // file pointer to sinogram data file
	double linear_density; // number of points on an edge of a cell
	int linear_elements;   //elements along N edges
	int current_angle = 0;
	pthread_t thr[64]; //available threads
	int thr_number;	   //number of threads
	int i, j, k;
	double *sinogram; //sinogram matrix
	double *point_positions_matrix;
	double *Area, *U, *I;
	double *comparison_I, comparison_index;
	int compare_frequency;
	int c = 0;
	clock_t t1, t2;

	get_info(&N, &angles, &thetamax, &projections);
	sinogram = (double *)malloc(angles * projections * sizeof(double)); //allocate matrix for sinogram data
	get_sinogram(&fp_to_sinogram, angles, projections, &sinogram);
	//show_sinogram(&sinogram, angles, projections);
	get_density(&linear_density);
	get_comparison(&comparison_index, &compare_frequency);
	linear_elements = lround(linear_density * (double)N);
	point_positions_matrix = (double *)malloc(2 * sizeof(double) * linear_elements * linear_elements);
	create_point_positions(&point_positions_matrix, linear_elements, linear_density, N);
	//write_points(&point_positions_matrix, linear_elements);
	get_thread_number(&thr_number);

	param p[thr_number];
	Area = (double *)malloc(sizeof(double) * angles * projections);
	U = (double *)malloc(sizeof(double) * N * N * projections);
	I = (double *)malloc(sizeof(double) * N * N);
	comparison_I = (double *)malloc(sizeof(double) * N * N);

	initialize_image(N, &I);

repeat:

	t1 = clock();
	if(current_angle >= angles)
		current_angle = 0;
	while (current_angle < angles)
	{

		for (i = 0; i < thr_number; i++)
		{
			(*(p + i)).jump = i;
			(*(p + i)).angle = current_angle;
			(*(p + i)).thr_number = thr_number;
			(*(p + i)).thetamax = thetamax;
			(*(p + i)).angles = angles;
			(*(p + i)).projections = projections;
			(*(p + i)).N = N;
			(*(p + i)).U = &U;
			(*(p + i)).Area = &Area;
			(*(p + i)).linear_elements = linear_elements;
			(*(p + i)).point_positions_matrix = &point_positions_matrix;
			pthread_create(thr + i, NULL, make_projection, (void *)(p + i));
		}
		for (i = 0; i < thr_number; i++)
		{
			pthread_join(*(thr + i), NULL);
		}
		normalize_projection(&U, N, projections, (double)(linear_density) * (double)(linear_density));

		//print_projection_matrices(&U, N, projections);

		contribution(projections, N, current_angle, &Area, &U, &sinogram, &I, (double)(linear_density * linear_density));
		if (current_angle % compare_frequency == 0 && current_angle != 0)
		{
			printf("COMPARISON\n");
			if (converged(I, &comparison_I, N, comparison_index))
				break;
		}

		printf("Current angle : %d\n", current_angle);
		current_angle++;
		//create_bmp(N, N, &I);
	}

	t2 = clock();
	printf("CPU time : %lf\n", (double)(t2 - t1) / (double)CLOCKS_PER_SEC);
	c = ask_for_repetition();
	if (c == 1)
		goto repeat;

	normalize_image(N, &I);
	output_image(N, &I);

	return 0;
}
