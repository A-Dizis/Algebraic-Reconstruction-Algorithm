#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

void get_info(int *N, int *angles, double *thetamax, int *projections);
void get_sinogram(FILE **fp_to_sinogram, int angles, int projections, double **sinogram);
void show_sinogram(double **sinogram, int angles, int projections); //for testing purposes
void get_density(double *linear_density);
void create_point_positions(double **point_positions_matrix, int linear_elements, double linear_density, int N);
void write_points(double **point_positions_matrix, int linear_elements); //for testing purposes
void get_thread_number(int *thr_number);
void initialize_image(int N, double **I);
void *make_projection(void *t);
void normalize_projection(double **U, int N, int projections, double surface_density);
void print_projection_matrices(double **U, int N, int projections); //for testing purposes
void contribution(int projections, int N, int current_angle, double **Area, double **U, double **sinogram, double **I, double surface_density);
int ask_for_repetition(void);
void normalize_image(int N, double **I);
void output_image(int N, double **I);

struct parameters
{
	int jump;
	int N;
	double thetamax;
	int angles;
	int projections;
	int angle;
	int linear_elements;
	int thr_number;
	double **Area;
	double **U;
	double **point_positions_matrix;
};

typedef struct parameters param;

int main(void)
{
	int N;				   //number of pixel analysis(NxN)
	int angles;			   //number of available angles
	double thetamax;	   //the max angle in which projections are taken(in degrees)
	int projections;	   //number of available projections
	FILE *fp_to_sinogram;  // file pointer to sinogram data file
	double linear_density; // number of points on an edge of a cell
	int linear_elements;   //elements along N edges
	int current_angle;
	pthread_t thr[64]; //available threads
	int thr_number;	//number of threads
	int i, j, k;
	double *sinogram; //sinogram matrix
	double *point_positions_matrix;
	double *Area, *U, *I;
	int c = 0;
	clock_t t1, t2;

	get_info(&N, &angles, &thetamax, &projections);
	sinogram = (double *)malloc(angles * projections * sizeof(double)); //allocate matrix for sinogram data
	get_sinogram(&fp_to_sinogram, angles, projections, &sinogram);
	//show_sinogram(&sinogram, angles, projections);
	get_density(&linear_density);
	linear_elements = lround(linear_density * (double)N);
	point_positions_matrix = (double *)malloc(2 * sizeof(double) * linear_elements * linear_elements);
	create_point_positions(&point_positions_matrix, linear_elements, linear_density, N);
	//write_points(&point_positions_matrix, linear_elements);
	get_thread_number(&thr_number);

	param p[thr_number];
	Area = (double *)malloc(sizeof(double) * angles * projections);
	U = (double *)malloc(sizeof(double) * N * N * projections);
	I = (double *)malloc(sizeof(double) * N * N);

	initialize_image(N, &I);

repeat:

	t1 = clock();
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
		
		
		printf("Current angle : %d\n", current_angle);
		current_angle++;
	}
	t2 = clock();
	printf("CPU time : %lf\n", (double)(t2-t1)/(double)CLOCKS_PER_SEC);
	c = ask_for_repetition();
	if (c == 1)
		goto repeat;

	normalize_image(N, &I);
	output_image(N, &I);

	return 0;
}

void get_info(int *N, int *angles, double *thetamax, int *projections)
{
	int aux = 0;
	char answer = 'r';

	do
	{
		printf("Available angles:\n");
		fflush(stdin);
		scanf("%d", &(*angles));
		getchar();
		printf("Last angle position in degrees:\n");
		fflush(stdin);
		scanf("%lf", &(*thetamax));
		getchar();
		printf("Available projections:\n");
		fflush(stdin);
		scanf("%d", &(*projections));
		getchar();

		aux = lround((double)(*projections) / sqrt(2.0));
		if (aux % 2 == 1)
			aux = aux - 1;

		printf("Pixel analysis N for NxN image: (Recommended %d)\n", aux);
		fflush(stdin);
		scanf("%d", &(*N));
		getchar();
		printf("You inserted:\n______________\n");
		printf("Angles : %d\nLast angle position : %lf degrees\nProjections : %d\nPixel analysis N for NxN image : %d\n_____________\nIs this ok? (y to confirm)\n", *angles, *thetamax, *projections, *N);
		fflush(stdin);
		answer = getchar();
		getchar();
		fflush(stdin);

	} while (answer != 'y');
}

void get_sinogram(FILE **fp_to_sinogram, int angles, int projections, double **sinogram)
{
	char c = 'n';
	int i, j;

	printf("We are retrieving sinogram from files: (please put your sinogram in the program folder under the name sinogram.dat)\n_____________\nWhen ready click r, then enter.\n");

	while (c != 'r')
	{
		fflush(stdin);
		c = getchar();
		getchar();
		fflush(stdin);
	}
	*fp_to_sinogram = fopen("sinogram.dat", "r");

	if (*fp_to_sinogram == NULL)
	{
		printf("Sinogram did not open\n");
		exit(-1);
	}
	else
		printf("Succeesful opening\n");

	for (i = 0; i < angles; i++)
	{
		for (j = 0; j < projections; j++)
		{
			fscanf(*fp_to_sinogram, " %lf", (*sinogram + i * projections + j));
		}
		fgetc(*fp_to_sinogram);
	}
}

void show_sinogram(double **sinogram, int angles, int projections)
{
	int i, j;

	for (i = 0; i < angles; i++)
	{
		for (j = 0; j < projections; j++)
		{
			printf("%.2lf\t", *(*sinogram + i * projections + j));
		}
		printf("\n");
	}
}

void get_density(double *linear_density)
{
	printf("Linear density of points when calculating projection matrices(better stay under 6) : \n");

	fflush(stdin);
	scanf("%lf", linear_density);
	getchar();
}

void create_point_positions(double **point_positions_matrix, int linear_elements, double linear_density, int N)
{
	int i, j;

	for (i = 0; i < linear_elements; i++) //x
	{
		for (j = 0; j < linear_elements; j++) //y
		{
			*((*point_positions_matrix) + 2 * i * linear_elements + 2 * j) = -((double)N / 2.0) + 0.5 / linear_density + (double)i / linear_density;	 //x
			*((*point_positions_matrix + 1) + 2 * i * linear_elements + 2 * j) = -((double)N / 2.0) + 0.5 / linear_density + (double)j / linear_density; //y
		}
	}
}

void write_points(double **point_positions_matrix, int linear_elements)
{
	FILE *fp;
	int i, j;

	fp = fopen("elements.dat", "w");

	for (i = 0; i < linear_elements; i++)
	{
		for (j = 0; j < linear_elements; j++)
		{
			fprintf(fp, "(%.2lf, %.2lf)\t", *((*point_positions_matrix) + 2 * i * linear_elements + 2 * j), *((*point_positions_matrix + 1) + 2 * i * linear_elements + 2 * j));
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void get_thread_number(int *thr_number)
{
	printf("Threads to be used (better stay below 64) : \n");

	do
	{
		fflush(stdin);
		scanf("%d", thr_number);
		getchar();
	} while (!(((*thr_number) > 0) && ((*thr_number) <= 64)));
	printf("Threads are set to be %d\n", *thr_number);
}

void initialize_image(int N, double **I)
{
	int i, j;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			*(*I + i * N + j) = 0.0;
}

void *make_projection(void *t)
{
	param p = *((param *)t);

	int i, j, k;
	double temp[11];
	int aux[10];
	double half_projections = (double)p.projections / 2.0;
	double half_N = (double)p.N / 2.0;
	double helper = 0;

	const double scal = 0.0174532925 * (double)p.thetamax * (double)p.angle / (double)p.angles;

	if ((scal < 1.57079632) || (scal > 3.1415926 && scal < 4.71238898))
		helper = fabs(cos(0.7853981634 - scal) * (double)p.N) / sqrt(2.0);
	else
		helper = fabs(cos(2.356194490 - scal) * (double)p.N) / sqrt(2.0);
	
	//printf("Thread %d Angle %d\n", p.jump, p.angle);

	for (k = p.jump; k < p.projections; k = k + p.thr_number)
	{

		if (fabs(k - half_projections) <= helper)
		{
			*(*(p.Area) + p.angle * p.projections + k) = 0.0;

			for (i = 0; i < p.N; i++)
				for (j = 0; j < p.N; j++)
					*(*(p.U) + k * p.N * p.N + p.N * i + j) = 0.0;

			temp[1] = tan(scal);
			temp[2] = ((double)k - half_projections) / cos(scal);
			temp[3] = 1.0 / cos(scal);

			temp[5] = cos(scal)/sin(scal);
			temp[6] = -((double)k - half_projections) / sin(scal);
			temp[7] = -1.0 / sin(scal);

			if (scal <= 0.7853981634 ||  scal > 5.49777871438)
			{
				for (i = 0; i <= ((p.linear_elements) * (p.linear_elements)) - 1; i++)
				{
					if ((*(*(p.point_positions_matrix) + 2 * i + 1) > *(*(p.point_positions_matrix) + 2 * i) * temp[1] + temp[2]) && (*(*(p.point_positions_matrix) + 2 * i + 1) < *(*(p.point_positions_matrix) + 2 * i) * temp[1] + temp[2] + temp[3]))
					{
						aux[0] = lround(*(*(p.point_positions_matrix) + 2 * i) + half_N);
						aux[1] = lround(*(*(p.point_positions_matrix) + 2 * i + 1) + half_N);

						*(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) = *(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) + 1.0;
						*(*(p.Area) + p.angle * p.projections + k) = *(*(p.Area) + p.angle * p.projections + k) + 1.0;
					}
				}
			}
			else if (scal <= 2.3561944902)
			{
				for (i = 0; i <= ((p.linear_elements) * (p.linear_elements)) - 1; i++)
				{
					if ((*(*(p.point_positions_matrix) + 2 * i) < *(*(p.point_positions_matrix) + 2 * i + 1) * temp[5] + temp[6]) && (*(*(p.point_positions_matrix) + 2 * i) > *(*(p.point_positions_matrix) + 2 * i + 1) * temp[5] + temp[6] + temp[7]))
					{
						aux[0] = lround(*(*(p.point_positions_matrix) + 2 * i) + half_N);
						aux[1] = lround(*(*(p.point_positions_matrix) + 2 * i + 1) + half_N);

						*(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) = *(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) + 1.0;
						*(*(p.Area) + p.angle * p.projections + k) = *(*(p.Area) + p.angle * p.projections + k) + 1.0;
					}
				}
			}
			else if (scal <= 3.926990817)
			{
				for (i = 0; i <= ((p.linear_elements) * (p.linear_elements)) - 1; i++)
				{
					if ((*(*(p.point_positions_matrix) + 2 * i + 1) < *(*(p.point_positions_matrix) + 2 * i) * temp[1] + temp[2]) && (*(*(p.point_positions_matrix) + 2 * i + 1) > *(*(p.point_positions_matrix) + 2 * i) * temp[1] + temp[2] + temp[3]))
					{
						aux[0] = lround(*(*(p.point_positions_matrix) + 2 * i) + half_N);
						aux[1] = lround(*(*(p.point_positions_matrix) + 2 * i + 1) + half_N);

						*(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) = *(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) + 1.0;
						*(*(p.Area) + p.angle * p.projections + k) = *(*(p.Area) + p.angle * p.projections + k) + 1.0;
					}
				}
			}
			else
			{
				for (i = 0; i <= ((p.linear_elements) * (p.linear_elements)) - 1; i++)
				{
					if ((*(*(p.point_positions_matrix) + 2 * i) > *(*(p.point_positions_matrix) + 2 * i + 1) * temp[5] + temp[6]) && (*(*(p.point_positions_matrix) + 2 * i) < *(*(p.point_positions_matrix) + 2 * i + 1) * temp[5] + temp[6] + temp[7]))
					{
						aux[0] = lround(*(*(p.point_positions_matrix) + 2 * i) + half_N);
						aux[1] = lround(*(*(p.point_positions_matrix) + 2 * i + 1) + half_N);

						*(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) = *(*(p.U) + k * (p.N) * (p.N) + aux[1] * p.N + aux[0]) + 1.0;
						*(*(p.Area) + p.angle * p.projections + k) = *(*(p.Area) + p.angle * p.projections + k) + 1.0;
					}
				}
			}
		}
			else
		{
			*(*(p.Area) + p.angle * p.projections + k) = 0.0;

			for (i = 0; i < p.N; i++)
			{
				for (j = 0; j < p.N; j++)
				{
					*(*(p.U) + k * p.N * p.N + p.N * i + j) = 0.0;
				}
			}
		}
	}

	pthread_exit(NULL);
}

void print_projection_matrices(double **U, int N, int projections)
{
	int i, j, k;

	FILE *fp;

	fp = fopen("projection_matrices.dat", "a");

	for (k = 0; k < projections; k++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				fprintf(fp, "%lf\t", *(*U + k * N * N + i * N + j));
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}
	fclose(fp);
}

void normalize_projection(double **U, int N, int projections, double surface_density)
{
	int i, j, k;
	for (k = 0; k < projections; k++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				*(*(U) + k * N * N + i * N + j) = *(*(U) + k * N * N + i * N + j) / surface_density;
			}
		}
	}
}

void contribution(int projections, int N, int current_angle, double **Area, double **U, double **sinogram, double **I, double surface_density)
{
	double temp[2];
	int i, j, k;

	for (k = 0; k < projections; k++)
	{
		temp[0] = 0.0;

		if (*(*(Area) + current_angle * projections + k) > 0.0)
		{
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)
					temp[0] = temp[0] + (*(*(I) + i * N + j)) * (*(*(U) + k * N * N + i * N + j));
			temp[0] = *(*(sinogram) + projections * current_angle + k) - temp[0] / (*(*Area + current_angle * projections + k));

			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
					if (*(*(U) + k * N * N + i * N + j) > 0.0)
					{
						temp[1] = *(*(I) + i * N + j) + temp[0] * (*(*(U) + k * N * N + i * N + j)) * surface_density;
						if (temp[1] > 0.0)
						{
							*(*(I) + i * N + j) = temp[1];
						}
						else
						{
							*(*(I) + i * N + j) = 0.0;
						}
					}
				}
			}
		}
	}
}

int ask_for_repetition(void)
{
	int c;

	printf("Repeat? (type 1 for yes, another digit for no) \n");
	fflush(stdin);
	scanf("%d", &c);
	getchar();

	return c;
}

void normalize_image(int N, double **I)
{
	double max;
	int i, j;

	max = 0.0;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			if (max < (*(*I + i * N + j)))
				max = *(*I + i * N + j);
		}
	}
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			*(*I + i * N + j) = *(*I + i * N + j) / max;
		}
	}
}

void output_image(int N, double **I)
{
	int i, j;
	FILE *fp;

	fp = fopen("output_image.dat", "w");

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(fp, "%lf\t", *(*I + i * N + j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
