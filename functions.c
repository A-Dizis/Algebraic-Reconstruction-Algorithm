
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

void get_comparison(double *comparison_index, int *compare_frequency)
{
	printf("Percentage of mean change between images: \n");

	fflush(stdin);
	scanf("%lf", comparison_index);
	getchar();

	printf("Compare images after # turns: \n");

	fflush(stdin);
	scanf("%d", compare_frequency);
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

int converged(double *I, double **comparison_I, int N, double comparison_index)
{
	double S = 0;
	double Q = 0;
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			S += fabs(*(I+i*N+j) - *(*comparison_I+i*N+j));
			Q += fabs(*(I+i*N+j)*comparison_index);
			*(*comparison_I+i*N+j) = *(I+i*N+j);
		}
	}
	printf("%lf\t%lf\n", S, Q);

	if(S <= Q)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	
}

int ask_for_repetition(void)
{
	char c;

	printf("Image has converged or a cycle has been completed. Continue? (type y for yes, another key for no) \n");
	fflush(stdin);
	if((c=getchar())!='y')
	{
		getchar();
		return 0;
	}
	else
	{
		getchar();
		return 1;
	}	
}

void normalize_image_bmp(int N, double **I, double *A)
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
			*(A + i * N + j) = *(*I + i * N + j) / max;
		}
	}
}

int create_bmp(int w, int h, double **red)
{
    FILE *f;
    unsigned char *img = NULL;
    int filesize = 54 + 3 * h * w; //w is your image width, h is image height, both int
    int x,y,r,g,b;

    img = (unsigned char *)malloc(3 * w * h);
    memset(img, 0, 3 * w * h);
    double *A = (double*)malloc(w*w*sizeof(double));

    normalize_image_bmp(w, red, A);

    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < h; j++)
        {
            x = i;
            y = (h - 1) - j;
            r = *(A + i*w +j) * 255;
            /*g = green[i][j] * 255;
            b = blue[i][j] * 255;*/
            if (r > 255)
                r = 255;
            /*if (g > 255)
                g = 255;
            if (b > 255)
                b = 255;*/
            img[(x + y * w) * 3 + 2] = (unsigned char)(r);
            img[(x + y * w) * 3 + 1] = (unsigned char)(r);//g
            img[(x + y * w) * 3 + 0] = (unsigned char)(r);//b
        }
    }

    unsigned char bmpfileheader[14] = {'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0};
    unsigned char bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0};
    unsigned char bmppad[3] = {0, 0, 0};

    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[4] = (unsigned char)(w);
    bmpinfoheader[5] = (unsigned char)(w >> 8);
    bmpinfoheader[6] = (unsigned char)(w >> 16);
    bmpinfoheader[7] = (unsigned char)(w >> 24);
    bmpinfoheader[8] = (unsigned char)(h);
    bmpinfoheader[9] = (unsigned char)(h >> 8);
    bmpinfoheader[10] = (unsigned char)(h >> 16);
    bmpinfoheader[11] = (unsigned char)(h >> 24);

    f = fopen("img.bmp", "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    for (int i = 0; i < h; i++)
    {
        fwrite(img + (w * (h - i - 1) * 3), 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }

    free(img);
    fclose(f);

    return 0;
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
	create_bmp(N, N, I);

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

