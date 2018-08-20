#ifndef __COMMON_H__
#define __COMMON_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdbool.h>

#define TOLERANCE (1e-4) 

double rtclock () {
	struct timezone Tzp;
	struct timeval Tp;
	int stat = gettimeofday (&Tp, &Tzp);
	if (stat != 0) 
		printf ("Error return from gettimeofday: %d", stat);
	return (Tp.tv_sec + Tp.tv_usec*1.0e-6);
}  


static double get_random () {
	return ((double)(rand())/(double)(RAND_MAX-1));
}

static double* getRandom1DArray (size_t width_y) {
	double (*a) = (double (*))malloc (sizeof (double) * width_y);
	int j;
	for (j = 0; j < width_y; j++)
		a[j] = get_random ();
	return (double*)a;
}

static double* getRandom2DArray (int width_y, int width_x) {
	double (*a)[width_x] = (double (*)[width_x]) malloc (sizeof (double) * width_y * width_x);
	int j, k;
	for (j = 0; j < width_y; j++) {
		for (k = 0; k < width_x; k++) {
			a[j][k] = get_random ();
		}
	}
	return (double*)a;
}

static double* getRandom3DArray (int height, int width_y, int width_x) {
	double (*a)[width_y][width_x] = (double (*)[width_y][width_x]) malloc (sizeof (double) * height * width_y * width_x);
	int i, j, k;
	for (i = 0; i < height; i++)
		for (j = 0; j < width_y; j++)
			for (k = 0; k < width_x; k++) {
				a[i][j][k] = get_random ();
			}
	return (double*)a;
}

static double* getRandom4DArray (int d, int height, int width_y, int width_x) {
	double (*a)[height][width_y][width_x] = (double (*)[height][width_y][width_x]) malloc (sizeof (double) * d * height * width_y * width_x);
	int t, i, j, k;
        for (t = 0; t < d; t++) 
	for (i = 0; i < height; i++)
		for (j = 0; j < width_y; j++)
			for (k = 0; k < width_x; k++) {
				a[t][i][j][k] = get_random ();
			}
	return (double*)a;
}

static double* getZero1DArray (int width_x) {
	double (*a) = (double *) malloc (sizeof (double) * width_x);
	memset((void*)a, 0, sizeof(double) * width_x);
	return (double*)a;
}

static double* getZero2DArray (int width_y, int width_x) {
	double (*a)[width_x] = (double (*)[width_x]) malloc (sizeof (double) * width_y * width_x);
	memset((void*)a, 0, sizeof(double) * width_y * width_x);
	return (double*)a;
}

static double* getZero4DArray (int d, int height, int width_y, int width_x) {
	double (*a)[height][width_y][width_x] = (double (*)[height][width_y][width_x]) malloc (sizeof (double) * d * height * width_y * width_x);
	memset((void*)a, 0, sizeof(double) * d * height * width_y * width_x);
	return (double*)a;
}

static double* getZero3DArray (int height, int width_y, int width_x) {
	double (*a)[width_y][width_x] = (double (*)[width_y][width_x]) malloc (sizeof (double) * height * width_y * width_x);
	memset((void*)a, 0, sizeof(double) * height * width_y * width_x);
	return (double*)a;
}

static double checkError1D (int width_x, const double *output, const double *reference, int x_lb, int x_ub) {
	double error = 0.0;
	double max_error = TOLERANCE;
	int max_j = 0, j;
	for (j = x_lb; j < x_ub; j++) { 
		double curr_error = output[j] - reference[j];
		curr_error = (curr_error < 0.0 ? -curr_error : curr_error);
		error += curr_error * curr_error;
		if (curr_error > max_error) {
			//printf ("Values at index (%d) differ : %.6lf and %.6lf\n", j, reference[j], output[j]);
			max_error = curr_error;
			max_j = j;
		}
	}
	//printf ("Max Error (%d) : %e\n", max_j, max_error);
	error = sqrt(error / (x_ub - x_lb));
	return error;
}

static double checkError2D (int width_x, int pad, const double *l_output, const double *l_reference, int y_lb, int y_ub, int x_lb, int x_ub) {
	const double (*output)[width_x+pad] = (const double (*)[width_x+pad])(l_output);
	const double (*reference)[width_x+pad] = (const double (*)[width_x+pad])(l_reference);
	double error = 0.0;
	double max_error = TOLERANCE;
	int max_k = 0, max_j = 0, j, k;
	for (j = y_lb; j < y_ub; j++) 
		for (k = x_lb; k < x_ub; k++) {
			//printf ("Values at index (%d,%d) are %.6lf and %.6lf\n", j, k, reference[j][k], output[j][k]);
			double curr_error = output[j][k] - reference[j][k];
			curr_error = (curr_error < 0.0 ? -curr_error : curr_error);
			error += curr_error * curr_error;
			if (curr_error > max_error) {
				//printf ("Values at index (%d,%d) differ : %.6lf and %.6lf\n", j, k, reference[j][k], output[j][k]);
				max_error = curr_error;
				max_k = k;
				max_j = j;
			}
		}
	//printf ("%e, ", max_error);
	error = sqrt(error / ( (y_ub - y_lb) * (x_ub - x_lb)));
	return error;
}

static double checkError3D (int width_y, int width_x, int pad, const double *l_output, const double *l_reference, int z_lb, int z_ub, int y_lb, int y_ub, int x_lb, int x_ub) {
	const double (*output)[width_y+pad][width_x+pad] = (const double (*)[width_y+pad][width_x+pad])(l_output);
	const double (*reference)[width_y+pad][width_x+pad] = (const double (*)[width_y+pad][width_x+pad])(l_reference);
	double error = 0.0;
	double max_error = TOLERANCE;
	int max_k = 0, max_j = 0, max_i = 0, i, j, k;
	//bool first_found = false;
	int num_error = 0;
	for (i = z_lb; i < z_ub; i++)
		for (j = y_lb; j < y_ub; j++)
			for (k = x_lb; k < x_ub; k++) {
				//printf ("real var1[%d][%d][%d] = %.6lf and %.6lf\n", i, j, k, reference[i][j][k], output[i][j][k]);
				double curr_error = output[i][j][k] - reference[i][j][k];
				curr_error = (curr_error < 0.0 ? -curr_error : curr_error);
				error += curr_error * curr_error;
				if (curr_error > max_error) {
					//if (!first_found) {
					//	printf ("Values at index (%d,%d,%d) differ : %.6lf and %.6lf\n", i, j, k, reference[i][j][k], output[i][j][k]);
					//	first_found = true;
					//}
					num_error++;
					max_error = curr_error;
					max_k = k;
					max_j = j;
					max_i = i;
				}
			}
	//printf ("%e, ", max_error);
	//printf ("Num error = %d\n", num_error);
	error = sqrt(error / ( (z_ub - z_lb) * (y_ub - y_lb) * (x_ub - x_lb)));
	return error;
}

static double checkError4D (int height, int width_y, int width_x, int pad, const double *l_output, const double *l_reference, int d_lb, int d_ub, int z_lb, int z_ub, int y_lb, int y_ub, int x_lb, int x_ub) {
	const double (*output)[height+pad][width_y+pad][width_x+pad] = (const double (*)[height+pad][width_y+pad][width_x+pad])(l_output);
	const double (*reference)[height+pad][width_y+pad][width_x+pad] = (const double (*)[height+pad][width_y+pad][width_x+pad])(l_reference);
	double error = 0.0;
	double max_error = TOLERANCE;
	int max_k = 0, max_j = 0, max_i = 0, i, j, k, t;
         for (t = d_lb; t < d_ub; t++) 
	for (i = z_lb; i < z_ub; i++)
		for (j = y_lb; j < y_ub; j++)
			for (k = x_lb; k < x_ub; k++) {
				//printf ("real var1[%d][%d][%d] = %.6lf and %.6lf\n", i, j, k, reference[i][j][k], output[i][j][k]);
				double curr_error = output[t][i][j][k] - reference[t][i][j][k];
				curr_error = (curr_error < 0.0 ? -curr_error : curr_error);
				error += curr_error * curr_error;
				if (curr_error > max_error) {
					//printf ("Values at index (%d,%d,%d,%d) differ : %.6lf and %.6lf\n", t, i, j, k, reference[t][i][j][k], output[t][i][j][k]);
					max_error = curr_error;
					max_k = k;
					max_j = j;
					max_i = i;
				}
			}
	//printf ("%e, ", max_error);
	error = sqrt(error / ( (d_ub - d_lb) * (z_ub - z_lb) * (y_ub - y_lb) * (x_ub - x_lb)));
	return error;
}



#endif
