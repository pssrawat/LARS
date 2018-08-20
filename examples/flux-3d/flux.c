#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void flux_calc_opt (double *vol_flux_x, double *vol_flux_y, double *vol_flux_z, double *xarea, double *yarea, double *zarea, double *xvel0, double *xvel1, double *yvel0, double *yvel1, double *zvel0, double *zvel1, double dt, int N);

void cacheflush (void) {
	int N = 26;
	size_t n = 1 << N;
	int i, t;
	double *a = malloc(sizeof(double)*(n));
	double sum = 0.0;

	for(i=0; i<n; ++i)
		a[i] = (i*i);

	for(t=0; t<16; ++t)
		for(i=0; i<n; ++i)
			sum += (a[i]);
	sum = sqrt (sum);
}

int main (void) {
	cacheflush ();
	int N = 258;
	double (*vol_flux_x_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*vol_flux_y_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*vol_flux_z_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*vol_flux_x)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*vol_flux_y)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*vol_flux_z)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);


	double (*xvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xvel1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double dt = 0.856;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						vol_flux_x[k][j][i] = 0.1 * 0.125 * dt * xarea[k][j][i] * (xvel0[k][j][i] + xvel0[k][j+1][i] + xvel0[k+1][j][i] + xvel0[k+1][j+1][i] + xvel1[k][j][i] + xvel1[k][j+1][i] + xvel1[k+1][j][i] + xvel1[k+1][j+1][i]);
					}
				}
			}
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						vol_flux_y[k][j][i] = 0.1 * 0.125 * dt * yarea[k][j][i] * (yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k+1][j][i] + yvel0[k+1][j][i+1] + yvel1[k][j][i] + yvel1[k][j][i+1] + yvel1[k+1][j][i] + yvel1[k+1][j][i+1]);
					}
				}
			}
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						vol_flux_z[k][j][i] = 0.1 * 0.125 * dt * zarea[k][j][i] * (zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k][j][i+1] + zvel0[k][j+1][i+1] + zvel1[k][j][i] + zvel1[k][j][i+1] + zvel1[k][j+1][i] + zvel1[k][j+1][i+1]);
					}
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						vol_flux_x[k][j][i] = 0.125 * dt * xarea[k][j][i] * (xvel0[k][j][i] + xvel0[k][j+1][i] + xvel0[k+1][j][i] + xvel0[k+1][j+1][i] + xvel1[k][j][i] + xvel1[k][j+1][i] + xvel1[k+1][j][i] + xvel1[k+1][j+1][i]);
					}
				}
			}
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						vol_flux_y[k][j][i] = 0.125 * dt * yarea[k][j][i] * (yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k+1][j][i] + yvel0[k+1][j][i+1] + yvel1[k][j][i] + yvel1[k][j][i+1] + yvel1[k+1][j][i] + yvel1[k+1][j][i+1]);
					}
				}
			}
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						vol_flux_z[k][j][i] = 0.125 * dt * zarea[k][j][i] * (zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k][j][i+1] + zvel0[k][j+1][i+1] + zvel1[k][j][i] + zvel1[k][j][i+1] + zvel1[k][j+1][i] + zvel1[k][j+1][i+1]);
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*30*5/(end_time - start_time)/1e9);

	flux_calc_opt ((double*)vol_flux_x_opt, (double*)vol_flux_y_opt, (double*)vol_flux_z_opt, (double*)xarea, (double*)yarea, (double*)zarea, (double*)xvel0, (double*)xvel1, (double*)yvel0, (double*)yvel1, (double*)zvel0, (double*)zvel1, dt, N);


	double error = checkError3D (N, N, 0, (double*)vol_flux_x_opt, (double*)vol_flux_x, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
	error = checkError3D (N, N, 0, (double*)vol_flux_y_opt, (double*)vol_flux_y, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
	error = checkError3D (N, N, 0, (double*)vol_flux_z_opt, (double*)vol_flux_z, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
}
