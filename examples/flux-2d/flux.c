#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void flux_calc_opt (double *vol_flux_x, double *vol_flux_y, double *xarea, double *yarea, double *xvel0, double *xvel1, double *yvel0, double *yvel1, double dt, int N);

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
	int N = 4098;
	double (*vol_flux_x_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*vol_flux_y_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*vol_flux_x)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*vol_flux_y)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);


	double (*xvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xvel1)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yvel1)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double dt = 0.856;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					vol_flux_x[j][i]=0.25*dt*xarea[j][i] *(xvel0[j][i]+xvel0[j+1][i]+xvel1[j][i]+xvel1[j+1][i]);
				}
			}
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					vol_flux_y[j][i]=0.25*dt*yarea[j][i] *(yvel0[j][i]+yvel0[j][i+1]+yvel1[j][i]+yvel1[j][i+1]);

				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<10; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					vol_flux_x[j][i]=0.25*dt*xarea[j][i] *(xvel0[j][i]+xvel0[j+1][i]+xvel1[j][i]+xvel1[j+1][i]);
				}
			}
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					vol_flux_y[j][i]=0.25*dt*yarea[j][i] *(yvel0[j][i]+yvel0[j][i+1]+yvel1[j][i]+yvel1[j][i+1]);

				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)4096*4096*12*10/(end_time - start_time)/1e9);

	flux_calc_opt ((double*)vol_flux_x_opt, (double*)vol_flux_y_opt, (double*)xarea, (double*)yarea, (double*)xvel0, (double*)xvel1, (double*)yvel0, (double*)yvel1, dt, N);

	double error = checkError2D (N, 0, (double*)vol_flux_x_opt, (double*)vol_flux_x, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
	error = checkError2D (N, 0, (double*)vol_flux_y_opt, (double*)vol_flux_y, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
}
