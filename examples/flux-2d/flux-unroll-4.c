#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void flux_calc_opt (double *t_vol_flux_x, double *t_vol_flux_y, double *t_xarea, double *t_yarea, double *t_xvel0, double *t_xvel1, double *t_yvel0, double *t_yvel1, double dt, int N) {
	double (*vol_flux_x)[4098] = (double (*)[4098]) t_vol_flux_x;
	double (*vol_flux_y)[4098] = (double (*)[4098]) t_vol_flux_y;
	double (*xarea)[4098] = (double (*)[4098]) t_xarea;
	double (*yarea)[4098] = (double (*)[4098]) t_yarea;
	double (*xvel0)[4098] = (double (*)[4098]) t_xvel0;
	double (*xvel1)[4098] = (double (*)[4098]) t_xvel1;
	double (*yvel0)[4098] = (double (*)[4098]) t_yvel0;
	double (*yvel1)[4098] = (double (*)[4098]) t_yvel1;

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
			for (j = 1; j < N-1; j+=4) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					vol_flux_x[j][i]=0.25*dt*xarea[j][i] *(xvel0[j][i]+xvel0[j+1][i]+xvel1[j][i]+xvel1[j+1][i]);
					vol_flux_x[j+1][i]=0.25*dt*xarea[j+1][i] *(xvel0[j+1][i]+xvel0[j+1+1][i]+xvel1[j+1][i]+xvel1[j+1+1][i]);
					vol_flux_x[j+2][i]=0.25*dt*xarea[j+2][i] *(xvel0[j+2][i]+xvel0[j+2+1][i]+xvel1[j+2][i]+xvel1[j+2+1][i]);
					vol_flux_x[j+3][i]=0.25*dt*xarea[j+3][i] *(xvel0[j+3][i]+xvel0[j+3+1][i]+xvel1[j+3][i]+xvel1[j+3+1][i]);
				}
			}
#pragma omp for private(i)
			for (j = 1; j < N-1; j+=4) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					vol_flux_y[j][i]=0.25*dt*yarea[j][i] *(yvel0[j][i]+yvel0[j][i+1]+yvel1[j][i]+yvel1[j][i+1]);
					vol_flux_y[j+1][i]=0.25*dt*yarea[j+1][i] *(yvel0[j+1][i]+yvel0[j+1][i+1]+yvel1[j+1][i]+yvel1[j+1][i+1]);
					vol_flux_y[j+2][i]=0.25*dt*yarea[j+2][i] *(yvel0[j+2][i]+yvel0[j+2][i+1]+yvel1[j+2][i]+yvel1[j+2][i+1]);
					vol_flux_y[j+3][i]=0.25*dt*yarea[j+3][i] *(yvel0[j+3][i]+yvel0[j+3][i+1]+yvel1[j+3][i]+yvel1[j+3][i+1]);
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("unroll: %6lf\n", (double)4096*4096*12*10/(end_time - start_time)/1e9);
}
