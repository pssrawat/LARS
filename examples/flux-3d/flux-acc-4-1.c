#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void flux_calc_opt (double *t_vol_flux_x, double *t_vol_flux_y, double *t_vol_flux_z, double *t_xarea, double *t_yarea, double *t_zarea, double *t_xvel0, double *t_xvel1, double *t_yvel0, double *t_yvel1, double *t_zvel0, double *t_zvel1, double dt, int N) {
	double (*vol_flux_x)[258][258] = (double (*)[258][258]) t_vol_flux_x;
	double (*vol_flux_y)[258][258] = (double (*)[258][258]) t_vol_flux_y;
	double (*vol_flux_z)[258][258] = (double (*)[258][258]) t_vol_flux_z;
	double (*xarea)[258][258] = (double (*)[258][258]) t_xarea;
	double (*yarea)[258][258] = (double (*)[258][258]) t_yarea;
	double (*zarea)[258][258] = (double (*)[258][258]) t_zarea;
	double (*xvel0)[258][258] = (double (*)[258][258]) t_xvel0;
	double (*xvel1)[258][258] = (double (*)[258][258]) t_xvel1;
	double (*yvel0)[258][258] = (double (*)[258][258]) t_yvel0;
	double (*yvel1)[258][258] = (double (*)[258][258]) t_yvel1;
	double (*zvel0)[258][258] = (double (*)[258][258]) t_zvel0;
	double (*zvel1)[258][258] = (double (*)[258][258]) t_zvel1;

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
			for (k = 1; k < N-1; k+=4) {
				for (j = 1; j < N-1; j++) {
					for (i = 1; i < N-1; i+=1) {
#pragma begin stencil1 unroll k=4,j=1,i=1 print-intrinsics false acc-size 1
						vol_flux_x[k][j][i] = 0.125 * dt * xarea[k][j][i] * (xvel0[k][j][i] + xvel0[k][j+1][i] + xvel0[k+1][j][i] + xvel0[k+1][j+1][i] + xvel1[k][j][i] + xvel1[k][j+1][i] + xvel1[k+1][j][i] + xvel1[k+1][j+1][i]);
						vol_flux_y[k][j][i] = 0.125 * dt * yarea[k][j][i] * (yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k+1][j][i] + yvel0[k+1][j][i+1] + yvel1[k][j][i] + yvel1[k][j][i+1] + yvel1[k+1][j][i] + yvel1[k+1][j][i+1]);
						vol_flux_z[k][j][i] = 0.125 * dt * zarea[k][j][i] * (zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k][j][i+1] + zvel0[k][j+1][i+1] + zvel1[k][j][i] + zvel1[k][j][i+1] + zvel1[k][j+1][i] + zvel1[k][j+1][i+1]);
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)256*256*256*30*5/(end_time - start_time)/1e9);
}
