#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void j3d7pt_opt (double * __restrict__ t_out, double * __restrict__ t_in, int N) {
	double (*out)[514][514] = (double (*)[514][514]) t_out;
	double (*in)[514][514] = (double (*)[514][514]) t_in;
	double a = 1.141;
	double b = 0.121;
	double h2inv = b*6.215;

	int t, i, j, k, kk, jj, ii;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma clang loop vectorize (enable) interleave(enable)
#pragma GCC ivdep
					for (i = 1; i < N-1; i++) {
						out[k][j][i] = 0.1 * a*in[k][j][i] - h2inv*(in[k][j][i+1] + in[k][j][i-1] + in[k][j+1][i] + in[k][j-1][i] + in[k+1][j][i] + in[k-1][j][i]) -in[k][j][i]*6.0;
					}
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<10; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k=1; k < N-1; k+=1) {
				for (j=1; j<N-1; j+=1) {
					for (i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll k=1,j=1,i=1 print-intrinsics true acc-size 1 
						out[k][j][i] = a*in[k][j][i] - h2inv*(in[k][j][i+1] + in[k][j][i-1] + in[k][j+1][i] + in[k][j-1][i] + in[k+1][j][i] + in[k-1][j][i]) -in[k][j][i]*6.0;
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
        printf ("opt: %6lf\n", (double)512*512*512*11*10/(end_time - start_time)/1e9);
}
