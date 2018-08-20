#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void poisson_opt (double *t_out, double *t_in, int N) {
	double (*in)[514][514] = (double (*)[514][514]) t_in;
	double (*out)[514][514] = (double (*)[514][514]) t_out;

	int i, j, k, t;
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
						out[k][j][i] = 0.1 * 2.666*in[k][j][i] -
							(0.166*in[k-1][j][i] + 0.166*in[k+1][j][i] + 0.166*in[k][j-1][i] + 0.166*in[k][j+1][i] + 0.166*in[k][j][i+1] + 0.166*in[k][j][i-1])-
							(0.0833*in[k-1][j-1][i] + 0.0833*in[k+1][j-1][i] + 0.0833*in[k-1][j+1][i] + 0.0833*in[k+1][j+1][i] +
							 0.0833*in[k-1][j][i-1] + 0.0833*in[k+1][j][i-1] + 0.0833*in[k][j-1][i-1] + 0.0833*in[k][j+1][i-1] +
							 0.0833*in[k-1][j][i+1] + 0.0833*in[k+1][j][i+1] + 0.0833*in[k][j-1][i+1] + 0.0833*in[k][j+1][i+1]);
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
			for (k = 1; k < N-1; k+=2) {
				for (j = 1; j < N-1; j+=4) {
					for (i = 1; i < N-1; i+=1) {
#pragma begin stencil1 unroll k=2,j=4,i=1 print-intrinsics false acc-size 1
						out[k][j][i] = 2.666*in[k][j][i] -
							(0.166*in[k-1][j][i] + 0.166*in[k+1][j][i] + 0.166*in[k][j-1][i] + 0.166*in[k][j+1][i] + 0.166*in[k][j][i+1] + 0.166*in[k][j][i-1])-
							(0.0833*in[k-1][j-1][i] + 0.0833*in[k+1][j-1][i] + 0.0833*in[k-1][j+1][i] + 0.0833*in[k+1][j+1][i] +
							 0.0833*in[k-1][j][i-1] + 0.0833*in[k+1][j][i-1] + 0.0833*in[k][j-1][i-1] + 0.0833*in[k][j+1][i-1] +
							 0.0833*in[k-1][j][i+1] + 0.0833*in[k+1][j][i+1] + 0.0833*in[k][j-1][i+1] + 0.0833*in[k][j+1][i+1]);
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)512*512*512*21*5/(end_time  - start_time)/1e9);
}
