#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void j3d27pt_opt (double * __restrict__ t_in, double * __restrict__ t_out, double * __restrict__ t_c, int N) {
	double (*in)[514][514] = (double (*)[514][514]) t_in;
	double (*out)[514][514] = (double (*)[514][514]) t_out;
	double (*c)[9] = (double (*)[9]) t_c;

	int t, i, j, k, kk, jj, ii;
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
						out[k][j][i] = 0.1 * c[0][0]*in[k-1][j-1][i-1] + c[0][1]*in[k-1][j-1][i] + c[0][2]*in[k-1][j-1][i+1] +
							c[1][0]*in[k-1][j][i-1] +   c[1][1]*in[k-1][j][i] +   c[1][2]*in[k-1][j][i+1] +
							c[2][0]*in[k-1][j+1][i-1] + c[2][1]*in[k-1][j+1][i] + c[2][2]*in[k-1][j+1][i+1] +
							c[3][0]*in[k][j-1][i-1] +   c[3][1]*in[k][j-1][i] +   c[3][2]*in[k][j-1][i+1] +
							c[4][0]*in[k][j][i-1] +     c[4][1]*in[k][j][i] +     c[4][2]*in[k][j][i+1] +
							c[5][0]*in[k][j+1][i-1] +   c[5][1]*in[k][j+1][i] +   c[5][2]*in[k][j+1][i+1] +
							c[6][0]*in[k+1][j-1][i-1] + c[6][1]*in[k+1][j-1][i] + c[6][2]*in[k+1][j-1][i+1] +
							c[7][0]*in[k+1][j][i-1] +   c[7][1]*in[k+1][j][i] +   c[7][2]*in[k+1][j][i+1] +
							c[8][0]*in[k+1][j+1][i-1] + c[8][1]*in[k+1][j+1][i] + c[8][2]*in[k+1][j+1][i+1];
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
			for (k=1; k < N-1; k+=2) {
				for (j=1; j<N-1; j+=2) {
					for (i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll k=2,j=2,i=1 print-intrinsics true acc-size 1
						out[k][j][i] = c[0][0]*in[k-1][j-1][i-1] + c[0][1]*in[k-1][j-1][i] + c[0][2]*in[k-1][j-1][i+1] +
							c[1][0]*in[k-1][j][i-1] +   c[1][1]*in[k-1][j][i] +   c[1][2]*in[k-1][j][i+1] +
							c[2][0]*in[k-1][j+1][i-1] + c[2][1]*in[k-1][j+1][i] + c[2][2]*in[k-1][j+1][i+1] +
							c[3][0]*in[k][j-1][i-1] +   c[3][1]*in[k][j-1][i] +   c[3][2]*in[k][j-1][i+1] +
							c[4][0]*in[k][j][i-1] +     c[4][1]*in[k][j][i] +     c[4][2]*in[k][j][i+1] +
							c[5][0]*in[k][j+1][i-1] +   c[5][1]*in[k][j+1][i] +   c[5][2]*in[k][j+1][i+1] +
							c[6][0]*in[k+1][j-1][i-1] + c[6][1]*in[k+1][j-1][i] + c[6][2]*in[k+1][j-1][i+1] +
							c[7][0]*in[k+1][j][i-1] +   c[7][1]*in[k+1][j][i] +   c[7][2]*in[k+1][j][i+1] +
							c[8][0]*in[k+1][j+1][i-1] + c[8][1]*in[k+1][j+1][i] + c[8][2]*in[k+1][j+1][i+1];
#pragma end stencil1
	
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)512*512*512*53*5/(end_time - start_time)/1e9);
}
