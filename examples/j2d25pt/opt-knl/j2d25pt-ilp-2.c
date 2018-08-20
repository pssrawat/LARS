#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void j2d25pt_opt (double *t_in, double *t_out, double *t_c, int N) {
	double (*in)[8196] = (double (*)[8196]) t_in;
	double (*out)[8196] = (double (*)[8196]) t_out;
	double (*c)[5] = (double (*)[5]) t_c;
	int t, i, j, j0, i0;
	double start_time, end_time;

	// Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 2; j < N-2; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 2; i < N-2; i++) {
					out[j][i] = 0.1 * c[0][0]*in[j-2][i-2] + c[0][1]*in[j-2][i-1] + c[0][2]*in[j-2][i] + c[0][3]*in[j-2][i+1] + c[0][4]*in[j-2][i+2] +
						c[1][0]*in[j-1][i-2] + c[1][1]*in[j-1][i-1] + c[1][2]*in[j-1][i] + c[1][3]*in[j-1][i+1] + c[1][4]*in[j-1][i+2] +
						c[2][0]*in[j][i-2] +   c[2][1]*in[j][i-1] +   c[2][2]*in[j][i] +   c[2][3]*in[j][i+1] +   c[2][4]*in[j][i+2] +
						c[3][0]*in[j+1][i-2] + c[3][1]*in[j+1][i-1] + c[3][2]*in[j+1][i] + c[3][3]*in[j+1][i+1] + c[3][4]*in[j+1][i+2] +
						c[4][0]*in[j+2][i-2] + c[4][1]*in[j+2][i-1] + c[4][2]*in[j+2][i] + c[4][3]*in[j+2][i+1] + c[4][4]*in[j+2][i+2];
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i,j)
			for (j = 2; j < N-2; j += 1) {
				for (i = 2; i < N-2; i+=8) {
#pragma begin stencil1 unroll j=1,i=1 print-intrinsics true acc-size 2 
		out[j][i] =  c[0][0]*in[j-2][i-2] + c[0][1]*in[j-2][i-1] + c[0][2]*in[j-2][i] + c[0][3]*in[j-2][i+1] + c[0][4]*in[j-2][i+2] +
						c[1][0]*in[j-1][i-2] + c[1][1]*in[j-1][i-1] + c[1][2]*in[j-1][i] + c[1][3]*in[j-1][i+1] + c[1][4]*in[j-1][i+2] +
						c[2][0]*in[j][i-2] +   c[2][1]*in[j][i-1] +   c[2][2]*in[j][i] +   c[2][3]*in[j][i+1] +   c[2][4]*in[j][i+2] +
						c[3][0]*in[j+1][i-2] + c[3][1]*in[j+1][i-1] + c[3][2]*in[j+1][i] + c[3][3]*in[j+1][i+1] + c[3][4]*in[j+1][i+2] +
						c[4][0]*in[j+2][i-2] + c[4][1]*in[j+2][i-1] + c[4][2]*in[j+2][i] + c[4][3]*in[j+2][i+1] + c[4][4]*in[j+2][i+2];
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %.6lf\n", (double)8192*8192*49*5/(end_time - start_time)/1e9);
}
