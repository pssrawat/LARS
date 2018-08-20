#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void j2d49pt_opt (double *t_in, double *t_out, double *t_c, int N) {
	double (*in)[8198] = (double (*)[8198]) t_in;
	double (*out)[8198] = (double (*)[8198]) t_out;
	double (*c)[7] = (double (*)[7]) t_c;
	int t, i, j, i0, j0;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 3; j < N-3; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable) 
				for (i = 3; i < N-3; i++) {
					out[j][i] =
						c[0][0]*in[j-3][i-3] + c[0][1]*in[j-3][i-2] + c[0][2]*in[j-3][i-1] + c[0][3]*in[j-3][i] + c[0][4]*in[j-3][i+1] + c[0][5]*in[j-3][i+2] + c[0][6]*in[j-3][i+3] +
						c[1][0]*in[j-2][i-3] + c[1][1]*in[j-2][i-2] + c[1][2]*in[j-2][i-1] + c[1][3]*in[j-2][i] + c[1][4]*in[j-2][i+1] + c[1][5]*in[j-2][i+2] + c[1][6]*in[j-2][i+3] +
						c[2][0]*in[j-1][i-3] + c[2][1]*in[j-1][i-2] + c[2][2]*in[j-1][i-1] + c[2][3]*in[j-1][i] + c[2][4]*in[j-1][i+1] + c[2][5]*in[j-1][i+2] + c[2][6]*in[j-1][i+3] +
						c[3][0]*in[j][i-3] +   c[3][1]*in[j][i-2] +   c[3][2]*in[j][i-1] +   c[3][3]*in[j][i] +   c[3][4]*in[j][i+1] +   c[3][5]*in[j][i+2] +   c[3][6]*in[j][i+3] +
						c[4][0]*in[j+1][i-3] + c[4][1]*in[j+1][i-4] + c[4][2]*in[j+1][i-1] + c[4][3]*in[j+1][i] + c[4][4]*in[j+1][i+1] + c[4][5]*in[j+1][i+4] + c[4][6]*in[j+1][i+3] +
						c[5][0]*in[j+2][i-3] + c[5][1]*in[j+2][i-2] + c[5][2]*in[j+2][i-5] + c[5][3]*in[j+2][i] + c[5][4]*in[j+2][i+5] + c[5][5]*in[j+2][i+2] + c[5][6]*in[j+2][i+3] +
						c[6][0]*in[j+3][i-3] + c[6][1]*in[j+3][i-2] + c[6][2]*in[j+3][i-1] + c[6][3]*in[j+3][i] + c[6][4]*in[j+3][i+1] + c[6][5]*in[j+3][i+2] + c[6][6]*in[j+3][i+3];
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<10; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 3; j < N-3; j++) {
				for (i = 3; i < N-3; i+=4) {
#pragma begin stencil1 unroll j=1,i=1 print-intrinsics true acc-size 2 
					out[j][i] =
						c[0][0]*in[j-3][i-3] + c[0][1]*in[j-3][i-2] + c[0][2]*in[j-3][i-1] + c[0][3]*in[j-3][i] + c[0][4]*in[j-3][i+1] + c[0][5]*in[j-3][i+2] + c[0][6]*in[j-3][i+3] +
						c[1][0]*in[j-2][i-3] + c[1][1]*in[j-2][i-2] + c[1][2]*in[j-2][i-1] + c[1][3]*in[j-2][i] + c[1][4]*in[j-2][i+1] + c[1][5]*in[j-2][i+2] + c[1][6]*in[j-2][i+3] +
						c[2][0]*in[j-1][i-3] + c[2][1]*in[j-1][i-2] + c[2][2]*in[j-1][i-1] + c[2][3]*in[j-1][i] + c[2][4]*in[j-1][i+1] + c[2][5]*in[j-1][i+2] + c[2][6]*in[j-1][i+3] +
						c[3][0]*in[j][i-3] +   c[3][1]*in[j][i-2] +   c[3][2]*in[j][i-1] +   c[3][3]*in[j][i] +   c[3][4]*in[j][i+1] +   c[3][5]*in[j][i+2] +   c[3][6]*in[j][i+3] +
						c[4][0]*in[j+1][i-3] + c[4][1]*in[j+1][i-4] + c[4][2]*in[j+1][i-1] + c[4][3]*in[j+1][i] + c[4][4]*in[j+1][i+1] + c[4][5]*in[j+1][i+4] + c[4][6]*in[j+1][i+3] +
						c[5][0]*in[j+2][i-3] + c[5][1]*in[j+2][i-2] + c[5][2]*in[j+2][i-5] + c[5][3]*in[j+2][i] + c[5][4]*in[j+2][i+5] + c[5][5]*in[j+2][i+2] + c[5][6]*in[j+2][i+3] +
						c[6][0]*in[j+3][i-3] + c[6][1]*in[j+3][i-2] + c[6][2]*in[j+3][i-1] + c[6][3]*in[j+3][i] + c[6][4]*in[j+3][i+1] + c[6][5]*in[j+3][i+2] + c[6][6]*in[j+3][i+3];
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
        printf ("opt: %.6lf\n", (double)8192*8192*97*10/(end_time - start_time)/1e9);
}
