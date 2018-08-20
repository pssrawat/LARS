#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void j2d64pt_opt (double *t_in, double *t_out, double *t_c, int N) {
	double (*in)[8200] = (double (*)[8200]) t_in;
	double (*out)[8200] = (double (*)[8200]) t_out;
	double (*c)[8] = (double (*)[8]) t_c;
	int t, i, j, i0, j0;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) { 
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 4; i < N-4; i++) {
					out[j][i] =  0.1 * 
						c[0][0]*in[j-4][i-4] +  c[0][1]*in[j-4][i-3] + c[0][2]*in[j-4][i-2] + c[0][3]*in[j-4][i-1] + c[0][4]*in[j-4][i+1] + c[0][5]*in[j-4][i+2] + c[0][6]*in[j-4][i+3] + c[0][7]*in[j-4][i+4] +
						c[1][0]*in[j-3][i-4] +  c[1][1]*in[j-3][i-3] + c[1][2]*in[j-3][i-2] + c[1][3]*in[j-3][i-1] + c[1][4]*in[j-3][i+1] + c[1][5]*in[j-3][i+2] + c[1][6]*in[j-3][i+3] + c[1][7]*in[j-3][i+4] +
						c[2][0]*in[j-2][i-4] +  c[2][1]*in[j-2][i-3] + c[2][2]*in[j-2][i-2] + c[2][3]*in[j-2][i-1] + c[2][4]*in[j-2][i+1] + c[2][5]*in[j-2][i+2] + c[2][6]*in[j-2][i+3] + c[2][7]*in[j-2][i+4] +
						c[3][0]*in[j-1][i-4] +  c[3][1]*in[j-1][i-3] + c[3][2]*in[j-1][i-2] + c[3][3]*in[j-1][i-1] + c[3][4]*in[j-1][i+1] + c[3][5]*in[j-1][i+2] + c[3][6]*in[j-1][i+3] + c[3][7]*in[j-1][i+4] +
						c[4][0]*in[j+1][i-4] +  c[4][1]*in[j+1][i-3] + c[4][2]*in[j+1][i-2] + c[4][3]*in[j+1][i-1] + c[4][4]*in[j+1][i+1] + c[4][5]*in[j+1][i+2] + c[4][6]*in[j+1][i+3] + c[4][7]*in[j+1][i+4] +
						c[5][0]*in[j+2][i-4] +  c[5][1]*in[j+2][i-3] + c[5][2]*in[j+2][i-2] + c[5][3]*in[j+2][i-1] + c[5][4]*in[j+2][i+1] + c[5][5]*in[j+2][i+2] + c[5][6]*in[j+2][i+3] + c[5][7]*in[j+2][i+4] +
						c[6][0]*in[j+3][i-4] +  c[6][1]*in[j+3][i-3] + c[6][2]*in[j+3][i-2] + c[6][3]*in[j+3][i-1] + c[6][4]*in[j+3][i+1] + c[6][5]*in[j+3][i+2] + c[6][6]*in[j+3][i+3] + c[6][7]*in[j+3][i+4] +
						c[7][0]*in[j+4][i-4] +  c[7][1]*in[j+4][i-3] + c[7][2]*in[j+4][i-2] + c[7][3]*in[j+4][i-1] + c[7][4]*in[j+4][i+1] + c[7][5]*in[j+4][i+2] + c[7][6]*in[j+4][i+3] + c[7][7]*in[j+4][i+4];
				}
			}
		}
	}


	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 4; j < N-4; j+=4) {
				for (i = 4; i < N-4; i+=4) {
#pragma begin stencil1 unroll j=4,i=1 print-intrinsics true acc-size 1
					out[j][i] = c[0][0]*in[j-4][i-4] +  c[0][1]*in[j-4][i-3] + c[0][2]*in[j-4][i-2] + c[0][3]*in[j-4][i-1] + c[0][4]*in[j-4][i+1] + c[0][5]*in[j-4][i+2] + c[0][6]*in[j-4][i+3] + c[0][7]*in[j-4][i+4] +
						c[1][0]*in[j-3][i-4] +  c[1][1]*in[j-3][i-3] + c[1][2]*in[j-3][i-2] + c[1][3]*in[j-3][i-1] + c[1][4]*in[j-3][i+1] + c[1][5]*in[j-3][i+2] + c[1][6]*in[j-3][i+3] + c[1][7]*in[j-3][i+4] +
						c[2][0]*in[j-2][i-4] +  c[2][1]*in[j-2][i-3] + c[2][2]*in[j-2][i-2] + c[2][3]*in[j-2][i-1] + c[2][4]*in[j-2][i+1] + c[2][5]*in[j-2][i+2] + c[2][6]*in[j-2][i+3] + c[2][7]*in[j-2][i+4] +
						c[3][0]*in[j-1][i-4] +  c[3][1]*in[j-1][i-3] + c[3][2]*in[j-1][i-2] + c[3][3]*in[j-1][i-1] + c[3][4]*in[j-1][i+1] + c[3][5]*in[j-1][i+2] + c[3][6]*in[j-1][i+3] + c[3][7]*in[j-1][i+4] +
						c[4][0]*in[j+1][i-4] +  c[4][1]*in[j+1][i-3] + c[4][2]*in[j+1][i-2] + c[4][3]*in[j+1][i-1] + c[4][4]*in[j+1][i+1] + c[4][5]*in[j+1][i+2] + c[4][6]*in[j+1][i+3] + c[4][7]*in[j+1][i+4] +
						c[5][0]*in[j+2][i-4] +  c[5][1]*in[j+2][i-3] + c[5][2]*in[j+2][i-2] + c[5][3]*in[j+2][i-1] + c[5][4]*in[j+2][i+1] + c[5][5]*in[j+2][i+2] + c[5][6]*in[j+2][i+3] + c[5][7]*in[j+2][i+4] +
						c[6][0]*in[j+3][i-4] +  c[6][1]*in[j+3][i-3] + c[6][2]*in[j+3][i-2] + c[6][3]*in[j+3][i-1] + c[6][4]*in[j+3][i+1] + c[6][5]*in[j+3][i+2] + c[6][6]*in[j+3][i+3] + c[6][7]*in[j+3][i+4] +
						c[7][0]*in[j+4][i-4] +  c[7][1]*in[j+4][i-3] + c[7][2]*in[j+4][i-2] + c[7][3]*in[j+4][i-1] + c[7][4]*in[j+4][i+1] + c[7][5]*in[j+4][i+2] + c[7][6]*in[j+4][i+3] + c[7][7]*in[j+4][i+4];
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %.6lf\n",  (double)8192*8192*127*5/(end_time - start_time)/1e9);
}
