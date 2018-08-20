#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void j2d81pt_opt (double *t_in, double *t_out, double *t_c, int N) {
	double (*in)[8200] = (double (*)[8200]) t_in;
	double (*out)[8200] = (double (*)[8200]) t_out;
	double (*c)[9] = (double (*)[9]) t_c;

	int t, i, j, i0, j0;
	double start_time, end_time;

	// Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 0; j < N-8; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 0; i < N-8; i++) {
					out[j][i] = 0.1 * 
						c[0][0]*in[j][i] +   c[0][1]*in[j][i+1] +   c[0][2]*in[j][i+2] +   c[0][3]*in[j][i+3] +   c[0][4]*in[j][i+4] +   c[0][5]*in[j][i+5] +   c[0][6]*in[j][i+6] +   c[0][7]*in[j][i+7] +   c[0][8]*in[j][i+8] +
						c[1][0]*in[j+1][i] + c[1][1]*in[j+1][i+1] + c[1][2]*in[j+1][i+2] + c[1][3]*in[j+1][i+3] + c[1][4]*in[j+1][i+4] + c[1][5]*in[j+1][i+5] + c[1][6]*in[j+1][i+6] + c[1][7]*in[j+1][i+7] + c[1][8]*in[j+1][i+8] +
						c[2][0]*in[j+2][i] + c[2][1]*in[j+2][i+1] + c[2][2]*in[j+2][i+2] + c[2][3]*in[j+2][i+3] + c[2][4]*in[j+2][i+4] + c[2][5]*in[j+2][i+5] + c[2][6]*in[j+2][i+6] + c[2][7]*in[j+2][i+7] + c[2][8]*in[j+2][i+8] +
						c[3][0]*in[j+3][i] + c[3][1]*in[j+3][i+1] + c[3][2]*in[j+3][i+2] + c[3][3]*in[j+3][i+3] + c[3][4]*in[j+3][i+4] + c[3][5]*in[j+3][i+5] + c[3][6]*in[j+3][i+6] + c[3][7]*in[j+3][i+7] + c[3][8]*in[j+3][i+8] +
						c[4][0]*in[j+4][i] + c[4][1]*in[j+4][i+1] + c[4][2]*in[j+4][i+2] + c[4][3]*in[j+4][i+3] + c[4][4]*in[j+4][i+4] + c[4][5]*in[j+4][i+5] + c[4][6]*in[j+4][i+6] + c[4][7]*in[j+4][i+7] + c[4][8]*in[j+4][i+8] +
						c[5][0]*in[j+5][i] + c[5][1]*in[j+5][i+1] + c[5][2]*in[j+5][i+2] + c[5][3]*in[j+5][i+3] + c[5][4]*in[j+5][i+4] + c[5][5]*in[j+5][i+5] + c[5][6]*in[j+5][i+6] + c[5][7]*in[j+5][i+7] + c[5][8]*in[j+5][i+8] +
						c[6][0]*in[j+6][i] + c[6][1]*in[j+6][i+1] + c[6][2]*in[j+6][i+2] + c[6][3]*in[j+6][i+3] + c[6][4]*in[j+6][i+4] + c[6][5]*in[j+6][i+5] + c[6][6]*in[j+6][i+6] + c[6][7]*in[j+6][i+7] + c[6][8]*in[j+6][i+8] +
						c[7][0]*in[j+7][i] + c[7][1]*in[j+7][i+1] + c[7][2]*in[j+7][i+2] + c[7][3]*in[j+7][i+3] + c[7][4]*in[j+7][i+4] + c[7][5]*in[j+7][i+5] + c[7][6]*in[j+7][i+6] + c[7][7]*in[j+7][i+7] + c[7][8]*in[j+7][i+8] +
						c[8][0]*in[j+8][i] + c[8][1]*in[j+8][i+1] + c[8][2]*in[j+8][i+2] + c[8][3]*in[j+8][i+3] + c[8][4]*in[j+8][i+4] + c[8][5]*in[j+8][i+5] + c[8][6]*in[j+8][i+6] + c[8][7]*in[j+8][i+7] + c[8][8]*in[j+8][i+8];
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 0; j < N-8; j+=2) {
				for (i = 0; i < N-8; i+=1) {
#pragma begin stencil1 unroll j=2,i=1 print-intrinsics false acc-size 1
					out[j][i] = 
						c[0][0]*in[j][i] +   c[0][1]*in[j][i+1] +   c[0][2]*in[j][i+2] +   c[0][3]*in[j][i+3] +   c[0][4]*in[j][i+4] +   c[0][5]*in[j][i+5] +   c[0][6]*in[j][i+6] +   c[0][7]*in[j][i+7] +   c[0][8]*in[j][i+8] +
						c[1][0]*in[j+1][i] + c[1][1]*in[j+1][i+1] + c[1][2]*in[j+1][i+2] + c[1][3]*in[j+1][i+3] + c[1][4]*in[j+1][i+4] + c[1][5]*in[j+1][i+5] + c[1][6]*in[j+1][i+6] + c[1][7]*in[j+1][i+7] + c[1][8]*in[j+1][i+8] +
						c[2][0]*in[j+2][i] + c[2][1]*in[j+2][i+1] + c[2][2]*in[j+2][i+2] + c[2][3]*in[j+2][i+3] + c[2][4]*in[j+2][i+4] + c[2][5]*in[j+2][i+5] + c[2][6]*in[j+2][i+6] + c[2][7]*in[j+2][i+7] + c[2][8]*in[j+2][i+8] +
						c[3][0]*in[j+3][i] + c[3][1]*in[j+3][i+1] + c[3][2]*in[j+3][i+2] + c[3][3]*in[j+3][i+3] + c[3][4]*in[j+3][i+4] + c[3][5]*in[j+3][i+5] + c[3][6]*in[j+3][i+6] + c[3][7]*in[j+3][i+7] + c[3][8]*in[j+3][i+8] +
						c[4][0]*in[j+4][i] + c[4][1]*in[j+4][i+1] + c[4][2]*in[j+4][i+2] + c[4][3]*in[j+4][i+3] + c[4][4]*in[j+4][i+4] + c[4][5]*in[j+4][i+5] + c[4][6]*in[j+4][i+6] + c[4][7]*in[j+4][i+7] + c[4][8]*in[j+4][i+8] +
						c[5][0]*in[j+5][i] + c[5][1]*in[j+5][i+1] + c[5][2]*in[j+5][i+2] + c[5][3]*in[j+5][i+3] + c[5][4]*in[j+5][i+4] + c[5][5]*in[j+5][i+5] + c[5][6]*in[j+5][i+6] + c[5][7]*in[j+5][i+7] + c[5][8]*in[j+5][i+8] +
						c[6][0]*in[j+6][i] + c[6][1]*in[j+6][i+1] + c[6][2]*in[j+6][i+2] + c[6][3]*in[j+6][i+3] + c[6][4]*in[j+6][i+4] + c[6][5]*in[j+6][i+5] + c[6][6]*in[j+6][i+6] + c[6][7]*in[j+6][i+7] + c[6][8]*in[j+6][i+8] +
						c[7][0]*in[j+7][i] + c[7][1]*in[j+7][i+1] + c[7][2]*in[j+7][i+2] + c[7][3]*in[j+7][i+3] + c[7][4]*in[j+7][i+4] + c[7][5]*in[j+7][i+5] + c[7][6]*in[j+7][i+6] + c[7][7]*in[j+7][i+7] + c[7][8]*in[j+7][i+8] +
						c[8][0]*in[j+8][i] + c[8][1]*in[j+8][i+1] + c[8][2]*in[j+8][i+2] + c[8][3]*in[j+8][i+3] + c[8][4]*in[j+8][i+4] + c[8][5]*in[j+8][i+5] + c[8][6]*in[j+8][i+6] + c[8][7]*in[j+8][i+7] + c[8][8]*in[j+8][i+8];
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %.6lf\n",  (double)8192*8192*161*5/(end_time - start_time)/1e9);
}
