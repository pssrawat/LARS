#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void j2d64pt_opt (double*, double*, double*, int);

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
	int N = 8200;
	double (*in)[8200] = (double (*)[8200]) getRandom2DArray (8200, 8200);
	double (*out_ref)[8200] = (double (*)[8200]) getZero2DArray (8200, 8200);
	double (*out_unroll)[8200] = (double (*)[8200]) getZero2DArray (8200, 8200);
	double (*out)[8200] = (double (*)[8200]) getZero2DArray (8200, 8200);
	double (*c)[8] = (double (*)[8]) getRandom2DArray (8, 8);

	int t, i, j;
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
					out_ref[j][i] =  0.1 * 
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
			for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 4; i < N-4; i++) {
					out_ref[j][i] =  c[0][0]*in[j-4][i-4] +  c[0][1]*in[j-4][i-3] + c[0][2]*in[j-4][i-2] + c[0][3]*in[j-4][i-1] + c[0][4]*in[j-4][i+1] + c[0][5]*in[j-4][i+2] + c[0][6]*in[j-4][i+3] + c[0][7]*in[j-4][i+4] +
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
	end_time = rtclock ();
        printf ("orig: %6lf\n", (double)8192*8192*127*5/(end_time - start_time)/1e9);

        j2d64pt_opt ((double*)in, (double*)out, (double*)c, N);

        double error = checkError2D (N, 0, (double*)out, (double*) out_ref, 4, N-4, 4, N-4);
        if (error > TOLERANCE)
                printf ("error %e\n", error);
}
