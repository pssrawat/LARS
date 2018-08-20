#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void j2d49pt_opt (double*, double*, double*, int);

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
	int N = 8198;
	double (*in)[8198] = (double (*)[8198]) getRandom2DArray (8198, 8198);
	double (*out_ref)[8198] = (double (*)[8198]) getZero2DArray (8198, 8198);
	double (*out)[8198] = (double (*)[8198]) getZero2DArray (8198, 8198);
	double (*c)[7] = (double (*)[7]) getRandom2DArray (7, 7);

	int t, i, j;
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
					out_ref[j][i] = 0.1 * 
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
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 3; j < N-3; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable) 
				for (i = 3; i < N-3; i++) {
					out_ref[j][i] =
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
        end_time = rtclock (); 
        printf ("orig: %6lf\n", (double)8192*8192*97*5/(end_time - start_time)/1e9);
                
        j2d49pt_opt ((double*)in, (double*)out, (double*)c, N);

        double error = checkError2D (N, 0, (double*)out, (double*) out_ref, 3, N-3, 3, N-3);
        if (error > TOLERANCE)
                printf ("error %e\n", error);
}
