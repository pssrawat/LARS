#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void j2d25pt_opt (double*, double*, double*, int);

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
	int N = 8196;
	double (*in)[8196] = (double (*)[8196]) getRandom2DArray (8196, 8196);
	double (*out_ref)[8196] = (double (*)[8196]) getZero2DArray (8196, 8196);
	double (*out)[8196] = (double (*)[8196]) getZero2DArray (8196, 8196);
	double (*c)[5] = (double (*)[5]) getRandom2DArray (5, 5);

	int t, i, j;
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
					out_ref[j][i] =  c[0][0]*in[j-2][i-2] + c[0][1]*in[j-2][i-1] + c[0][2]*in[j-2][i] + c[0][3]*in[j-2][i+1] + c[0][4]*in[j-2][i+2] +
						c[1][0]*in[j-1][i-2] + c[1][1]*in[j-1][i-1] + c[1][2]*in[j-1][i] + c[1][3]*in[j-1][i+1] + c[1][4]*in[j-1][i+2] +
						c[2][0]*in[j][i-2] +   c[2][1]*in[j][i-1] +   c[2][2]*in[j][i] +   c[2][3]*in[j][i+1] +   c[2][4]*in[j][i+2] +
						c[3][0]*in[j+1][i-2] + c[3][1]*in[j+1][i-1] + c[3][2]*in[j+1][i] + c[3][3]*in[j+1][i+1] + c[3][4]*in[j+1][i+2] +
						c[4][0]*in[j+2][i-2] + c[4][1]*in[j+2][i-1] + c[4][2]*in[j+2][i] + c[4][3]*in[j+2][i+1] + c[4][4]*in[j+2][i+2];
				}
			}
		}
	}

	// TIME ORIG
	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 2; j < N-2; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 2; i < N-2; i++) {
					out_ref[j][i] =  c[0][0]*in[j-2][i-2] + c[0][1]*in[j-2][i-1] + c[0][2]*in[j-2][i] + c[0][3]*in[j-2][i+1] + c[0][4]*in[j-2][i+2] +
						c[1][0]*in[j-1][i-2] + c[1][1]*in[j-1][i-1] + c[1][2]*in[j-1][i] + c[1][3]*in[j-1][i+1] + c[1][4]*in[j-1][i+2] +
						c[2][0]*in[j][i-2] +   c[2][1]*in[j][i-1] +   c[2][2]*in[j][i] +   c[2][3]*in[j][i+1] +   c[2][4]*in[j][i+2] +
						c[3][0]*in[j+1][i-2] + c[3][1]*in[j+1][i-1] + c[3][2]*in[j+1][i] + c[3][3]*in[j+1][i+1] + c[3][4]*in[j+1][i+2] +
						c[4][0]*in[j+2][i-2] + c[4][1]*in[j+2][i-1] + c[4][2]*in[j+2][i] + c[4][3]*in[j+2][i+1] + c[4][4]*in[j+2][i+2];
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)8192*8192*49*5/(end_time - start_time)/1e9);

	j2d25pt_opt ((double*)in, (double*)out, (double*)c, N);

	double error = checkError2D (N, 0, (double*)out, (double*) out_ref, 2, N-2, 2, N-2);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
