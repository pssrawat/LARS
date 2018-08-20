#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void poisson_opt (double*, double*, int);

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
	int N = 514;
	double (*in)[514][514] = (double (*)[514][514]) getRandom3DArray (514, 514, 514);
	double (*out_ref)[514][514] = (double (*)[514][514]) getZero3DArray (514, 514, 514);
	double (*out)[514][514] = (double (*)[514][514]) getZero3DArray (514, 514, 514);

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
						out_ref[k][j][i] = 0.1 * 2.666*in[k][j][i] -
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
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma clang loop vectorize (enable) interleave(enable)
#pragma GCC ivdep
					for (i = 1; i < N-1; i++) {
						out_ref[k][j][i] = 2.666*in[k][j][i] -
							(0.166*in[k-1][j][i] + 0.166*in[k+1][j][i] + 0.166*in[k][j-1][i] + 0.166*in[k][j+1][i] + 0.166*in[k][j][i+1] + 0.166*in[k][j][i-1])-
							(0.0833*in[k-1][j-1][i] + 0.0833*in[k+1][j-1][i] + 0.0833*in[k-1][j+1][i] + 0.0833*in[k+1][j+1][i] +
							 0.0833*in[k-1][j][i-1] + 0.0833*in[k+1][j][i-1] + 0.0833*in[k][j-1][i-1] + 0.0833*in[k][j+1][i-1] +
							 0.0833*in[k-1][j][i+1] + 0.0833*in[k+1][j][i+1] + 0.0833*in[k][j-1][i+1] + 0.0833*in[k][j+1][i+1]);
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)512*512*512*21*5/(end_time  - start_time)/1e9);

	poisson_opt ((double*)out, (double*)in, N);

	double error = checkError3D (N, N, 0, (double*)out, (double*)out_ref, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
}
