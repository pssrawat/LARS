#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void helmholtz_opt (double *out, double *alpha, double *x, double *beta_i, double *beta_j, double *beta_k, int N);

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
	double (*alpha)[514+37][514+37] = (double (*)[514+37][514+37]) getRandom3DArray (514+37, 514+37, 514+37);
	double (*x)[514+37][514+37] = (double (*)[514+37][514+37]) getRandom3DArray (514+37, 514+37, 514+37);
	double (*beta_i)[514+37][514+37] = (double (*)[514+37][514+37]) getRandom3DArray (514+37, 514+37, 514+37);
	double (*beta_j)[514+37][514+37] = (double (*)[514+37][514+37]) getRandom3DArray (514+37, 514+37, 514+37);
	double (*beta_k)[514+37][514+37] = (double (*)[514+37][514+37]) getRandom3DArray (514+37, 514+37, 514+37);
	double (*out_ref)[514+37][514+37] = (double (*)[514+37][514+37]) getZero3DArray (514+37, 514+37, 514+37);
	double (*out)[514+37][514+37] = (double (*)[514+37][514+37]) getZero3DArray (514+37, 514+37, 514+37);
	double a = 1.141;
	double b = 0.121;
	double h2inv = 6.215;

	int i, j, k, t;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma clang loop vectorize (enable) interleave(enable)
#pragma GCC ivdep
					for (i = 1; i < N-1; i++) {
						out_ref[k][j][i] = 0.1 * a*alpha[k][j][i]*x[k][j][i] - b*h2inv*(
								beta_i[k][j][i+1] * (x[k][j][i+1] - x[k][j][i]) + 
								beta_i[k][j][i] * (x[k][j][i-1] - x[k][j][i]) + 
								beta_j[k][j+1][i] * (x[k][j+1][i] - x[k][j][i]) + 
								beta_j[k][j][i] * (x[k][j-1][i] - x[k][j][i]) + 
								beta_k[k+1][j][i] * (x[k+1][j][i] - x[k][j][i]) + 
								beta_k[k][j][i] * (x[k-1][j][i] - x[k][j][i]));
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
			for (k = 1; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma clang loop vectorize (enable) interleave(enable)
#pragma GCC ivdep
					for (i = 1; i < N-1; i++) {
						out_ref[k][j][i] = a*alpha[k][j][i]*x[k][j][i] - b*h2inv*(
								beta_i[k][j][i+1] * (x[k][j][i+1] - x[k][j][i]) + 
								beta_i[k][j][i] * (x[k][j][i-1] - x[k][j][i]) + 
								beta_j[k][j+1][i] * (x[k][j+1][i] - x[k][j][i]) + 
								beta_j[k][j][i] * (x[k][j-1][i] - x[k][j][i]) + 
								beta_k[k+1][j][i] * (x[k+1][j][i] - x[k][j][i]) + 
								beta_k[k][j][i] * (x[k-1][j][i] - x[k][j][i]));
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)512*512*512*22*5/(end_time  - start_time)/1e9);

	helmholtz_opt ((double*)out, (double*)alpha, (double*)x, (double*)beta_i, (double*)beta_j, (double*)beta_k, N);

	double error = checkError3D (N, N, 37, (double*)out, (double*)out_ref, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
	printf("error %e\n",error);
}
