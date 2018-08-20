#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void cheby_opt (double * __restrict__ t_Ac, double * __restrict__ t_Ap, double * __restrict__ t_Dinv, double * __restrict__ t_RHS, double * __restrict__ t_out, double c1, double c2, double h2inv, int N) {
	double (*Ac)[514][514] = (double (*)[514][514]) t_Ac;
	double (*Ap)[514][514] = (double (*)[514][514]) t_Ap;
	double (*Dinv)[514][514] = (double (*)[514][514]) t_Dinv;
	double (*RHS)[514][514] = (double (*)[514][514]) t_RHS;
	double (*out)[514][514] = (double (*)[514][514]) t_out;

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
						out[k][j][i] = 0.1 * Ac[k][j][i] + c1 * (Ac[k][j][i] + Ap[k][j][i]) + c2 * Dinv[k][j][i] * (RHS[k][j][i] + (Ac[k][j][i] + h2inv * (0.03 * (Ac[k-1][j-1][i-1] + Ac[k-1][j-1][i+1] + Ac[k-1][j+1][i-1] + Ac[k-1][j+1][i+1] + Ac[k+1][j-1][i-1] + Ac[k+1][j-1][i+1] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1]) + 0.1 * (Ac[k-1][j-1][i] + Ac[k-1][j][i-1] + Ac[k-1][j][i+1] + Ac[k-1][j+1][i] + Ac[k][j-1][i-1] + Ac[k][j-1][i+1] + Ac[k][j+1][i-1] + Ac[k][j+1][i+1] + Ac[k+1][j-1][i] + Ac[k+1][j][i-1] + Ac[k+1][j][i+1] + Ac[k+1][j+1][i]) + 0.46 * (Ac[k-1][j][i] + Ac[k][j-1][i] + Ac[k][j][i-1] + Ac[k][j][i+1] + Ac[k][j+1][i] + Ac[k+1][j][i]) + 4.26 * Ac[k][j][i])));
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
					for (i = 1; i < N-1; i+=4) {
#pragma begin stencil1 unroll k=2,j=2,i=1 print-intrinsics true acc-size 1
					out[k][j][i] = Ac[k][j][i] + c1 * (Ac[k][j][i] + Ap[k][j][i]) + c2 * Dinv[k][j][i] * (RHS[k][j][i] + (Ac[k][j][i] + h2inv * (0.03 * (Ac[k-1][j-1][i-1] + Ac[k-1][j-1][i+1] + Ac[k-1][j+1][i-1] + Ac[k-1][j+1][i+1] + Ac[k+1][j-1][i-1] + Ac[k+1][j-1][i+1] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1]) + 0.1 * (Ac[k-1][j-1][i] + Ac[k-1][j][i-1] + Ac[k-1][j][i+1] + Ac[k-1][j+1][i] + Ac[k][j-1][i-1] + Ac[k][j-1][i+1] + Ac[k][j+1][i-1] + Ac[k][j+1][i+1] + Ac[k+1][j-1][i] + Ac[k+1][j][i-1] + Ac[k+1][j][i+1] + Ac[k+1][j+1][i]) + 0.46 * (Ac[k-1][j][i] + Ac[k][j-1][i] + Ac[k][j][i-1] + Ac[k][j][i+1] + Ac[k][j+1][i] + Ac[k+1][j][i]) + 4.26 * Ac[k][j][i])));
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)512*512*512*39*5/(end_time - start_time)/1e9);
}
