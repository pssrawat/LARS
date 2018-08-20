#include<stdio.h>
#include<stdlib.h>
#define EIGENVALUE 2.175

extern double rtclock (void);

void cheby_opt (double * __restrict__ t_Ac, double * __restrict__ t_Ap, double * __restrict__ t_Dinv, double * __restrict__ t_RHS, double * __restrict__ t_out, double c1, double c2, double h2inv, int N) {
	double (*Ac)[514][514] = (double (*)[514][514]) t_Ac;
	double (*Ap)[514][514] = (double (*)[514][514]) t_Ap;
	double (*Dinv)[514][514] = (double (*)[514][514]) t_Dinv;
	double (*RHS)[514][514] = (double (*)[514][514]) t_RHS;
	double (*out)[514][514] = (double (*)[514][514]) t_out;
	int t, i, j, k;                         
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
			for (k = 1; k < N-1; k+=2) {
				for (j = 1; j < N-1; j+=4) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						out[k][j][i] = Ac[k][j][i] + c1 * (Ac[k][j][i] + Ap[k][j][i]) + c2 * Dinv[k][j][i] * (RHS[k][j][i] + (Ac[k][j][i] + h2inv * (0.03 * (Ac[k-1][j-1][i-1] + Ac[k-1][j-1][i+1] + Ac[k-1][j+1][i-1] + Ac[k-1][j+1][i+1] + Ac[k+1][j-1][i-1] + Ac[k+1][j-1][i+1] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1]) + 0.1 * (Ac[k-1][j-1][i] + Ac[k-1][j][i-1] + Ac[k-1][j][i+1] + Ac[k-1][j+1][i] + Ac[k][j-1][i-1] + Ac[k][j-1][i+1] + Ac[k][j+1][i-1] + Ac[k][j+1][i+1] + Ac[k+1][j-1][i] + Ac[k+1][j][i-1] + Ac[k+1][j][i+1] + Ac[k+1][j+1][i]) + 0.46 * (Ac[k-1][j][i] + Ac[k][j-1][i] + Ac[k][j][i-1] + Ac[k][j][i+1] + Ac[k][j+1][i] + Ac[k+1][j][i]) + 4.26 * Ac[k][j][i])));
						out[k][j+1][i] = Ac[k][j+1][i] + c1 * (Ac[k][j+1][i] + Ap[k][j+1][i]) + c2 * Dinv[k][j+1][i] * (RHS[k][j+1][i] + (Ac[k][j+1][i] + h2inv * (0.03 * (Ac[k-1][j+1-1][i-1] + Ac[k-1][j+1-1][i+1] + Ac[k-1][j+1+1][i-1] + Ac[k-1][j+1+1][i+1] + Ac[k+1][j+1-1][i-1] + Ac[k+1][j+1-1][i+1] + Ac[k+1][j+1+1][i-1] + Ac[k+1][j+1+1][i+1]) + 0.1 * (Ac[k-1][j+1-1][i] + Ac[k-1][j+1][i-1] + Ac[k-1][j+1][i+1] + Ac[k-1][j+1+1][i] + Ac[k][j+1-1][i-1] + Ac[k][j+1-1][i+1] + Ac[k][j+1+1][i-1] + Ac[k][j+1+1][i+1] + Ac[k+1][j+1-1][i] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1] + Ac[k+1][j+1+1][i]) + 0.46 * (Ac[k-1][j+1][i] + Ac[k][j+1-1][i] + Ac[k][j+1][i-1] + Ac[k][j+1][i+1] + Ac[k][j+1+1][i] + Ac[k+1][j+1][i]) + 4.26 * Ac[k][j+1][i])));
						out[k][j+2][i] = Ac[k][j+2][i] + c1 * (Ac[k][j+2][i] + Ap[k][j+2][i]) + c2 * Dinv[k][j+2][i] * (RHS[k][j+2][i] + (Ac[k][j+2][i] + h2inv * (0.03 * (Ac[k-1][j+2-1][i-1] + Ac[k-1][j+2-1][i+1] + Ac[k-1][j+2+1][i-1] + Ac[k-1][j+2+1][i+1] + Ac[k+1][j+2-1][i-1] + Ac[k+1][j+2-1][i+1] + Ac[k+1][j+2+1][i-1] + Ac[k+1][j+2+1][i+1]) + 0.1 * (Ac[k-1][j+2-1][i] + Ac[k-1][j+2][i-1] + Ac[k-1][j+2][i+1] + Ac[k-1][j+2+1][i] + Ac[k][j+2-1][i-1] + Ac[k][j+2-1][i+1] + Ac[k][j+2+1][i-1] + Ac[k][j+2+1][i+1] + Ac[k+1][j+2-1][i] + Ac[k+1][j+2][i-1] + Ac[k+1][j+2][i+1] + Ac[k+1][j+2+1][i]) + 0.46 * (Ac[k-1][j+2][i] + Ac[k][j+2-1][i] + Ac[k][j+2][i-1] + Ac[k][j+2][i+1] + Ac[k][j+2+1][i] + Ac[k+1][j+2][i]) + 4.26 * Ac[k][j+2][i])));
						out[k][j+3][i] = Ac[k][j+3][i] + c1 * (Ac[k][j+3][i] + Ap[k][j+3][i]) + c2 * Dinv[k][j+3][i] * (RHS[k][j+3][i] + (Ac[k][j+3][i] + h2inv * (0.03 * (Ac[k-1][j+3-1][i-1] + Ac[k-1][j+3-1][i+1] + Ac[k-1][j+3+1][i-1] + Ac[k-1][j+3+1][i+1] + Ac[k+1][j+3-1][i-1] + Ac[k+1][j+3-1][i+1] + Ac[k+1][j+3+1][i-1] + Ac[k+1][j+3+1][i+1]) + 0.1 * (Ac[k-1][j+3-1][i] + Ac[k-1][j+3][i-1] + Ac[k-1][j+3][i+1] + Ac[k-1][j+3+1][i] + Ac[k][j+3-1][i-1] + Ac[k][j+3-1][i+1] + Ac[k][j+3+1][i-1] + Ac[k][j+3+1][i+1] + Ac[k+1][j+3-1][i] + Ac[k+1][j+3][i-1] + Ac[k+1][j+3][i+1] + Ac[k+1][j+3+1][i]) + 0.46 * (Ac[k-1][j+3][i] + Ac[k][j+3-1][i] + Ac[k][j+3][i-1] + Ac[k][j+3][i+1] + Ac[k][j+3+1][i] + Ac[k+1][j+3][i]) + 4.26 * Ac[k][j+3][i])));

						out[k+1][j][i] = Ac[k+1][j][i] + c1 * (Ac[k+1][j][i] + Ap[k+1][j][i]) + c2 * Dinv[k+1][j][i] * (RHS[k+1][j][i] + (Ac[k+1][j][i] + h2inv * (0.03 * (Ac[k+1-1][j-1][i-1] + Ac[k+1-1][j-1][i+1] + Ac[k+1-1][j+1][i-1] + Ac[k+1-1][j+1][i+1] + Ac[k+1+1][j-1][i-1] + Ac[k+1+1][j-1][i+1] + Ac[k+1+1][j+1][i-1] + Ac[k+1+1][j+1][i+1]) + 0.1 * (Ac[k+1-1][j-1][i] + Ac[k+1-1][j][i-1] + Ac[k+1-1][j][i+1] + Ac[k+1-1][j+1][i] + Ac[k+1][j-1][i-1] + Ac[k+1][j-1][i+1] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1] + Ac[k+1+1][j-1][i] + Ac[k+1+1][j][i-1] + Ac[k+1+1][j][i+1] + Ac[k+1+1][j+1][i]) + 0.46 * (Ac[k+1-1][j][i] + Ac[k+1][j-1][i] + Ac[k+1][j][i-1] + Ac[k+1][j][i+1] + Ac[k+1][j+1][i] + Ac[k+1+1][j][i]) + 4.26 * Ac[k+1][j][i])));
						out[k+1][j+1][i] = Ac[k+1][j+1][i] + c1 * (Ac[k+1][j+1][i] + Ap[k+1][j+1][i]) + c2 * Dinv[k+1][j+1][i] * (RHS[k+1][j+1][i] + (Ac[k+1][j+1][i] + h2inv * (0.03 * (Ac[k+1-1][j+1-1][i-1] + Ac[k+1-1][j+1-1][i+1] + Ac[k+1-1][j+1+1][i-1] + Ac[k+1-1][j+1+1][i+1] + Ac[k+1+1][j+1-1][i-1] + Ac[k+1+1][j+1-1][i+1] + Ac[k+1+1][j+1+1][i-1] + Ac[k+1+1][j+1+1][i+1]) + 0.1 * (Ac[k+1-1][j+1-1][i] + Ac[k+1-1][j+1][i-1] + Ac[k+1-1][j+1][i+1] + Ac[k+1-1][j+1+1][i] + Ac[k+1][j+1-1][i-1] + Ac[k+1][j+1-1][i+1] + Ac[k+1][j+1+1][i-1] + Ac[k+1][j+1+1][i+1] + Ac[k+1+1][j+1-1][i] + Ac[k+1+1][j+1][i-1] + Ac[k+1+1][j+1][i+1] + Ac[k+1+1][j+1+1][i]) + 0.46 * (Ac[k+1-1][j+1][i] + Ac[k+1][j+1-1][i] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1] + Ac[k+1][j+1+1][i] + Ac[k+1+1][j+1][i]) + 4.26 * Ac[k+1][j+1][i])));
						out[k+1][j+2][i] = Ac[k+1][j+2][i] + c1 * (Ac[k+1][j+2][i] + Ap[k+1][j+2][i]) + c2 * Dinv[k+1][j+2][i] * (RHS[k+1][j+2][i] + (Ac[k+1][j+2][i] + h2inv * (0.03 * (Ac[k+1-1][j+2-1][i-1] + Ac[k+1-1][j+2-1][i+1] + Ac[k+1-1][j+2+1][i-1] + Ac[k+1-1][j+2+1][i+1] + Ac[k+1+1][j+2-1][i-1] + Ac[k+1+1][j+2-1][i+1] + Ac[k+1+1][j+2+1][i-1] + Ac[k+1+1][j+2+1][i+1]) + 0.1 * (Ac[k+1-1][j+2-1][i] + Ac[k+1-1][j+2][i-1] + Ac[k+1-1][j+2][i+1] + Ac[k+1-1][j+2+1][i] + Ac[k+1][j+2-1][i-1] + Ac[k+1][j+2-1][i+1] + Ac[k+1][j+2+1][i-1] + Ac[k+1][j+2+1][i+1] + Ac[k+1+1][j+2-1][i] + Ac[k+1+1][j+2][i-1] + Ac[k+1+1][j+2][i+1] + Ac[k+1+1][j+2+1][i]) + 0.46 * (Ac[k+1-1][j+2][i] + Ac[k+1][j+2-1][i] + Ac[k+1][j+2][i-1] + Ac[k+1][j+2][i+1] + Ac[k+1][j+2+1][i] + Ac[k+1+1][j+2][i]) + 4.26 * Ac[k+1][j+2][i])));
						out[k+1][j+3][i] = Ac[k+1][j+3][i] + c1 * (Ac[k+1][j+3][i] + Ap[k+1][j+3][i]) + c2 * Dinv[k+1][j+3][i] * (RHS[k+1][j+3][i] + (Ac[k+1][j+3][i] + h2inv * (0.03 * (Ac[k+1-1][j+3-1][i-1] + Ac[k+1-1][j+3-1][i+1] + Ac[k+1-1][j+3+1][i-1] + Ac[k+1-1][j+3+1][i+1] + Ac[k+1+1][j+3-1][i-1] + Ac[k+1+1][j+3-1][i+1] + Ac[k+1+1][j+3+1][i-1] + Ac[k+1+1][j+3+1][i+1]) + 0.1 * (Ac[k+1-1][j+3-1][i] + Ac[k+1-1][j+3][i-1] + Ac[k+1-1][j+3][i+1] + Ac[k+1-1][j+3+1][i] + Ac[k+1][j+3-1][i-1] + Ac[k+1][j+3-1][i+1] + Ac[k+1][j+3+1][i-1] + Ac[k+1][j+3+1][i+1] + Ac[k+1+1][j+3-1][i] + Ac[k+1+1][j+3][i-1] + Ac[k+1+1][j+3][i+1] + Ac[k+1+1][j+3+1][i]) + 0.46 * (Ac[k+1-1][j+3][i] + Ac[k+1][j+3-1][i] + Ac[k+1][j+3][i-1] + Ac[k+1][j+3][i+1] + Ac[k+1][j+3+1][i] + Ac[k+1+1][j+3][i]) + 4.26 * Ac[k+1][j+3][i])));
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("unroll: %6lf\n", (double)512*512*512*39*5/(end_time - start_time)/1e9);
}
