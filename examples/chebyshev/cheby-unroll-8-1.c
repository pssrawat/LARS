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
			for (k = 1; k < N-1; k+=8) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						out[k][j][i] = Ac[k][j][i] + c1 * (Ac[k][j][i] + Ap[k][j][i]) + c2 * Dinv[k][j][i] * (RHS[k][j][i] + (Ac[k][j][i] + h2inv * (0.03f * (Ac[k-1][j-1][i-1] + Ac[k-1][j-1][i+1] + Ac[k-1][j+1][i-1] + Ac[k-1][j+1][i+1] + Ac[k+1][j-1][i-1] + Ac[k+1][j-1][i+1] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1]) + 0.1f * (Ac[k-1][j-1][i] + Ac[k-1][j][i-1] + Ac[k-1][j][i+1] + Ac[k-1][j+1][i] + Ac[k][j-1][i-1] + Ac[k][j-1][i+1] + Ac[k][j+1][i-1] + Ac[k][j+1][i+1] + Ac[k+1][j-1][i] + Ac[k+1][j][i-1] + Ac[k+1][j][i+1] + Ac[k+1][j+1][i]) + 0.46f * (Ac[k-1][j][i] + Ac[k][j-1][i] + Ac[k][j][i-1] + Ac[k][j][i+1] + Ac[k][j+1][i] + Ac[k+1][j][i]) + 4.26f * Ac[k][j][i])));
						out[k+1][j][i] = Ac[k+1][j][i] + c1 * (Ac[k+1][j][i] + Ap[k+1][j][i]) + c2 * Dinv[k+1][j][i] * (RHS[k+1][j][i] + (Ac[k+1][j][i] + h2inv * (0.03f * (Ac[k+1-1][j-1][i-1] + Ac[k+1-1][j-1][i+1] + Ac[k+1-1][j+1][i-1] + Ac[k+1-1][j+1][i+1] + Ac[k+1+1][j-1][i-1] + Ac[k+1+1][j-1][i+1] + Ac[k+1+1][j+1][i-1] + Ac[k+1+1][j+1][i+1]) + 0.1f * (Ac[k+1-1][j-1][i] + Ac[k+1-1][j][i-1] + Ac[k+1-1][j][i+1] + Ac[k+1-1][j+1][i] + Ac[k+1][j-1][i-1] + Ac[k+1][j-1][i+1] + Ac[k+1][j+1][i-1] + Ac[k+1][j+1][i+1] + Ac[k+1+1][j-1][i] + Ac[k+1+1][j][i-1] + Ac[k+1+1][j][i+1] + Ac[k+1+1][j+1][i]) + 0.46f * (Ac[k+1-1][j][i] + Ac[k+1][j-1][i] + Ac[k+1][j][i-1] + Ac[k+1][j][i+1] + Ac[k+1][j+1][i] + Ac[k+1+1][j][i]) + 4.26f * Ac[k+1][j][i])));
						out[k+2][j][i] = Ac[k+2][j][i] + c1 * (Ac[k+2][j][i] + Ap[k+2][j][i]) + c2 * Dinv[k+2][j][i] * (RHS[k+2][j][i] + (Ac[k+2][j][i] + h2inv * (0.03f * (Ac[k+2-1][j-1][i-1] + Ac[k+2-1][j-1][i+1] + Ac[k+2-1][j+1][i-1] + Ac[k+2-1][j+1][i+1] + Ac[k+2+1][j-1][i-1] + Ac[k+2+1][j-1][i+1] + Ac[k+2+1][j+1][i-1] + Ac[k+2+1][j+1][i+1]) + 0.1f * (Ac[k+2-1][j-1][i] + Ac[k+2-1][j][i-1] + Ac[k+2-1][j][i+1] + Ac[k+2-1][j+1][i] + Ac[k+2][j-1][i-1] + Ac[k+2][j-1][i+1] + Ac[k+2][j+1][i-1] + Ac[k+2][j+1][i+1] + Ac[k+2+1][j-1][i] + Ac[k+2+1][j][i-1] + Ac[k+2+1][j][i+1] + Ac[k+2+1][j+1][i]) + 0.46f * (Ac[k+2-1][j][i] + Ac[k+2][j-1][i] + Ac[k+2][j][i-1] + Ac[k+2][j][i+1] + Ac[k+2][j+1][i] + Ac[k+2+1][j][i]) + 4.26f * Ac[k+2][j][i])));
						out[k+3][j][i] = Ac[k+3][j][i] + c1 * (Ac[k+3][j][i] + Ap[k+3][j][i]) + c2 * Dinv[k+3][j][i] * (RHS[k+3][j][i] + (Ac[k+3][j][i] + h2inv * (0.03f * (Ac[k+3-1][j-1][i-1] + Ac[k+3-1][j-1][i+1] + Ac[k+3-1][j+1][i-1] + Ac[k+3-1][j+1][i+1] + Ac[k+3+1][j-1][i-1] + Ac[k+3+1][j-1][i+1] + Ac[k+3+1][j+1][i-1] + Ac[k+3+1][j+1][i+1]) + 0.1f * (Ac[k+3-1][j-1][i] + Ac[k+3-1][j][i-1] + Ac[k+3-1][j][i+1] + Ac[k+3-1][j+1][i] + Ac[k+3][j-1][i-1] + Ac[k+3][j-1][i+1] + Ac[k+3][j+1][i-1] + Ac[k+3][j+1][i+1] + Ac[k+3+1][j-1][i] + Ac[k+3+1][j][i-1] + Ac[k+3+1][j][i+1] + Ac[k+3+1][j+1][i]) + 0.46f * (Ac[k+3-1][j][i] + Ac[k+3][j-1][i] + Ac[k+3][j][i-1] + Ac[k+3][j][i+1] + Ac[k+3][j+1][i] + Ac[k+3+1][j][i]) + 4.26f * Ac[k+3][j][i])));
						out[k+4][j][i] = Ac[k+4][j][i] + c1 * (Ac[k+4][j][i] + Ap[k+4][j][i]) + c2 * Dinv[k+4][j][i] * (RHS[k+4][j][i] + (Ac[k+4][j][i] + h2inv * (0.03f * (Ac[k+4-1][j-1][i-1] + Ac[k+4-1][j-1][i+1] + Ac[k+4-1][j+1][i-1] + Ac[k+4-1][j+1][i+1] + Ac[k+4+1][j-1][i-1] + Ac[k+4+1][j-1][i+1] + Ac[k+4+1][j+1][i-1] + Ac[k+4+1][j+1][i+1]) + 0.1f * (Ac[k+4-1][j-1][i] + Ac[k+4-1][j][i-1] + Ac[k+4-1][j][i+1] + Ac[k+4-1][j+1][i] + Ac[k+4][j-1][i-1] + Ac[k+4][j-1][i+1] + Ac[k+4][j+1][i-1] + Ac[k+4][j+1][i+1] + Ac[k+4+1][j-1][i] + Ac[k+4+1][j][i-1] + Ac[k+4+1][j][i+1] + Ac[k+4+1][j+1][i]) + 0.46f * (Ac[k+4-1][j][i] + Ac[k+4][j-1][i] + Ac[k+4][j][i-1] + Ac[k+4][j][i+1] + Ac[k+4][j+1][i] + Ac[k+4+1][j][i]) + 4.26f * Ac[k+4][j][i])));
						out[k+5][j][i] = Ac[k+5][j][i] + c1 * (Ac[k+5][j][i] + Ap[k+5][j][i]) + c2 * Dinv[k+5][j][i] * (RHS[k+5][j][i] + (Ac[k+5][j][i] + h2inv * (0.03f * (Ac[k+5-1][j-1][i-1] + Ac[k+5-1][j-1][i+1] + Ac[k+5-1][j+1][i-1] + Ac[k+5-1][j+1][i+1] + Ac[k+5+1][j-1][i-1] + Ac[k+5+1][j-1][i+1] + Ac[k+5+1][j+1][i-1] + Ac[k+5+1][j+1][i+1]) + 0.1f * (Ac[k+5-1][j-1][i] + Ac[k+5-1][j][i-1] + Ac[k+5-1][j][i+1] + Ac[k+5-1][j+1][i] + Ac[k+5][j-1][i-1] + Ac[k+5][j-1][i+1] + Ac[k+5][j+1][i-1] + Ac[k+5][j+1][i+1] + Ac[k+5+1][j-1][i] + Ac[k+5+1][j][i-1] + Ac[k+5+1][j][i+1] + Ac[k+5+1][j+1][i]) + 0.46f * (Ac[k+5-1][j][i] + Ac[k+5][j-1][i] + Ac[k+5][j][i-1] + Ac[k+5][j][i+1] + Ac[k+5][j+1][i] + Ac[k+5+1][j][i]) + 4.26f * Ac[k+5][j][i])));
						out[k+6][j][i] = Ac[k+6][j][i] + c1 * (Ac[k+6][j][i] + Ap[k+6][j][i]) + c2 * Dinv[k+6][j][i] * (RHS[k+6][j][i] + (Ac[k+6][j][i] + h2inv * (0.03f * (Ac[k+6-1][j-1][i-1] + Ac[k+6-1][j-1][i+1] + Ac[k+6-1][j+1][i-1] + Ac[k+6-1][j+1][i+1] + Ac[k+6+1][j-1][i-1] + Ac[k+6+1][j-1][i+1] + Ac[k+6+1][j+1][i-1] + Ac[k+6+1][j+1][i+1]) + 0.1f * (Ac[k+6-1][j-1][i] + Ac[k+6-1][j][i-1] + Ac[k+6-1][j][i+1] + Ac[k+6-1][j+1][i] + Ac[k+6][j-1][i-1] + Ac[k+6][j-1][i+1] + Ac[k+6][j+1][i-1] + Ac[k+6][j+1][i+1] + Ac[k+6+1][j-1][i] + Ac[k+6+1][j][i-1] + Ac[k+6+1][j][i+1] + Ac[k+6+1][j+1][i]) + 0.46f * (Ac[k+6-1][j][i] + Ac[k+6][j-1][i] + Ac[k+6][j][i-1] + Ac[k+6][j][i+1] + Ac[k+6][j+1][i] + Ac[k+6+1][j][i]) + 4.26f * Ac[k+6][j][i])));
						out[k+7][j][i] = Ac[k+7][j][i] + c1 * (Ac[k+7][j][i] + Ap[k+7][j][i]) + c2 * Dinv[k+7][j][i] * (RHS[k+7][j][i] + (Ac[k+7][j][i] + h2inv * (0.03f * (Ac[k+7-1][j-1][i-1] + Ac[k+7-1][j-1][i+1] + Ac[k+7-1][j+1][i-1] + Ac[k+7-1][j+1][i+1] + Ac[k+7+1][j-1][i-1] + Ac[k+7+1][j-1][i+1] + Ac[k+7+1][j+1][i-1] + Ac[k+7+1][j+1][i+1]) + 0.1f * (Ac[k+7-1][j-1][i] + Ac[k+7-1][j][i-1] + Ac[k+7-1][j][i+1] + Ac[k+7-1][j+1][i] + Ac[k+7][j-1][i-1] + Ac[k+7][j-1][i+1] + Ac[k+7][j+1][i-1] + Ac[k+7][j+1][i+1] + Ac[k+7+1][j-1][i] + Ac[k+7+1][j][i-1] + Ac[k+7+1][j][i+1] + Ac[k+7+1][j+1][i]) + 0.46f * (Ac[k+7-1][j][i] + Ac[k+7][j-1][i] + Ac[k+7][j][i-1] + Ac[k+7][j][i+1] + Ac[k+7][j+1][i] + Ac[k+7+1][j][i]) + 4.26f * Ac[k+7][j][i])));
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("unroll: %6lf\n", (double)512*512*512*39*5/(end_time - start_time)/1e9);
}
