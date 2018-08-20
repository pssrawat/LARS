#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void j3d27pt_opt (double *t_in, double *t_out, double *c, int N) {
        double (*in)[514][514] = (double (*)[514][514]) t_in;
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
						out[k][j][i] = 0.1 * c[0] * in[k][j][i] + 
							c[1] * (in[k-1][j][i] + in[k+1][j][i] + in[k][j-1][i] + in[k][j+1][i] + in[k][j][i-1] + in[k][j][i+1]) + 
							c[2] * (in[k-1][j-1][i-1] + in[k-1][j-1][i+1] + in[k-1][j+1][i-1] + in[k-1][j+1][i+1] + in[k+1][j-1][i-1] + in[k+1][j-1][i+1] + in[k+1][j+1][i-1] + in[k+1][j+1][i+1]) + 
							c[3] * (in[k-1][j-1][i] + in[k-1][j][i-1] + in[k-1][j][i+1] + in[k-1][j+1][i] + in[k][j-1][i-1] + in[k][j-1][i+1] + in[k][j+1][i-1] + in[k][j+1][i+1] + in[k+1][j-1][i] + in[k+1][j][i-1] + in[k+1][j][i+1] + in[k+1][j+1][i]);
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
			for (k = 1; k < N-1; k+=4) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						out[k][j][i] = c[0] * in[k][j][i] + c[1] * (in[k-1][j][i] + in[k+1][j][i] + in[k][j-1][i] + in[k][j+1][i] + in[k][j][i-1] + in[k][j][i+1]) + c[2] * (in[k-1][j-1][i-1] + in[k-1][j-1][i+1] + in[k-1][j+1][i-1] + in[k-1][j+1][i+1] + in[k+1][j-1][i-1] + in[k+1][j-1][i+1] + in[k+1][j+1][i-1] + in[k+1][j+1][i+1]) + c[3] * (in[k-1][j-1][i] + in[k-1][j][i-1] + in[k-1][j][i+1] + in[k-1][j+1][i] + in[k][j-1][i-1] + in[k][j-1][i+1] + in[k][j+1][i-1] + in[k][j+1][i+1] + in[k+1][j-1][i] + in[k+1][j][i-1] + in[k+1][j][i+1] + in[k+1][j+1][i]);
						out[k+1][j][i] = c[0] * in[k+1][j][i] + c[1] * (in[k+1-1][j][i] + in[k+1+1][j][i] + in[k+1][j-1][i] + in[k+1][j+1][i] + in[k+1][j][i-1] + in[k+1][j][i+1]) + c[2] * (in[k+1-1][j-1][i-1] + in[k+1-1][j-1][i+1] + in[k+1-1][j+1][i-1] + in[k+1-1][j+1][i+1] + in[k+1+1][j-1][i-1] + in[k+1+1][j-1][i+1] + in[k+1+1][j+1][i-1] + in[k+1+1][j+1][i+1]) + c[3] * (in[k+1-1][j-1][i] + in[k+1-1][j][i-1] + in[k+1-1][j][i+1] + in[k+1-1][j+1][i] + in[k+1][j-1][i-1] + in[k+1][j-1][i+1] + in[k+1][j+1][i-1] + in[k+1][j+1][i+1] + in[k+1+1][j-1][i] + in[k+1+1][j][i-1] + in[k+1+1][j][i+1] + in[k+1+1][j+1][i]);
						out[k+2][j][i] = c[0] * in[k+2][j][i] + c[1] * (in[k+2-1][j][i] + in[k+2+1][j][i] + in[k+2][j-1][i] + in[k+2][j+1][i] + in[k+2][j][i-1] + in[k+2][j][i+1]) + c[2] * (in[k+2-1][j-1][i-1] + in[k+2-1][j-1][i+1] + in[k+2-1][j+1][i-1] + in[k+2-1][j+1][i+1] + in[k+2+1][j-1][i-1] + in[k+2+1][j-1][i+1] + in[k+2+1][j+1][i-1] + in[k+2+1][j+1][i+1]) + c[3] * (in[k+2-1][j-1][i] + in[k+2-1][j][i-1] + in[k+2-1][j][i+1] + in[k+2-1][j+1][i] + in[k+2][j-1][i-1] + in[k+2][j-1][i+1] + in[k+2][j+1][i-1] + in[k+2][j+1][i+1] + in[k+2+1][j-1][i] + in[k+2+1][j][i-1] + in[k+2+1][j][i+1] + in[k+2+1][j+1][i]);
						out[k+3][j][i] = c[0] * in[k+3][j][i] + c[1] * (in[k+3-1][j][i] + in[k+3+1][j][i] + in[k+3][j-1][i] + in[k+3][j+1][i] + in[k+3][j][i-1] + in[k+3][j][i+1]) + c[2] * (in[k+3-1][j-1][i-1] + in[k+3-1][j-1][i+1] + in[k+3-1][j+1][i-1] + in[k+3-1][j+1][i+1] + in[k+3+1][j-1][i-1] + in[k+3+1][j-1][i+1] + in[k+3+1][j+1][i-1] + in[k+3+1][j+1][i+1]) + c[3] * (in[k+3-1][j-1][i] + in[k+3-1][j][i-1] + in[k+3-1][j][i+1] + in[k+3-1][j+1][i] + in[k+3][j-1][i-1] + in[k+3][j-1][i+1] + in[k+3][j+1][i-1] + in[k+3][j+1][i+1] + in[k+3+1][j-1][i] + in[k+3+1][j][i-1] + in[k+3+1][j][i+1] + in[k+3+1][j+1][i]);
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("unroll: %lf\n", (double)512*512*512*30*5/(end_time - start_time)/1e9);
}
