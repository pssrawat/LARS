#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void j3d7pt_opt (double * __restrict__ t_out, double * __restrict__ t_in, int N) {
	double (*out)[514][514] = (double (*)[514][514]) t_out;
	double (*in)[514][514] = (double (*)[514][514]) t_in;
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
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma clang loop vectorize (enable) interleave(enable)
#pragma GCC ivdep
					for (i = 1; i < N-1; i++) {
						out[k][j][i] = 0.1 * a*in[k][j][i] - b*h2inv*(in[k][j][i+1] + in[k][j][i-1] + in[k][j+1][i] + in[k][j-1][i] + in[k+1][j][i] + in[k-1][j][i]) -in[k][j][i]*6.0;
					}
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<10; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k+=4) {
				for (j = 1; j < N-1; j+=2) {
#pragma clang loop vectorize (enable) interleave(enable)
#pragma GCC ivdep
					for (i = 1; i < N-1; i++) {
						out[k][j][i] = a*in[k][j][i] - b*h2inv*(in[k][j][i+1] + in[k][j][i-1] + in[k][j+1][i] + in[k][j-1][i] + in[k+1][j][i] + in[k-1][j][i]) -in[k][j][i]*6.0;
						out[k][j+1][i] = a*in[k][j+1][i] - b*h2inv*(in[k][j+1][i+1] + in[k][j+1][i-1] + in[k][j+1+1][i] + in[k][j+1-1][i] + in[k+1][j+1][i] + in[k-1][j+1][i]) -in[k][j+1][i]*6.0;
						out[k+1][j][i] = a*in[k+1][j][i] - b*h2inv*(in[k+1][j][i+1] + in[k+1][j][i-1] + in[k+1][j+1][i] + in[k+1][j-1][i] + in[k+1+1][j][i] + in[k+1-1][j][i]) -in[k+1][j][i]*6.0;
						out[k+1][j+1][i] = a*in[k+1][j+1][i] - b*h2inv*(in[k+1][j+1][i+1] + in[k+1][j+1][i-1] + in[k+1][j+1+1][i] + in[k+1][j+1-1][i] + in[k+1+1][j+1][i] + in[k+1-1][j+1][i]) -in[k+1][j+1][i]*6.0;
						out[k+2][j][i] = a*in[k+2][j][i] - b*h2inv*(in[k+2][j][i+1] + in[k+2][j][i-1] + in[k+2][j+1][i] + in[k+2][j-1][i] + in[k+2+1][j][i] + in[k+2-1][j][i]) -in[k+2][j][i]*6.0;
						out[k+2][j+1][i] = a*in[k+2][j+1][i] - b*h2inv*(in[k+2][j+1][i+1] + in[k+2][j+1][i-1] + in[k+2][j+1+1][i] + in[k+2][j+1-1][i] + in[k+2+1][j+1][i] + in[k+2-1][j+1][i]) -in[k+2][j+1][i]*6.0;
						out[k+3][j][i] = a*in[k+3][j][i] - b*h2inv*(in[k+3][j][i+1] + in[k+3][j][i-1] + in[k+3][j+1][i] + in[k+3][j-1][i] + in[k+3+1][j][i] + in[k+3-1][j][i]) -in[k+3][j][i]*6.0;
						out[k+3][j+1][i] = a*in[k+3][j+1][i] - b*h2inv*(in[k+3][j+1][i+1] + in[k+3][j+1][i-1] + in[k+3][j+1+1][i] + in[k+3][j+1-1][i] + in[k+3+1][j+1][i] + in[k+3-1][j+1][i]) -in[k+3][j+1][i]*6.0;
					}
				}
			}
		}
	}
	end_time = rtclock ();
        printf ("unroll: %6lf\n", (double)512*512*512*11*10/(end_time - start_time)/1e9);
}
