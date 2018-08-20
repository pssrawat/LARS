#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void j3d125pt_opt (double*, double*, double*, int);

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
	int N = 516;
	double (*in)[516][516] = (double (*)[516][516]) getRandom3DArray (516, 516, 516);
	double (*out_ref)[516][516] = (double (*)[516][516]) getZero3DArray (516, 516, 516);
	double (*out)[516][516] = (double (*)[516][516]) getZero3DArray (516, 516, 516);
	double (*c)[5] = (double (*)[5]) getRandom2DArray (25, 5);

	int t, i, j, k, j0, k0, i0;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 2; k < N-2; k++) {
				for (j = 2; j < N-2; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 2; i < N-2; i++) {
						out_ref[k][j][i] = 0.1 *
							c[0][0]*in[k-2][j-2][i-2] + c[0][1]*in[k-2][j-2][i-1] + c[0][2]*in[k-2][j-2][i] + c[0][3]*in[k-2][j-2][i+1] + c[0][4]*in[k-2][j-2][i+2] +
							c[1][0]*in[k-2][j-1][i-2] + c[1][1]*in[k-2][j-1][i-1] + c[1][2]*in[k-2][j-1][i] + c[1][3]*in[k-2][j-1][i+1] + c[1][4]*in[k-2][j-1][i+2] +
							c[2][0]*in[k-2][j][i-2] +   c[2][1]*in[k-2][j][i-1] +   c[2][2]*in[k-2][j][i] +   c[2][3]*in[k-2][j][i+1] +   c[2][4]*in[k-2][j][i+2] +
							c[3][0]*in[k-2][j+1][i-2] + c[3][1]*in[k-2][j+1][i-1] + c[3][2]*in[k-2][j+1][i] + c[3][3]*in[k-2][j+1][i+1] + c[3][4]*in[k-2][j+1][i+2] +
							c[4][0]*in[k-2][j+2][i-2] + c[4][1]*in[k-2][j+2][i-1] + c[4][2]*in[k-2][j+2][i] + c[4][3]*in[k-2][j+2][i+1] + c[4][4]*in[k-2][j+2][i+2] +

							c[5][0]*in[k-1][j-2][i-2] + c[5][1]*in[k-1][j-2][i-1] + c[5][2]*in[k-1][j-2][i] + c[5][3]*in[k-1][j-2][i+1] + c[5][4]*in[k-1][j-2][i+2] +
							c[6][0]*in[k-1][j-1][i-2] + c[6][1]*in[k-1][j-1][i-1] + c[6][2]*in[k-1][j-1][i] + c[6][3]*in[k-1][j-1][i+1] + c[6][4]*in[k-1][j-1][i+2] +
							c[7][0]*in[k-1][j][i-2] +   c[7][1]*in[k-1][j][i-1] +   c[7][2]*in[k-1][j][i] +   c[7][3]*in[k-1][j][i+1] +   c[7][4]*in[k-1][j][i+2] +
							c[8][0]*in[k-1][j+1][i-2] + c[8][1]*in[k-1][j+1][i-1] + c[8][2]*in[k-1][j+1][i] + c[8][3]*in[k-1][j+1][i+1] + c[8][4]*in[k-1][j+1][i+2] +
							c[9][0]*in[k-1][j+2][i-2] + c[9][1]*in[k-1][j+2][i-1] + c[9][2]*in[k-1][j+2][i] + c[9][3]*in[k-1][j+2][i+1] + c[9][4]*in[k-1][j+2][i+2] +

							c[10][0]*in[k  ][j-2][i-2] + c[10][1]*in[k  ][j-2][i-1] + c[10][2]*in[k  ][j-2][i] + c[10][3]*in[k  ][j-2][i+1] + c[10][4]*in[k  ][j-2][i+2] +
							c[11][0]*in[k  ][j-1][i-2] + c[11][1]*in[k  ][j-1][i-1] + c[11][2]*in[k  ][j-1][i] + c[11][3]*in[k  ][j-1][i+1] + c[11][4]*in[k  ][j-1][i+2] +
							c[11][0]*in[k  ][j][i-2] +   c[11][1]*in[k  ][j][i-1] +   c[11][2]*in[k  ][j][i] +   c[11][3]*in[k  ][j][i+1] +   c[11][4]*in[k  ][j][i+2] +
							c[13][0]*in[k  ][j+1][i-2] + c[13][1]*in[k  ][j+1][i-1] + c[13][2]*in[k  ][j+1][i] + c[13][3]*in[k  ][j+1][i+1] + c[13][4]*in[k  ][j+1][i+2] +
							c[14][0]*in[k  ][j+2][i-2] + c[14][1]*in[k  ][j+2][i-1] + c[14][2]*in[k  ][j+2][i] + c[14][3]*in[k  ][j+2][i+1] + c[14][4]*in[k  ][j+2][i+2] +

							c[15][0]*in[k+1][j-2][i-2] + c[15][1]*in[k+1][j-2][i-1] + c[15][2]*in[k+1][j-2][i] + c[15][3]*in[k+1][j-2][i+1] + c[15][4]*in[k+1][j-2][i+2] +
							c[16][0]*in[k+1][j-1][i-2] + c[16][1]*in[k+1][j-1][i-1] + c[16][2]*in[k+1][j-1][i] + c[16][3]*in[k+1][j-1][i+1] + c[16][4]*in[k+1][j-1][i+2] +
							c[17][0]*in[k+1][j][i-2] +   c[17][1]*in[k+1][j][i-1] +   c[17][2]*in[k+1][j][i] +   c[17][3]*in[k+1][j][i+1] +   c[17][4]*in[k+1][j][i+2] +
							c[18][0]*in[k+1][j+1][i-2] + c[18][1]*in[k+1][j+1][i-1] + c[18][2]*in[k+1][j+1][i] + c[18][3]*in[k+1][j+1][i+1] + c[18][4]*in[k+1][j+1][i+2] +
							c[19][0]*in[k+1][j+2][i-2] + c[19][1]*in[k+1][j+2][i-1] + c[19][2]*in[k+1][j+2][i] + c[19][3]*in[k+1][j+2][i+1] + c[19][4]*in[k+1][j+2][i+2] +

							c[20][0]*in[k+2][j-2][i-2] + c[20][1]*in[k+2][j-2][i-1] + c[20][2]*in[k+2][j-2][i] + c[20][3]*in[k+2][j-2][i+1] + c[20][4]*in[k+2][j-2][i+2] +
							c[21][0]*in[k+2][j-1][i-2] + c[21][1]*in[k+2][j-1][i-1] + c[21][2]*in[k+2][j-1][i] + c[21][3]*in[k+2][j-1][i+1] + c[21][4]*in[k+2][j-1][i+2] +
							c[22][0]*in[k+2][j][i-2] +   c[22][1]*in[k+2][j][i-1] +   c[22][2]*in[k+2][j][i] +   c[22][3]*in[k+2][j][i+1] +   c[22][4]*in[k+2][j][i+2] +
							c[23][0]*in[k+2][j+1][i-2] + c[23][1]*in[k+2][j+1][i-1] + c[23][2]*in[k+2][j+1][i] + c[23][3]*in[k+2][j+1][i+1] + c[23][4]*in[k+2][j+1][i+2] +
							c[24][0]*in[k+2][j+2][i-2] + c[24][1]*in[k+2][j+2][i-1] + c[24][2]*in[k+2][j+2][i] + c[24][3]*in[k+2][j+2][i+1] + c[24][4]*in[k+2][j+2][i+2];
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
			for (k = 2; k < N-2; k++) {
				for (j = 2; j < N-2; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 2; i < N-2; i++) {
						out_ref[k][j][i] =
							c[0][0]*in[k-2][j-2][i-2] + c[0][1]*in[k-2][j-2][i-1] + c[0][2]*in[k-2][j-2][i] + c[0][3]*in[k-2][j-2][i+1] + c[0][4]*in[k-2][j-2][i+2] +
							c[1][0]*in[k-2][j-1][i-2] + c[1][1]*in[k-2][j-1][i-1] + c[1][2]*in[k-2][j-1][i] + c[1][3]*in[k-2][j-1][i+1] + c[1][4]*in[k-2][j-1][i+2] +
							c[2][0]*in[k-2][j][i-2] +   c[2][1]*in[k-2][j][i-1] +   c[2][2]*in[k-2][j][i] +   c[2][3]*in[k-2][j][i+1] +   c[2][4]*in[k-2][j][i+2] +
							c[3][0]*in[k-2][j+1][i-2] + c[3][1]*in[k-2][j+1][i-1] + c[3][2]*in[k-2][j+1][i] + c[3][3]*in[k-2][j+1][i+1] + c[3][4]*in[k-2][j+1][i+2] +
							c[4][0]*in[k-2][j+2][i-2] + c[4][1]*in[k-2][j+2][i-1] + c[4][2]*in[k-2][j+2][i] + c[4][3]*in[k-2][j+2][i+1] + c[4][4]*in[k-2][j+2][i+2] +

							c[5][0]*in[k-1][j-2][i-2] + c[5][1]*in[k-1][j-2][i-1] + c[5][2]*in[k-1][j-2][i] + c[5][3]*in[k-1][j-2][i+1] + c[5][4]*in[k-1][j-2][i+2] +
							c[6][0]*in[k-1][j-1][i-2] + c[6][1]*in[k-1][j-1][i-1] + c[6][2]*in[k-1][j-1][i] + c[6][3]*in[k-1][j-1][i+1] + c[6][4]*in[k-1][j-1][i+2] +
							c[7][0]*in[k-1][j][i-2] +   c[7][1]*in[k-1][j][i-1] +   c[7][2]*in[k-1][j][i] +   c[7][3]*in[k-1][j][i+1] +   c[7][4]*in[k-1][j][i+2] +
							c[8][0]*in[k-1][j+1][i-2] + c[8][1]*in[k-1][j+1][i-1] + c[8][2]*in[k-1][j+1][i] + c[8][3]*in[k-1][j+1][i+1] + c[8][4]*in[k-1][j+1][i+2] +
							c[9][0]*in[k-1][j+2][i-2] + c[9][1]*in[k-1][j+2][i-1] + c[9][2]*in[k-1][j+2][i] + c[9][3]*in[k-1][j+2][i+1] + c[9][4]*in[k-1][j+2][i+2] +

							c[10][0]*in[k  ][j-2][i-2] + c[10][1]*in[k  ][j-2][i-1] + c[10][2]*in[k  ][j-2][i] + c[10][3]*in[k  ][j-2][i+1] + c[10][4]*in[k  ][j-2][i+2] +
							c[11][0]*in[k  ][j-1][i-2] + c[11][1]*in[k  ][j-1][i-1] + c[11][2]*in[k  ][j-1][i] + c[11][3]*in[k  ][j-1][i+1] + c[11][4]*in[k  ][j-1][i+2] +
							c[11][0]*in[k  ][j][i-2] +   c[11][1]*in[k  ][j][i-1] +   c[11][2]*in[k  ][j][i] +   c[11][3]*in[k  ][j][i+1] +   c[11][4]*in[k  ][j][i+2] +
							c[13][0]*in[k  ][j+1][i-2] + c[13][1]*in[k  ][j+1][i-1] + c[13][2]*in[k  ][j+1][i] + c[13][3]*in[k  ][j+1][i+1] + c[13][4]*in[k  ][j+1][i+2] +
							c[14][0]*in[k  ][j+2][i-2] + c[14][1]*in[k  ][j+2][i-1] + c[14][2]*in[k  ][j+2][i] + c[14][3]*in[k  ][j+2][i+1] + c[14][4]*in[k  ][j+2][i+2] +

							c[15][0]*in[k+1][j-2][i-2] + c[15][1]*in[k+1][j-2][i-1] + c[15][2]*in[k+1][j-2][i] + c[15][3]*in[k+1][j-2][i+1] + c[15][4]*in[k+1][j-2][i+2] +
							c[16][0]*in[k+1][j-1][i-2] + c[16][1]*in[k+1][j-1][i-1] + c[16][2]*in[k+1][j-1][i] + c[16][3]*in[k+1][j-1][i+1] + c[16][4]*in[k+1][j-1][i+2] +
							c[17][0]*in[k+1][j][i-2] +   c[17][1]*in[k+1][j][i-1] +   c[17][2]*in[k+1][j][i] +   c[17][3]*in[k+1][j][i+1] +   c[17][4]*in[k+1][j][i+2] +
							c[18][0]*in[k+1][j+1][i-2] + c[18][1]*in[k+1][j+1][i-1] + c[18][2]*in[k+1][j+1][i] + c[18][3]*in[k+1][j+1][i+1] + c[18][4]*in[k+1][j+1][i+2] +
							c[19][0]*in[k+1][j+2][i-2] + c[19][1]*in[k+1][j+2][i-1] + c[19][2]*in[k+1][j+2][i] + c[19][3]*in[k+1][j+2][i+1] + c[19][4]*in[k+1][j+2][i+2] +

							c[20][0]*in[k+2][j-2][i-2] + c[20][1]*in[k+2][j-2][i-1] + c[20][2]*in[k+2][j-2][i] + c[20][3]*in[k+2][j-2][i+1] + c[20][4]*in[k+2][j-2][i+2] +
							c[21][0]*in[k+2][j-1][i-2] + c[21][1]*in[k+2][j-1][i-1] + c[21][2]*in[k+2][j-1][i] + c[21][3]*in[k+2][j-1][i+1] + c[21][4]*in[k+2][j-1][i+2] +
							c[22][0]*in[k+2][j][i-2] +   c[22][1]*in[k+2][j][i-1] +   c[22][2]*in[k+2][j][i] +   c[22][3]*in[k+2][j][i+1] +   c[22][4]*in[k+2][j][i+2] +
							c[23][0]*in[k+2][j+1][i-2] + c[23][1]*in[k+2][j+1][i-1] + c[23][2]*in[k+2][j+1][i] + c[23][3]*in[k+2][j+1][i+1] + c[23][4]*in[k+2][j+1][i+2] +
							c[24][0]*in[k+2][j+2][i-2] + c[24][1]*in[k+2][j+2][i-1] + c[24][2]*in[k+2][j+2][i] + c[24][3]*in[k+2][j+2][i+1] + c[24][4]*in[k+2][j+2][i+2];
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)512*512*512*249*5/(end_time - start_time)/1e9);

	j3d125pt_opt ((double*)in, (double*)out, (double*)c, N);

	double error = checkError3D (N, N, 0, (double*)out, (double*) out_ref, 2, N-2, 2, N-2, 2, N-2);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
