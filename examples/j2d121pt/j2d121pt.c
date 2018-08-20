#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void j2d121pt_opt (double*, double*, double*, int);

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
	int N = 8202;
	double (*in)[8202] = (double (*)[8202]) getRandom2DArray (8202, 8202);
	double (*out_ref)[8202] = (double (*)[8202]) getZero2DArray (8202, 8202);
	double (*out)[8202] = (double (*)[8202]) getZero2DArray (8202, 8202);
	double (*c)[11] = (double (*)[11]) getRandom2DArray (11, 11);

	int t, i, j;
	double start_time, end_time;

	// Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 5; j < N-5; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 5; i < N-5; i++) {
					out_ref[j][i] = 0.1 * 
						c[0][0]*in[j-5][i-5]  +   c[0][1]*in[j-5][i-4]  +  c[0][2]*in[j-5][i-3]  +  c[0][3]*in[j-5][i-2]  -  c[0][4]*in[j-5][i-1]  + c[0][5]*in[j-5][i]  - c[0][6]*in[j-5][i+1]  + c[0][7]*in[j-5][i+2]  - c[0][8]*in[j-5][i+3]  + c[0][9]*in[j-5][i+4]  + c[0][10]*in[j-5][i+5] -
						c[1][0]*in[j-4][i-5]  +   c[1][1]*in[j-4][i-4]  +  c[1][2]*in[j-4][i-3]  +  c[1][3]*in[j-4][i-2]  -  c[1][4]*in[j-4][i-1]  + c[1][5]*in[j-4][i]  - c[1][6]*in[j-4][i+1]  + c[1][7]*in[j-4][i+2]  - c[1][8]*in[j-4][i+3]  + c[1][9]*in[j-4][i+4]  + c[1][10]*in[j-4][i+5] -
						c[2][0]*in[j-3][i-5]  +   c[2][1]*in[j-3][i-4]  +  c[2][2]*in[j-3][i-3]  +  c[2][3]*in[j-3][i-2]  -  c[2][4]*in[j-3][i-1]  + c[2][5]*in[j-3][i]  - c[2][6]*in[j-3][i+1]  + c[2][7]*in[j-3][i+2]  - c[2][8]*in[j-3][i+3]  + c[2][9]*in[j-3][i+4]  + c[2][10]*in[j-3][i+5] -
						c[3][0]*in[j-2][i-5]  +   c[3][1]*in[j-2][i-4]  +  c[3][2]*in[j-2][i-3]  +  c[3][3]*in[j-2][i-2]  -  c[3][4]*in[j-2][i-1]  + c[3][5]*in[j-2][i]  - c[3][6]*in[j-2][i+1]  + c[3][7]*in[j-2][i+2]  - c[3][8]*in[j-2][i+3]  + c[3][9]*in[j-2][i+4]  + c[3][10]*in[j-2][i+5] -
						c[4][0]*in[j-1][i-5]  +   c[4][1]*in[j-1][i-4]  +  c[4][2]*in[j-1][i-3]  +  c[4][3]*in[j-1][i-2]  -  c[4][4]*in[j-1][i-1]  + c[4][5]*in[j-1][i]  - c[4][6]*in[j-1][i+1]  + c[4][7]*in[j-1][i+2]  - c[4][8]*in[j-1][i+3]  + c[4][9]*in[j-1][i+4]  + c[4][10]*in[j-1][i+5] -
						c[5][0]*in[j][i-5] +      c[5][1]*in[j][i-4] +     c[5][2]*in[j][i-3] +     c[5][3]*in[j][i-2] -     c[5][4]*in[j  ][i-1] +    c[5][5]*in[j][i] -    c[5][6]*in[j][i+1] +    c[5][7]*in[j][i+2] -    c[5][8]*in[j][i+3] +    c[5][9]*in[j][i+4] +    c[5][10]*in[j][i+5] -
						c[6][0]*in[j+1][i-5]  +   c[6][1]*in[j+1][i-4]  +  c[6][2]*in[j+1][i-3]  +  c[6][3]*in[j+1][i-2]  -  c[6][4]*in[j+1][i-1]  + c[6][5]*in[j+1][i]  - c[6][6]*in[j+1][i+1]  + c[6][7]*in[j+1][i+2]  - c[6][8]*in[j+1][i+3]  + c[6][9]*in[j+1][i+4]  + c[6][10]*in[j+1][i+5] -
						c[7][0]*in[j+2][i-5]  +   c[7][1]*in[j+2][i-4]  +  c[7][2]*in[j+2][i-3]  +  c[7][3]*in[j+2][i-2]  -  c[7][4]*in[j+2][i-1]  + c[7][5]*in[j+2][i]  - c[7][6]*in[j+2][i+1]  + c[7][7]*in[j+2][i+2]  - c[7][8]*in[j+2][i+3]  + c[7][9]*in[j+2][i+4]  + c[7][10]*in[j+2][i+5] -
						c[8][0]*in[j+3][i-5]  +   c[8][1]*in[j+3][i-4]  +  c[8][2]*in[j+3][i-3]  +  c[8][3]*in[j+3][i-2]  -  c[8][4]*in[j+3][i-1]  + c[8][5]*in[j+3][i]  - c[8][6]*in[j+3][i+1]  + c[8][7]*in[j+3][i+2]  - c[8][8]*in[j+3][i+3]  + c[8][9]*in[j+3][i+4]  + c[8][10]*in[j+3][i+5] -
						c[9][0]*in[j+4][i-5]  +   c[9][1]*in[j+4][i-4]  +  c[9][2]*in[j+4][i-3]  +  c[9][3]*in[j+4][i-2]  -  c[9][4]*in[j+4][i-1]  + c[9][5]*in[j+4][i]  - c[9][6]*in[j+4][i+1]  + c[9][7]*in[j+4][i+2]  - c[9][8]*in[j+4][i+3]  + c[9][9]*in[j+4][i+4]  + c[9][10]*in[j+4][i+5] -
						c[10][0]*in[j+5][i-5] +  c[10][1]*in[j+5][i-4] +  c[10][2]*in[j+5][i-3] +  c[10][3]*in[j+5][i-2] -  c[10][4]*in[j+5][i-1] +  c[10][5]*in[j+5][i] -  c[10][6]*in[j+5][i+1] +  c[10][7]*in[j+5][i+2] -  c[10][8]*in[j+5][i+3] +  c[10][9]*in[j+5][i+4] +  c[10][10]*in[j+5][i+5];
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 5; j < N-5; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 5; i < N-5; i++) {
					out_ref[j][i] =
						c[0][0]*in[j-5][i-5]  +   c[0][1]*in[j-5][i-4]  +  c[0][2]*in[j-5][i-3]  +  c[0][3]*in[j-5][i-2]  -  c[0][4]*in[j-5][i-1]  + c[0][5]*in[j-5][i]  - c[0][6]*in[j-5][i+1]  + c[0][7]*in[j-5][i+2]  - c[0][8]*in[j-5][i+3]  + c[0][9]*in[j-5][i+4]  + c[0][10]*in[j-5][i+5] -
						c[1][0]*in[j-4][i-5]  +   c[1][1]*in[j-4][i-4]  +  c[1][2]*in[j-4][i-3]  +  c[1][3]*in[j-4][i-2]  -  c[1][4]*in[j-4][i-1]  + c[1][5]*in[j-4][i]  - c[1][6]*in[j-4][i+1]  + c[1][7]*in[j-4][i+2]  - c[1][8]*in[j-4][i+3]  + c[1][9]*in[j-4][i+4]  + c[1][10]*in[j-4][i+5] -
						c[2][0]*in[j-3][i-5]  +   c[2][1]*in[j-3][i-4]  +  c[2][2]*in[j-3][i-3]  +  c[2][3]*in[j-3][i-2]  -  c[2][4]*in[j-3][i-1]  + c[2][5]*in[j-3][i]  - c[2][6]*in[j-3][i+1]  + c[2][7]*in[j-3][i+2]  - c[2][8]*in[j-3][i+3]  + c[2][9]*in[j-3][i+4]  + c[2][10]*in[j-3][i+5] -
						c[3][0]*in[j-2][i-5]  +   c[3][1]*in[j-2][i-4]  +  c[3][2]*in[j-2][i-3]  +  c[3][3]*in[j-2][i-2]  -  c[3][4]*in[j-2][i-1]  + c[3][5]*in[j-2][i]  - c[3][6]*in[j-2][i+1]  + c[3][7]*in[j-2][i+2]  - c[3][8]*in[j-2][i+3]  + c[3][9]*in[j-2][i+4]  + c[3][10]*in[j-2][i+5] -
						c[4][0]*in[j-1][i-5]  +   c[4][1]*in[j-1][i-4]  +  c[4][2]*in[j-1][i-3]  +  c[4][3]*in[j-1][i-2]  -  c[4][4]*in[j-1][i-1]  + c[4][5]*in[j-1][i]  - c[4][6]*in[j-1][i+1]  + c[4][7]*in[j-1][i+2]  - c[4][8]*in[j-1][i+3]  + c[4][9]*in[j-1][i+4]  + c[4][10]*in[j-1][i+5] -
						c[5][0]*in[j][i-5] +      c[5][1]*in[j][i-4] +     c[5][2]*in[j][i-3] +     c[5][3]*in[j][i-2] -     c[5][4]*in[j  ][i-1] +    c[5][5]*in[j][i] -    c[5][6]*in[j][i+1] +    c[5][7]*in[j][i+2] -    c[5][8]*in[j][i+3] +    c[5][9]*in[j][i+4] +    c[5][10]*in[j][i+5] -
						c[6][0]*in[j+1][i-5]  +   c[6][1]*in[j+1][i-4]  +  c[6][2]*in[j+1][i-3]  +  c[6][3]*in[j+1][i-2]  -  c[6][4]*in[j+1][i-1]  + c[6][5]*in[j+1][i]  - c[6][6]*in[j+1][i+1]  + c[6][7]*in[j+1][i+2]  - c[6][8]*in[j+1][i+3]  + c[6][9]*in[j+1][i+4]  + c[6][10]*in[j+1][i+5] -
						c[7][0]*in[j+2][i-5]  +   c[7][1]*in[j+2][i-4]  +  c[7][2]*in[j+2][i-3]  +  c[7][3]*in[j+2][i-2]  -  c[7][4]*in[j+2][i-1]  + c[7][5]*in[j+2][i]  - c[7][6]*in[j+2][i+1]  + c[7][7]*in[j+2][i+2]  - c[7][8]*in[j+2][i+3]  + c[7][9]*in[j+2][i+4]  + c[7][10]*in[j+2][i+5] -
						c[8][0]*in[j+3][i-5]  +   c[8][1]*in[j+3][i-4]  +  c[8][2]*in[j+3][i-3]  +  c[8][3]*in[j+3][i-2]  -  c[8][4]*in[j+3][i-1]  + c[8][5]*in[j+3][i]  - c[8][6]*in[j+3][i+1]  + c[8][7]*in[j+3][i+2]  - c[8][8]*in[j+3][i+3]  + c[8][9]*in[j+3][i+4]  + c[8][10]*in[j+3][i+5] -
						c[9][0]*in[j+4][i-5]  +   c[9][1]*in[j+4][i-4]  +  c[9][2]*in[j+4][i-3]  +  c[9][3]*in[j+4][i-2]  -  c[9][4]*in[j+4][i-1]  + c[9][5]*in[j+4][i]  - c[9][6]*in[j+4][i+1]  + c[9][7]*in[j+4][i+2]  - c[9][8]*in[j+4][i+3]  + c[9][9]*in[j+4][i+4]  + c[9][10]*in[j+4][i+5] -
						c[10][0]*in[j+5][i-5] +  c[10][1]*in[j+5][i-4] +  c[10][2]*in[j+5][i-3] +  c[10][3]*in[j+5][i-2] -  c[10][4]*in[j+5][i-1] +  c[10][5]*in[j+5][i] -  c[10][6]*in[j+5][i+1] +  c[10][7]*in[j+5][i+2] -  c[10][8]*in[j+5][i+3] +  c[10][9]*in[j+5][i+4] +  c[10][10]*in[j+5][i+5];
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)8192*8192*241*5/(end_time - start_time)/1e9);

	j2d121pt_opt ((double*)in, (double*)out, (double*)c, N);

	double error = checkError2D (N, 0, (double*)out, (double*) out_ref, 5, N-5, 5, N-5);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
