#include "../common/common.h"
#include <stdio.h>
#include<stdlib.h>

#define tf (3d0/4)
#define i6 (1d0/6)
#define i144 (1d0/144.0)
#define i12 (1d0/12)
#define d4a (2d0/3)
#define d4b (-(1d0/12))


extern void host_code (double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int);

void flush_llc() {
	int size_in_mb = 64;
	//  fprintf(stderr, "Flushing last level cache... ");
	// Seed random # gen
	struct timezone Tzp;
	struct timeval Tp;
	int stat;
	stat = gettimeofday (&Tp, &Tzp);
	if (stat != 0) fprintf(stderr, "Error return from gettimeofday: %d",stat);
	srand(Tp.tv_usec);

	// Allocate a big array, set first element to random #, traverse it twice
	int i,j;
	double *flush = (double*)malloc(size_in_mb*1024*128*sizeof(double));
	flush[0] = (rand() % 128) * (((double)rand()) / RAND_MAX) + 1;
	for (i = 0; i < 2; i++) {
		for (j = 1; j < size_in_mb*1024*128; j++) {
			flush[j] = flush[j-1]*1.00000000000000001;
		}
	}
}

int main(int argc, char** argv) {
	int N = 324; 

	double (*u_0)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*u_1)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*u_2)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*mu)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*la)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double *strx = (double *) getRandom1DArray(324);
	double *stry = (double *) getRandom1DArray(324);
	double *strz = (double *) getRandom1DArray(324);
	double (*uacc_opt_0)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*uacc_opt_1)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*uacc_opt_2)[324][324] = (double (*)[324][324]) getRandom3DArray (324, 324, 324);
	double (*uacc_0)[324][324] = (double (*)[324][324]) getZero3DArray(324, 324, 324);
	double (*uacc_1)[324][324] = (double (*)[324][324]) getZero3DArray(324, 324, 324);
	double (*uacc_2)[324][324] = (double (*)[324][324]) getZero3DArray(324, 324, 324);
	memcpy(uacc_0, uacc_opt_0, sizeof(double)*324*324*324);
	memcpy(uacc_1, uacc_opt_1, sizeof(double)*324*324*324);
	memcpy(uacc_2, uacc_opt_2, sizeof(double)*324*324*324);

	int t, i, j, k;
	double start_time, end_time;

	double a1 = 1.2;
	double h = 1.0/(N-1);
	double cof = 1e0 / ( h *  h);


	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel for private (j,i)
		for (k = 2; k < N-2; k++) {
			for (j = 2; j < N-2; j++) {
#pragma GCC ivdep 
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 2; i < N-2; i++) {
					double mux1 = mu[k][j][i-1] * strx[i-1] - 3e0 / 4 * (mu[k][j][i] * strx[i] + mu[k][j][i-2] * strx[i-2]);
					double mux2 = mu[k][j][i-2] * strx[i-2] + mu[k][j][i+1] * strx[i+1] + 3 * (mu[k][j][i] * strx[i] + mu[k][j][i-1] * strx[i-1]);
					double mux3 = mu[k][j][i-1] * strx[i-1] + mu[k][j][i+2] * strx[i+2] + 3 * (mu[k][j][i+1] * strx[i+1] + mu[k][j][i] * strx[i]);
					double mux4 = mu[k][j][i+1] * strx[i+1] - 3e0 / 4 * (mu[k][j][i] * strx[i] + mu[k][j][i+2] * strx[i+2]);
					double muy1 = mu[k][j-1][i] * stry[j-1] - 3e0 / 4 * (mu[k][j][i] * stry[j] + mu[k][j-2][i] * stry[j-2]);
					double muy2 = mu[k][j-2][i] * stry[j-2] + mu[k][j+1][i] * stry[j+1] + 3 * (mu[k][j][i] * stry[j] + mu[k][j-1][i] * stry[j-1]);
					double muy3 = mu[k][j-1][i] * stry[j-1] + mu[k][j+2][i] * stry[j+2] + 3 * (mu[k][j+1][i] * stry[j+1] + mu[k][j][i] * stry[j]);
					double muy4 = mu[k][j+1][i] * stry[j+1] - 3e0 / 4 * (mu[k][j][i] * stry[j] + mu[k][j+2][i] * stry[j+2]);
					double muz1 = mu[k-1][j][i] * strz[k-1] - 3e0 / 4 * (mu[k][j][i] * strz[k] + mu[k-2][j][i] * strz[k-2]);
					double muz2 = mu[k-2][j][i] * strz[k-2] + mu[k+1][j][i] * strz[k+1] + 3 * (mu[k][j][i] * strz[k] + mu[k-1][j][i] * strz[k-1]);
					double muz3 = mu[k-1][j][i] * strz[k-1] + mu[k+2][j][i] * strz[k+2] + 3 * (mu[k+1][j][i] * strz[k+1] + mu[k][j][i] * strz[k]);
					double muz4 = mu[k+1][j][i] * strz[k+1] - 3e0 / 4 * (mu[k][j][i] * strz[k] + mu[k+2][j][i] * strz[k+2]);

					double r1 = 1e0 / 6 * (strx[i] * ((2 * mux1 + la[k][j][i-1] * strx[i-1] - 3e0 / 4 * (la[k][j][i] * strx[i] + la[k][j][i-2] * strx[i-2])) * (u_0[k][j][i-2] - u_0[k][j][i]) + (2 * mux2 + la[k][j][i-2] * strx[i-2] + la[k][j][i+1] * strx[i+1] + 3 * (la[k][j][i] * strx[i] + la[k][j][i-1] * strx[i-1])) * (u_0[k][j][i-1] - u_0[k][j][i]) + (2 * mux3 + la[k][j][i-1] * strx[i-1] + la[k][j][i+2] * strx[i+2] + 3 * (la[k][j][i+1] * strx[i+1] + la[k][j][i] * strx[i])) * (u_0[k][j][i+1] - u_0[k][j][i]) + (2 * mux4 + la[k][j][i+1] * strx[i+1] - 3e0 / 4 * (la[k][j][i] * strx[i] + la[k][j][i+2] * strx[i+2])) * (u_0[k][j][i+2] - u_0[k][j][i])) + stry[j] * (muy1 * (u_0[k][j-2][i] - u_0[k][j][i]) + muy2 * (u_0[k][j-1][i] - u_0[k][j][i]) + muy3 * (u_0[k][j+1][i] - u_0[k][j][i]) + muy4 * (u_0[k][j+2][i] - u_0[k][j][i])) + strz[k] * (muz1 * (u_0[k-2][j][i] - u_0[k][j][i]) + muz2 * (u_0[k-1][j][i] - u_0[k][j][i]) + muz3 * (u_0[k+1][j][i] - u_0[k][j][i]) + muz4 * (u_0[k+2][j][i] - u_0[k][j][i])));
					double r2 = 1e0 / 6 * (strx[i] * (mux1 * (u_1[k][j][i-2] - u_1[k][j][i]) + mux2 * (u_1[k][j][i-1] - u_1[k][j][i]) + mux3 * (u_1[k][j][i+1] - u_1[k][j][i]) + mux4 * (u_1[k][j][i+2] - u_1[k][j][i])) + stry[j] * ((2 * muy1 + la[k][j-1][i] * stry[j-1] - 3e0 / 4 * (la[k][j][i] * stry[j] + la[k][j-2][i] * stry[j-2])) * (u_1[k][j-2][i] - u_1[k][j][i]) + (2 * muy2 + la[k][j-2][i] * stry[j-2] + la[k][j+1][i] * stry[j+1] + 3 * (la[k][j][i] * stry[j] + la[k][j-1][i] * stry[j-1])) * (u_1[k][j-1][i] - u_1[k][j][i]) + (2 * muy3 + la[k][j-1][i] * stry[j-1] + la[k][j+2][i] * stry[j+2] + 3 * (la[k][j+1][i] * stry[j+1] + la[k][j][i] * stry[j])) * (u_1[k][j+1][i] - u_1[k][j][i]) + (2 * muy4 + la[k][j+1][i] * stry[j+1] - 3e0 / 4 * (la[k][j][i] * stry[j] + la[k][j+2][i] * stry[j+2])) * (u_1[k][j+2][i] - u_1[k][j][i])) + strz[k] * (muz1 * (u_1[k-2][j][i] - u_1[k][j][i]) + muz2 * (u_1[k-1][j][i] - u_1[k][j][i]) + muz3 * (u_1[k+1][j][i] - u_1[k][j][i]) + muz4 * (u_1[k+2][j][i] - u_1[k][j][i])));
					double r3 = 1e0 / 6 * (strx[i] * (mux1 * (u_2[k][j][i-2] - u_2[k][j][i]) + mux2 * (u_2[k][j][i-1] - u_2[k][j][i]) + mux3 * (u_2[k][j][i+1] - u_2[k][j][i]) + mux4 * (u_2[k][j][i+2] - u_2[k][j][i])) + stry[j] * (muy1 * (u_2[k][j-2][i] - u_2[k][j][i]) + muy2 * (u_2[k][j-1][i] - u_2[k][j][i]) + muy3 * (u_2[k][j+1][i] - u_2[k][j][i]) + muy4 * (u_2[k][j+2][i] - u_2[k][j][i])) + strz[k] * ((2 * muz1 + la[k-1][j][i] * strz[k-1] - 3e0 / 4 * (la[k][j][i] * strz[k] + la[k-2][j][i] * strz[k-2])) * (u_2[k-2][j][i] - u_2[k][j][i]) + (2 * muz2 + la[k-2][j][i] * strz[k-2] + la[k+1][j][i] * strz[k+1] + 3 * (la[k][j][i] * strz[k] + la[k-1][j][i] * strz[k-1])) * (u_2[k-1][j][i] - u_2[k][j][i]) + (2 * muz3 + la[k-1][j][i] * strz[k-1] + la[k+2][j][i] * strz[k+2] + 3 * (la[k+1][j][i] * strz[k+1] + la[k][j][i] * strz[k])) * (u_2[k+1][j][i] - u_2[k][j][i]) + (2 * muz4 + la[k+1][j][i] * strz[k+1] - 3e0 / 4 * (la[k][j][i] * strz[k] + la[k+2][j][i] * strz[k+2])) * (u_2[k+2][j][i] - u_2[k][j][i])));

					r1 += strx[i] * stry[j] * (1e0 / 144) * (la[k][j][i-2] * (u_1[k][j-2][i-2] - u_1[k][j+2][i-2] + 8 * (-u_1[k][j-1][i-2] + u_1[k][j+1][i-2])) - 8 * (la[k][j][i-1] * (u_1[k][j-2][i-1] - u_1[k][j+2][i-1] + 8 * (-u_1[k][j-1][i-1] + u_1[k][j+1][i-1]))) + 8 * (la[k][j][i+1] * (u_1[k][j-2][i+1] - u_1[k][j+2][i+1] + 8 * (-u_1[k][j-1][i+1] + u_1[k][j+1][i+1]))) - (la[k][j][i+2] * (u_1[k][j-2][i+2] - u_1[k][j+2][i+2] + 8 * (-u_1[k][j-1][i+2] + u_1[k][j+1][i+2])))) + strx[i] * strz[k] * (1e0 / 144) * (la[k][j][i-2] * (u_2[k-2][j][i-2] - u_2[k+2][j][i-2] + 8 * (-u_2[k-1][j][i-2] + u_2[k+1][j][i-2])) - 8 * (la[k][j][i-1] * (u_2[k-2][j][i-1] - u_2[k+2][j][i-1] + 8 * (-u_2[k-1][j][i-1] + u_2[k+1][j][i-1]))) + 8 * (la[k][j][i+1] * (u_2[k-2][j][i+1] - u_2[k+2][j][i+1] + 8 * (-u_2[k-1][j][i+1] + u_2[k+1][j][i+1]))) - (la[k][j][i+2] * (u_2[k-2][j][i+2] - u_2[k+2][j][i+2] + 8 * (-u_2[k-1][j][i+2] + u_2[k+1][j][i+2])))) + strx[i] * stry[j] * (1e0 / 144) * (mu[k][j-2][i] * (u_1[k][j-2][i-2] - u_1[k][j-2][i+2] + 8 * (-u_1[k][j-2][i-1] + u_1[k][j-2][i+1])) - 8 * (mu[k][j-1][i] * (u_1[k][j-1][i-2] - u_1[k][j-1][i+2] + 8 * (-u_1[k][j-1][i-1] + u_1[k][j-1][i+1]))) + 8 * (mu[k][j+1][i] * (u_1[k][j+1][i-2] - u_1[k][j+1][i+2] + 8 * (-u_1[k][j+1][i-1] + u_1[k][j+1][i+1]))) - (mu[k][j+2][i] * (u_1[k][j+2][i-2] - u_1[k][j+2][i+2] + 8 * (-u_1[k][j+2][i-1] + u_1[k][j+2][i+1])))) + strx[i] * strz[k] * (1e0 / 144) * (mu[k-2][j][i] * (u_2[k-2][j][i-2] - u_2[k-2][j][i+2] + 8 * (-u_2[k-2][j][i-1] + u_2[k-2][j][i+1])) - 8 * (mu[k-1][j][i] * (u_2[k-1][j][i-2] - u_2[k-1][j][i+2] + 8 * (-u_2[k-1][j][i-1] + u_2[k-1][j][i+1]))) + 8 * (mu[k+1][j][i] * (u_2[k+1][j][i-2] - u_2[k+1][j][i+2] + 8 * (-u_2[k+1][j][i-1] + u_2[k+1][j][i+1]))) - (mu[k+2][j][i] * (u_2[k+2][j][i-2] - u_2[k+2][j][i+2] + 8 * (-u_2[k+2][j][i-1] + u_2[k+2][j][i+1]))));
					r2 += strx[i] * stry[j] * (1e0 / 144) * (mu[k][j][i-2] * (u_0[k][j-2][i-2] - u_0[k][j+2][i-2] + 8 * (-u_0[k][j-1][i-2] + u_0[k][j+1][i-2])) - 8 * (mu[k][j][i-1] * (u_0[k][j-2][i-1] - u_0[k][j+2][i-1] + 8 * (-u_0[k][j-1][i-1] + u_0[k][j+1][i-1]))) + 8 * (mu[k][j][i+1] * (u_0[k][j-2][i+1] - u_0[k][j+2][i+1] + 8 * (-u_0[k][j-1][i+1] + u_0[k][j+1][i+1]))) - (mu[k][j][i+2] * (u_0[k][j-2][i+2] - u_0[k][j+2][i+2] + 8 * (-u_0[k][j-1][i+2] + u_0[k][j+1][i+2])))) + strx[i] * stry[j] * (1e0 / 144) * (la[k][j-2][i] * (u_0[k][j-2][i-2] - u_0[k][j-2][i+2] + 8 * (-u_0[k][j-2][i-1] + u_0[k][j-2][i+1])) - 8 * (la[k][j-1][i] * (u_0[k][j-1][i-2] - u_0[k][j-1][i+2] + 8 * (-u_0[k][j-1][i-1] + u_0[k][j-1][i+1]))) + 8 * (la[k][j+1][i] * (u_0[k][j+1][i-2] - u_0[k][j+1][i+2] + 8 * (-u_0[k][j+1][i-1] + u_0[k][j+1][i+1]))) - (la[k][j+2][i] * (u_0[k][j+2][i-2] - u_0[k][j+2][i+2] + 8 * (-u_0[k][j+2][i-1] + u_0[k][j+2][i+1])))) + stry[j] * strz[k] * (1e0 / 144) * (la[k][j-2][i] * (u_2[k-2][j-2][i] - u_2[k+2][j-2][i] + 8 * (-u_2[k-1][j-2][i] + u_2[k+1][j-2][i])) - 8 * (la[k][j-1][i] * (u_2[k-2][j-1][i] - u_2[k+2][j-1][i] + 8 * (-u_2[k-1][j-1][i] + u_2[k+1][j-1][i]))) + 8 * (la[k][j+1][i] * (u_2[k-2][j+1][i] - u_2[k+2][j+1][i] + 8 * (-u_2[k-1][j+1][i] + u_2[k+1][j+1][i]))) - (la[k][j+2][i] * (u_2[k-2][j+2][i] - u_2[k+2][j+2][i] + 8 * (-u_2[k-1][j+2][i] + u_2[k+1][j+2][i])))) + stry[j] * strz[k] * (1e0 / 144) * (mu[k-2][j][i] * (u_2[k-2][j-2][i] - u_2[k-2][j+2][i] + 8 * (-u_2[k-2][j-1][i] + u_2[k-2][j+1][i])) - 8 * (mu[k-1][j][i] * (u_2[k-1][j-2][i] - u_2[k-1][j+2][i] + 8 * (-u_2[k-1][j-1][i] + u_2[k-1][j+1][i]))) + 8 * (mu[k+1][j][i] * (u_2[k+1][j-2][i] - u_2[k+1][j+2][i] + 8 * (-u_2[k+1][j-1][i] + u_2[k+1][j+1][i]))) - (mu[k+2][j][i] * (u_2[k+2][j-2][i] - u_2[k+2][j+2][i] + 8 * (-u_2[k+2][j-1][i] + u_2[k+2][j+1][i]))));
					r3 += strx[i] * strz[k] * (1e0 / 144) * (mu[k][j][i-2] * (u_0[k-2][j][i-2] - u_0[k+2][j][i-2] + 8 * (-u_0[k-1][j][i-2] + u_0[k+1][j][i-2])) - 8 * (mu[k][j][i-1] * (u_0[k-2][j][i-1] - u_0[k+2][j][i-1] + 8 * (-u_0[k-1][j][i-1] + u_0[k+1][j][i-1]))) + 8 * (mu[k][j][i+1] * (u_0[k-2][j][i+1] - u_0[k+2][j][i+1] + 8 * (-u_0[k-1][j][i+1] + u_0[k+1][j][i+1]))) - (mu[k][j][i+2] * (u_0[k-2][j][i+2] - u_0[k+2][j][i+2] + 8 * (-u_0[k-1][j][i+2] + u_0[k+1][j][i+2])))) + stry[j] * strz[k] * (1e0 / 144) * (mu[k][j-2][i] * (u_1[k-2][j-2][i] - u_1[k+2][j-2][i] + 8 * (-u_1[k-1][j-2][i] + u_1[k+1][j-2][i])) - 8 * (mu[k][j-1][i] * (u_1[k-2][j-1][i] - u_1[k+2][j-1][i] + 8 * (-u_1[k-1][j-1][i] + u_1[k+1][j-1][i]))) + 8 * (mu[k][j+1][i] * (u_1[k-2][j+1][i] - u_1[k+2][j+1][i] + 8 * (-u_1[k-1][j+1][i] + u_1[k+1][j+1][i]))) - (mu[k][j+2][i] * (u_1[k-2][j+2][i] - u_1[k+2][j+2][i] + 8 * (-u_1[k-1][j+2][i] + u_1[k+1][j+2][i])))) + strx[i] * strz[k] * (1e0 / 144) * (la[k-2][j][i] * (u_0[k-2][j][i-2] - u_0[k-2][j][i+2] + 8 * (-u_0[k-2][j][i-1] + u_0[k-2][j][i+1])) - 8 * (la[k-1][j][i] * (u_0[k-1][j][i-2] - u_0[k-1][j][i+2] + 8 * (-u_0[k-1][j][i-1] + u_0[k-1][j][i+1]))) + 8 * (la[k+1][j][i] * (u_0[k+1][j][i-2] - u_0[k+1][j][i+2] + 8 * (-u_0[k+1][j][i-1] + u_0[k+1][j][i+1]))) - (la[k+2][j][i] * (u_0[k+2][j][i-2] - u_0[k+2][j][i+2] + 8 * (-u_0[k+2][j][i-1] + u_0[k+2][j][i+1])))) + stry[j] * strz[k] * (1e0 / 144) * (la[k-2][j][i] * (u_1[k-2][j-2][i] - u_1[k-2][j+2][i] + 8 * (-u_1[k-2][j-1][i] + u_1[k-2][j+1][i])) - 8 * (la[k-1][j][i] * (u_1[k-1][j-2][i] - u_1[k-1][j+2][i] + 8 * (-u_1[k-1][j-1][i] + u_1[k-1][j+1][i]))) + 8 * (la[k+1][j][i] * (u_1[k+1][j-2][i] - u_1[k+1][j+2][i] + 8 * (-u_1[k+1][j-1][i] + u_1[k+1][j+1][i]))) - (la[k+2][j][i] * (u_1[k+2][j-2][i] - u_1[k+2][j+2][i] + 8 * (-u_1[k+2][j-1][i] + u_1[k+2][j+1][i]))));

					uacc_0[k][j][i] = 0.1 * a1 * uacc_0[k][j][i] + cof * r1;
					uacc_1[k][j][i] = 0.1 * a1 * uacc_1[k][j][i] + cof * r2;
					uacc_2[k][j][i] = 0.1 * a1 * uacc_2[k][j][i] + cof * r3;
				}
			}
		}
	}

	flush_llc ();

	start_time = rtclock ();
	for (t=0; t<2; t++) {
#pragma omp parallel for private (j,i)
		for (k = 2; k < N-2; k++) {
			for (j = 2; j < N-2; j++) {
#pragma GCC ivdep 
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 2; i < N-2; i++) {
					double mux1 = mu[k][j][i-1] * strx[i-1] - 3e0 / 4 * (mu[k][j][i] * strx[i] + mu[k][j][i-2] * strx[i-2]);
					double mux2 = mu[k][j][i-2] * strx[i-2] + mu[k][j][i+1] * strx[i+1] + 3 * (mu[k][j][i] * strx[i] + mu[k][j][i-1] * strx[i-1]);
					double mux3 = mu[k][j][i-1] * strx[i-1] + mu[k][j][i+2] * strx[i+2] + 3 * (mu[k][j][i+1] * strx[i+1] + mu[k][j][i] * strx[i]);
					double mux4 = mu[k][j][i+1] * strx[i+1] - 3e0 / 4 * (mu[k][j][i] * strx[i] + mu[k][j][i+2] * strx[i+2]);
					double muy1 = mu[k][j-1][i] * stry[j-1] - 3e0 / 4 * (mu[k][j][i] * stry[j] + mu[k][j-2][i] * stry[j-2]);
					double muy2 = mu[k][j-2][i] * stry[j-2] + mu[k][j+1][i] * stry[j+1] + 3 * (mu[k][j][i] * stry[j] + mu[k][j-1][i] * stry[j-1]);
					double muy3 = mu[k][j-1][i] * stry[j-1] + mu[k][j+2][i] * stry[j+2] + 3 * (mu[k][j+1][i] * stry[j+1] + mu[k][j][i] * stry[j]);
					double muy4 = mu[k][j+1][i] * stry[j+1] - 3e0 / 4 * (mu[k][j][i] * stry[j] + mu[k][j+2][i] * stry[j+2]);
					double muz1 = mu[k-1][j][i] * strz[k-1] - 3e0 / 4 * (mu[k][j][i] * strz[k] + mu[k-2][j][i] * strz[k-2]);
					double muz2 = mu[k-2][j][i] * strz[k-2] + mu[k+1][j][i] * strz[k+1] + 3 * (mu[k][j][i] * strz[k] + mu[k-1][j][i] * strz[k-1]);
					double muz3 = mu[k-1][j][i] * strz[k-1] + mu[k+2][j][i] * strz[k+2] + 3 * (mu[k+1][j][i] * strz[k+1] + mu[k][j][i] * strz[k]);
					double muz4 = mu[k+1][j][i] * strz[k+1] - 3e0 / 4 * (mu[k][j][i] * strz[k] + mu[k+2][j][i] * strz[k+2]);

					double r1 = 1e0 / 6 * (strx[i] * ((2 * mux1 + la[k][j][i-1] * strx[i-1] - 3e0 / 4 * (la[k][j][i] * strx[i] + la[k][j][i-2] * strx[i-2])) * (u_0[k][j][i-2] - u_0[k][j][i]) + (2 * mux2 + la[k][j][i-2] * strx[i-2] + la[k][j][i+1] * strx[i+1] + 3 * (la[k][j][i] * strx[i] + la[k][j][i-1] * strx[i-1])) * (u_0[k][j][i-1] - u_0[k][j][i]) + (2 * mux3 + la[k][j][i-1] * strx[i-1] + la[k][j][i+2] * strx[i+2] + 3 * (la[k][j][i+1] * strx[i+1] + la[k][j][i] * strx[i])) * (u_0[k][j][i+1] - u_0[k][j][i]) + (2 * mux4 + la[k][j][i+1] * strx[i+1] - 3e0 / 4 * (la[k][j][i] * strx[i] + la[k][j][i+2] * strx[i+2])) * (u_0[k][j][i+2] - u_0[k][j][i])) + stry[j] * (muy1 * (u_0[k][j-2][i] - u_0[k][j][i]) + muy2 * (u_0[k][j-1][i] - u_0[k][j][i]) + muy3 * (u_0[k][j+1][i] - u_0[k][j][i]) + muy4 * (u_0[k][j+2][i] - u_0[k][j][i])) + strz[k] * (muz1 * (u_0[k-2][j][i] - u_0[k][j][i]) + muz2 * (u_0[k-1][j][i] - u_0[k][j][i]) + muz3 * (u_0[k+1][j][i] - u_0[k][j][i]) + muz4 * (u_0[k+2][j][i] - u_0[k][j][i])));
					double r2 = 1e0 / 6 * (strx[i] * (mux1 * (u_1[k][j][i-2] - u_1[k][j][i]) + mux2 * (u_1[k][j][i-1] - u_1[k][j][i]) + mux3 * (u_1[k][j][i+1] - u_1[k][j][i]) + mux4 * (u_1[k][j][i+2] - u_1[k][j][i])) + stry[j] * ((2 * muy1 + la[k][j-1][i] * stry[j-1] - 3e0 / 4 * (la[k][j][i] * stry[j] + la[k][j-2][i] * stry[j-2])) * (u_1[k][j-2][i] - u_1[k][j][i]) + (2 * muy2 + la[k][j-2][i] * stry[j-2] + la[k][j+1][i] * stry[j+1] + 3 * (la[k][j][i] * stry[j] + la[k][j-1][i] * stry[j-1])) * (u_1[k][j-1][i] - u_1[k][j][i]) + (2 * muy3 + la[k][j-1][i] * stry[j-1] + la[k][j+2][i] * stry[j+2] + 3 * (la[k][j+1][i] * stry[j+1] + la[k][j][i] * stry[j])) * (u_1[k][j+1][i] - u_1[k][j][i]) + (2 * muy4 + la[k][j+1][i] * stry[j+1] - 3e0 / 4 * (la[k][j][i] * stry[j] + la[k][j+2][i] * stry[j+2])) * (u_1[k][j+2][i] - u_1[k][j][i])) + strz[k] * (muz1 * (u_1[k-2][j][i] - u_1[k][j][i]) + muz2 * (u_1[k-1][j][i] - u_1[k][j][i]) + muz3 * (u_1[k+1][j][i] - u_1[k][j][i]) + muz4 * (u_1[k+2][j][i] - u_1[k][j][i])));
					double r3 = 1e0 / 6 * (strx[i] * (mux1 * (u_2[k][j][i-2] - u_2[k][j][i]) + mux2 * (u_2[k][j][i-1] - u_2[k][j][i]) + mux3 * (u_2[k][j][i+1] - u_2[k][j][i]) + mux4 * (u_2[k][j][i+2] - u_2[k][j][i])) + stry[j] * (muy1 * (u_2[k][j-2][i] - u_2[k][j][i]) + muy2 * (u_2[k][j-1][i] - u_2[k][j][i]) + muy3 * (u_2[k][j+1][i] - u_2[k][j][i]) + muy4 * (u_2[k][j+2][i] - u_2[k][j][i])) + strz[k] * ((2 * muz1 + la[k-1][j][i] * strz[k-1] - 3e0 / 4 * (la[k][j][i] * strz[k] + la[k-2][j][i] * strz[k-2])) * (u_2[k-2][j][i] - u_2[k][j][i]) + (2 * muz2 + la[k-2][j][i] * strz[k-2] + la[k+1][j][i] * strz[k+1] + 3 * (la[k][j][i] * strz[k] + la[k-1][j][i] * strz[k-1])) * (u_2[k-1][j][i] - u_2[k][j][i]) + (2 * muz3 + la[k-1][j][i] * strz[k-1] + la[k+2][j][i] * strz[k+2] + 3 * (la[k+1][j][i] * strz[k+1] + la[k][j][i] * strz[k])) * (u_2[k+1][j][i] - u_2[k][j][i]) + (2 * muz4 + la[k+1][j][i] * strz[k+1] - 3e0 / 4 * (la[k][j][i] * strz[k] + la[k+2][j][i] * strz[k+2])) * (u_2[k+2][j][i] - u_2[k][j][i])));

					r1 += strx[i] * stry[j] * (1e0 / 144) * (la[k][j][i-2] * (u_1[k][j-2][i-2] - u_1[k][j+2][i-2] + 8 * (-u_1[k][j-1][i-2] + u_1[k][j+1][i-2])) - 8 * (la[k][j][i-1] * (u_1[k][j-2][i-1] - u_1[k][j+2][i-1] + 8 * (-u_1[k][j-1][i-1] + u_1[k][j+1][i-1]))) + 8 * (la[k][j][i+1] * (u_1[k][j-2][i+1] - u_1[k][j+2][i+1] + 8 * (-u_1[k][j-1][i+1] + u_1[k][j+1][i+1]))) - (la[k][j][i+2] * (u_1[k][j-2][i+2] - u_1[k][j+2][i+2] + 8 * (-u_1[k][j-1][i+2] + u_1[k][j+1][i+2])))) + strx[i] * strz[k] * (1e0 / 144) * (la[k][j][i-2] * (u_2[k-2][j][i-2] - u_2[k+2][j][i-2] + 8 * (-u_2[k-1][j][i-2] + u_2[k+1][j][i-2])) - 8 * (la[k][j][i-1] * (u_2[k-2][j][i-1] - u_2[k+2][j][i-1] + 8 * (-u_2[k-1][j][i-1] + u_2[k+1][j][i-1]))) + 8 * (la[k][j][i+1] * (u_2[k-2][j][i+1] - u_2[k+2][j][i+1] + 8 * (-u_2[k-1][j][i+1] + u_2[k+1][j][i+1]))) - (la[k][j][i+2] * (u_2[k-2][j][i+2] - u_2[k+2][j][i+2] + 8 * (-u_2[k-1][j][i+2] + u_2[k+1][j][i+2])))) + strx[i] * stry[j] * (1e0 / 144) * (mu[k][j-2][i] * (u_1[k][j-2][i-2] - u_1[k][j-2][i+2] + 8 * (-u_1[k][j-2][i-1] + u_1[k][j-2][i+1])) - 8 * (mu[k][j-1][i] * (u_1[k][j-1][i-2] - u_1[k][j-1][i+2] + 8 * (-u_1[k][j-1][i-1] + u_1[k][j-1][i+1]))) + 8 * (mu[k][j+1][i] * (u_1[k][j+1][i-2] - u_1[k][j+1][i+2] + 8 * (-u_1[k][j+1][i-1] + u_1[k][j+1][i+1]))) - (mu[k][j+2][i] * (u_1[k][j+2][i-2] - u_1[k][j+2][i+2] + 8 * (-u_1[k][j+2][i-1] + u_1[k][j+2][i+1])))) + strx[i] * strz[k] * (1e0 / 144) * (mu[k-2][j][i] * (u_2[k-2][j][i-2] - u_2[k-2][j][i+2] + 8 * (-u_2[k-2][j][i-1] + u_2[k-2][j][i+1])) - 8 * (mu[k-1][j][i] * (u_2[k-1][j][i-2] - u_2[k-1][j][i+2] + 8 * (-u_2[k-1][j][i-1] + u_2[k-1][j][i+1]))) + 8 * (mu[k+1][j][i] * (u_2[k+1][j][i-2] - u_2[k+1][j][i+2] + 8 * (-u_2[k+1][j][i-1] + u_2[k+1][j][i+1]))) - (mu[k+2][j][i] * (u_2[k+2][j][i-2] - u_2[k+2][j][i+2] + 8 * (-u_2[k+2][j][i-1] + u_2[k+2][j][i+1]))));
					r2 += strx[i] * stry[j] * (1e0 / 144) * (mu[k][j][i-2] * (u_0[k][j-2][i-2] - u_0[k][j+2][i-2] + 8 * (-u_0[k][j-1][i-2] + u_0[k][j+1][i-2])) - 8 * (mu[k][j][i-1] * (u_0[k][j-2][i-1] - u_0[k][j+2][i-1] + 8 * (-u_0[k][j-1][i-1] + u_0[k][j+1][i-1]))) + 8 * (mu[k][j][i+1] * (u_0[k][j-2][i+1] - u_0[k][j+2][i+1] + 8 * (-u_0[k][j-1][i+1] + u_0[k][j+1][i+1]))) - (mu[k][j][i+2] * (u_0[k][j-2][i+2] - u_0[k][j+2][i+2] + 8 * (-u_0[k][j-1][i+2] + u_0[k][j+1][i+2])))) + strx[i] * stry[j] * (1e0 / 144) * (la[k][j-2][i] * (u_0[k][j-2][i-2] - u_0[k][j-2][i+2] + 8 * (-u_0[k][j-2][i-1] + u_0[k][j-2][i+1])) - 8 * (la[k][j-1][i] * (u_0[k][j-1][i-2] - u_0[k][j-1][i+2] + 8 * (-u_0[k][j-1][i-1] + u_0[k][j-1][i+1]))) + 8 * (la[k][j+1][i] * (u_0[k][j+1][i-2] - u_0[k][j+1][i+2] + 8 * (-u_0[k][j+1][i-1] + u_0[k][j+1][i+1]))) - (la[k][j+2][i] * (u_0[k][j+2][i-2] - u_0[k][j+2][i+2] + 8 * (-u_0[k][j+2][i-1] + u_0[k][j+2][i+1])))) + stry[j] * strz[k] * (1e0 / 144) * (la[k][j-2][i] * (u_2[k-2][j-2][i] - u_2[k+2][j-2][i] + 8 * (-u_2[k-1][j-2][i] + u_2[k+1][j-2][i])) - 8 * (la[k][j-1][i] * (u_2[k-2][j-1][i] - u_2[k+2][j-1][i] + 8 * (-u_2[k-1][j-1][i] + u_2[k+1][j-1][i]))) + 8 * (la[k][j+1][i] * (u_2[k-2][j+1][i] - u_2[k+2][j+1][i] + 8 * (-u_2[k-1][j+1][i] + u_2[k+1][j+1][i]))) - (la[k][j+2][i] * (u_2[k-2][j+2][i] - u_2[k+2][j+2][i] + 8 * (-u_2[k-1][j+2][i] + u_2[k+1][j+2][i])))) + stry[j] * strz[k] * (1e0 / 144) * (mu[k-2][j][i] * (u_2[k-2][j-2][i] - u_2[k-2][j+2][i] + 8 * (-u_2[k-2][j-1][i] + u_2[k-2][j+1][i])) - 8 * (mu[k-1][j][i] * (u_2[k-1][j-2][i] - u_2[k-1][j+2][i] + 8 * (-u_2[k-1][j-1][i] + u_2[k-1][j+1][i]))) + 8 * (mu[k+1][j][i] * (u_2[k+1][j-2][i] - u_2[k+1][j+2][i] + 8 * (-u_2[k+1][j-1][i] + u_2[k+1][j+1][i]))) - (mu[k+2][j][i] * (u_2[k+2][j-2][i] - u_2[k+2][j+2][i] + 8 * (-u_2[k+2][j-1][i] + u_2[k+2][j+1][i]))));
					r3 += strx[i] * strz[k] * (1e0 / 144) * (mu[k][j][i-2] * (u_0[k-2][j][i-2] - u_0[k+2][j][i-2] + 8 * (-u_0[k-1][j][i-2] + u_0[k+1][j][i-2])) - 8 * (mu[k][j][i-1] * (u_0[k-2][j][i-1] - u_0[k+2][j][i-1] + 8 * (-u_0[k-1][j][i-1] + u_0[k+1][j][i-1]))) + 8 * (mu[k][j][i+1] * (u_0[k-2][j][i+1] - u_0[k+2][j][i+1] + 8 * (-u_0[k-1][j][i+1] + u_0[k+1][j][i+1]))) - (mu[k][j][i+2] * (u_0[k-2][j][i+2] - u_0[k+2][j][i+2] + 8 * (-u_0[k-1][j][i+2] + u_0[k+1][j][i+2])))) + stry[j] * strz[k] * (1e0 / 144) * (mu[k][j-2][i] * (u_1[k-2][j-2][i] - u_1[k+2][j-2][i] + 8 * (-u_1[k-1][j-2][i] + u_1[k+1][j-2][i])) - 8 * (mu[k][j-1][i] * (u_1[k-2][j-1][i] - u_1[k+2][j-1][i] + 8 * (-u_1[k-1][j-1][i] + u_1[k+1][j-1][i]))) + 8 * (mu[k][j+1][i] * (u_1[k-2][j+1][i] - u_1[k+2][j+1][i] + 8 * (-u_1[k-1][j+1][i] + u_1[k+1][j+1][i]))) - (mu[k][j+2][i] * (u_1[k-2][j+2][i] - u_1[k+2][j+2][i] + 8 * (-u_1[k-1][j+2][i] + u_1[k+1][j+2][i])))) + strx[i] * strz[k] * (1e0 / 144) * (la[k-2][j][i] * (u_0[k-2][j][i-2] - u_0[k-2][j][i+2] + 8 * (-u_0[k-2][j][i-1] + u_0[k-2][j][i+1])) - 8 * (la[k-1][j][i] * (u_0[k-1][j][i-2] - u_0[k-1][j][i+2] + 8 * (-u_0[k-1][j][i-1] + u_0[k-1][j][i+1]))) + 8 * (la[k+1][j][i] * (u_0[k+1][j][i-2] - u_0[k+1][j][i+2] + 8 * (-u_0[k+1][j][i-1] + u_0[k+1][j][i+1]))) - (la[k+2][j][i] * (u_0[k+2][j][i-2] - u_0[k+2][j][i+2] + 8 * (-u_0[k+2][j][i-1] + u_0[k+2][j][i+1])))) + stry[j] * strz[k] * (1e0 / 144) * (la[k-2][j][i] * (u_1[k-2][j-2][i] - u_1[k-2][j+2][i] + 8 * (-u_1[k-2][j-1][i] + u_1[k-2][j+1][i])) - 8 * (la[k-1][j][i] * (u_1[k-1][j-2][i] - u_1[k-1][j+2][i] + 8 * (-u_1[k-1][j-1][i] + u_1[k-1][j+1][i]))) + 8 * (la[k+1][j][i] * (u_1[k+1][j-2][i] - u_1[k+1][j+2][i] + 8 * (-u_1[k+1][j-1][i] + u_1[k+1][j+1][i]))) - (la[k+2][j][i] * (u_1[k+2][j-2][i] - u_1[k+2][j+2][i] + 8 * (-u_1[k+2][j-1][i] + u_1[k+2][j+1][i]))));

					uacc_0[k][j][i] = a1 * uacc_0[k][j][i] + cof * r1;
					uacc_1[k][j][i] = a1 * uacc_1[k][j][i] + cof * r2;
					uacc_2[k][j][i] = a1 * uacc_2[k][j][i] + cof * r3;
				}
			}
		}
	}
	end_time = rtclock ();
        printf ("orig: %6lf\n", (double)320*320*320*687*2/(end_time - start_time)/1e9);
	flush_llc ();
	host_code ((double*)uacc_opt_0, (double*)uacc_opt_1, (double*)uacc_opt_2, (double*)u_0, (double*)u_1, (double*)u_2, (double*)mu, (double*)la, (double*)strx, (double*)stry, (double*)strz, N);

	double error_0 = checkError3D (N, N, 0, (double*)uacc_0, (double*)uacc_opt_0, 2, N-2, 2, N-2, 2, N-2);
	if (error_0 > TOLERANCE)
		printf ("error %e\n", error_0);
	double error_1 = checkError3D (N, N, 0, (double*)uacc_1, (double*)uacc_opt_1, 2, N-2, 2, N-2, 2, N-2);
	if (error_1 > TOLERANCE)
		printf ("error %e\n", error_1);
	double error_2 = checkError3D (N, N, 0, (double*)uacc_2, (double*)uacc_opt_2, 2, N-2, 2, N-2, 2, N-2);
	if (error_2 > TOLERANCE)
		printf ("error %e\n", error_2);

}
