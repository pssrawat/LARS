#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void ideal_gas_opt (double* soundspeed, double* pressure, double* density, double* energy, int N);

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
	int N = 258;
	double (*soundspeed_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*soundspeed)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*pressure)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*density)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*energy)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);

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
						double v=1.0/density[k][j][i];
						pressure[k][j][i]=(1.4-1.0)*density[k][j][i]*energy[k][j][i];
						double pressurebyenergy=(1.4-1.0)*density[k][j][i];
						double pressurebyvolume=-density[k][j][i]*pressure[k][j][i];
						double sound_speed_squared=v*v*(pressure[k][j][i]*pressurebyenergy-pressurebyvolume);
						soundspeed[k][j][i]=0.1*sqrt(sound_speed_squared);
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
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						double v=1.0/density[k][j][i];
						pressure[k][j][i]=(1.4-1.0)*density[k][j][i]*energy[k][j][i];
						double pressurebyenergy=(1.4-1.0)*density[k][j][i];
						double pressurebyvolume=-density[k][j][i]*pressure[k][j][i];
						double sound_speed_squared=v*v*(pressure[k][j][i]*pressurebyenergy-pressurebyvolume);
						soundspeed[k][j][i]=sqrt(sound_speed_squared);
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*12*10/(end_time - start_time)/1e9);

	ideal_gas_opt ((double*)soundspeed_opt, (double*)pressure, (double*)density, (double*)energy, N);

	double error = checkError3D (N, N, 0, (double*)soundspeed_opt, (double*)soundspeed, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e, ",error);
}
