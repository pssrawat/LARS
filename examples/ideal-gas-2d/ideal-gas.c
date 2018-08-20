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
	int N = 4098;
	double (*soundspeed_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*soundspeed)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*pressure)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*density)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*energy)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
        			double v=1.0/density[j][i];
        			pressure[j][i]=(1.4-1.0)*density[j][i]*energy[j][i];
        			double pressurebyenergy=(1.4-1.0)*density[j][i];
        			double pressurebyvolume=-density[j][i]*pressure[j][i];
        			double sound_speed_squared=v*v*(pressure[j][i]*pressurebyenergy-pressurebyvolume);
        			soundspeed[j][i]=0.1*sqrt(sound_speed_squared);
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<10; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
        			double v=1.0/density[j][i];
        			pressure[j][i]=(1.4-1.0)*density[j][i]*energy[j][i];
        			double pressurebyenergy=(1.4-1.0)*density[j][i];
        			double pressurebyvolume=-density[j][i]*pressure[j][i];
        			double sound_speed_squared=v*v*(pressure[j][i]*pressurebyenergy-pressurebyvolume);
        			soundspeed[j][i]=sqrt(sound_speed_squared);
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)4096*4096*12*10/(end_time - start_time)/1e9);

	ideal_gas_opt ((double*)soundspeed_opt, (double*)pressure, (double*)density, (double*)energy, N);

	double error = checkError2D (N, 0, (double*)soundspeed_opt, (double*)soundspeed, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n", error);
}
