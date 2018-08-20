#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern double rtclock (void);

void ideal_gas_opt (double* t_soundspeed, double* t_pressure, double* t_density, double* t_energy, int N) { 
	double (*soundspeed)[4098] = (double (*)[4098]) t_soundspeed;
	double (*pressure)[4098] = (double (*)[4098]) t_pressure;
	double (*density)[4098] = (double (*)[4098]) t_density;
	double (*energy)[4098] = (double (*)[4098]) t_energy;

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
			for (j = 1; j < N-1; j+=4) {
				for (i = 1; i < N-1; i++) {
#pragma begin stencil1 unroll j=4,i=1 print-intrinsics false acc-size 1
					v=1.0/density[j][i];
					pressure[j][i]=(1.4-1.0)*density[j][i]*energy[j][i];
					pressurebyenergy=(1.4-1.0)*density[j][i];
					pressurebyvolume=-density[j][i]*pressure[j][i];
					sound_speed_squared=v*v*(pressure[j][i]*pressurebyenergy-pressurebyvolume);
					soundspeed[j][i]=sqrt(sound_speed_squared);
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)4096*4096*12*10/(end_time - start_time)/1e9);
}
