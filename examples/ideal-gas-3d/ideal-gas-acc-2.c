#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern double rtclock (void);

void ideal_gas_opt (double* t_soundspeed, double* t_pressure, double* t_density, double* t_energy, int N) { 
	double (*soundspeed)[258][258] = (double (*)[258][258]) t_soundspeed;
	double (*pressure)[258][258] = (double (*)[258][258]) t_pressure;
	double (*density)[258][258] = (double (*)[258][258]) t_density;
	double (*energy)[258][258] = (double (*)[258][258]) t_energy;

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
			for (k = 1; k < N-1; k+=2) {
				for (j = 1; j < N-1; j++) {
					for (i = 1; i < N-1; i+=1) {
#pragma begin stencil1 unroll k=2,j=1,i=1 print-intrinsics false acc-size 1
						v=1.0/density[k][j][i];
						pressure[k][j][i]=(1.4-1.0)*density[k][j][i]*energy[k][j][i];
						pressurebyenergy=(1.4-1.0)*density[k][j][i];
						pressurebyvolume=-density[k][j][i]*pressure[k][j][i];
						sound_speed_squared=v*v*(pressure[k][j][i]*pressurebyenergy-pressurebyvolume);
						soundspeed[k][j][i]=sqrt(sound_speed_squared);
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)256*256*256*12*10/(end_time - start_time)/1e9);
}
