#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void pdv_opt (double* t_density1, double* t_energy1, double* t_density0, double* t_energy0, double* t_volume, double *t_pressure, double *t_viscosity, double* t_xarea, double* t_yarea, double* t_xvel0, double* t_xvel1, double* t_yvel0, double* t_yvel1, int N) {
	double (*density1)[4098] = (double (*)[4098]) t_density1;
	double (*energy1)[4098] = (double (*)[4098]) t_energy1;
	double (*density0)[4098] = (double (*)[4098]) t_density0;
	double (*energy0)[4098] = (double (*)[4098]) t_energy0;
	double (*volume)[4098] = (double (*)[4098]) t_volume;
	double (*pressure)[4098] = (double (*)[4098]) t_pressure;
	double (*viscosity)[4098] = (double (*)[4098]) t_viscosity;
	double (*xarea)[4098] = (double (*)[4098]) t_xarea;
	double (*yarea)[4098] = (double (*)[4098]) t_yarea;
	double (*xvel0)[4098] = (double (*)[4098]) t_xvel0;
	double (*xvel1)[4098] = (double (*)[4098]) t_xvel1;
	double (*yvel0)[4098] = (double (*)[4098]) t_yvel0;
	double (*yvel1)[4098] = (double (*)[4098]) t_yvel1;
	double dt = 0.013;

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
					double left_flux =   (xarea[j][i]*(xvel0[j][i]+xvel0[j+1][i] +xvel0[j][i]+xvel0[j+1][i]))*0.25*dt*0.5;
					double right_flux =  (xarea[j][i+1]*(xvel0[j][i+1]+xvel0[j+1][i+1]+xvel0[j][i+1]+xvel0[j+1][i+1]))*0.25*dt*0.5;
					double bottom_flux = (yarea[j][i]*(yvel0[j][i]+yvel0[j][i+1] +yvel0[j][i]+yvel0[j][i+1]))*0.25*dt*0.5;
					double top_flux =    (yarea[j+1][i]*(yvel0[j+1][i]+yvel0[j+1][i+1]+yvel0[j+1][i]+yvel0[j+1][i+1]))*0.25*dt*0.5;
					double total_flux = right_flux-left_flux+top_flux-bottom_flux;
					double volume_change_s = volume[j][i]/(volume[j][i]+total_flux);
					double min_cell_volume= volume[j][i]+right_flux-left_flux+top_flux-bottom_flux + volume[j][i]+right_flux-left_flux + volume[j][i]+top_flux-bottom_flux;
					double recip_volume=1.0/volume[j][i];
					double energy_change=(pressure[j][i]/density0[j][i]+viscosity[j][i]/density0[j][i])*total_flux*recip_volume;
					energy1[j][i]=0.1*energy0[j][i]-energy_change;
					density1[j][i]=0.1*density0[j][i]*volume_change_s;
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
				for (i = 1; i < N-1; i+=4) {
#pragma begin stencil1 unroll j=4,i=1 print-intrinsics true acc-size 1
					left_flux =   (xarea[j][i]*(xvel0[j][i]+xvel0[j+1][i] +xvel0[j][i]+xvel0[j+1][i]))*0.25*dt*0.5;
					right_flux =  (xarea[j][i+1]*(xvel0[j][i+1]+xvel0[j+1][i+1]+xvel0[j][i+1]+xvel0[j+1][i+1]))*0.25*dt*0.5;
					bottom_flux = (yarea[j][i]*(yvel0[j][i]+yvel0[j][i+1] +yvel0[j][i]+yvel0[j][i+1]))*0.25*dt*0.5;
					top_flux =    (yarea[j+1][i]*(yvel0[j+1][i]+yvel0[j+1][i+1]+yvel0[j+1][i]+yvel0[j+1][i+1]))*0.25*dt*0.5;
					total_flux = right_flux-left_flux+top_flux-bottom_flux;
					volume_change_s = volume[j][i]/(volume[j][i]+total_flux);
					min_cell_volume= volume[j][i]+right_flux-left_flux+top_flux-bottom_flux + volume[j][i]+right_flux-left_flux + volume[j][i]+top_flux-bottom_flux;
					recip_volume=1.0/volume[j][i];
					energy_change=(pressure[j][i]/density0[j][i]+viscosity[j][i]/density0[j][i])*total_flux*recip_volume;
					energy1[j][i]=energy0[j][i]-energy_change;
					density1[j][i]=density0[j][i]*volume_change_s;
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)4096*4096*51*10/(end_time - start_time)/1e9);
}
