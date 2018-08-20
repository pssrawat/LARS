#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void pdv_opt (double* t_density1, double* t_energy1, double* t_density0, double* t_energy0, double* t_volume, double *t_pressure, double *t_viscosity, double* t_xarea, double* t_yarea, double* t_zarea, double* t_xvel0, double* t_xvel1, double* t_yvel0, double* t_yvel1, double* t_zvel0, double* t_zvel1, int N) {
	double (*density1)[258][258] = (double (*)[258][258]) t_density1;
	double (*energy1)[258][258] = (double (*)[258][258]) t_energy1;
	double (*density0)[258][258] = (double (*)[258][258]) t_density0;
	double (*energy0)[258][258] = (double (*)[258][258]) t_energy0;
	double (*volume)[258][258] = (double (*)[258][258]) t_volume;
	double (*pressure)[258][258] = (double (*)[258][258]) t_pressure;
	double (*viscosity)[258][258] = (double (*)[258][258]) t_viscosity;
	double (*xarea)[258][258] = (double (*)[258][258]) t_xarea;
	double (*yarea)[258][258] = (double (*)[258][258]) t_yarea;
	double (*zarea)[258][258] = (double (*)[258][258]) t_zarea;
	double (*xvel0)[258][258] = (double (*)[258][258]) t_xvel0;
	double (*xvel1)[258][258] = (double (*)[258][258]) t_xvel1;
	double (*yvel0)[258][258] = (double (*)[258][258]) t_yvel0;
	double (*yvel1)[258][258] = (double (*)[258][258]) t_yvel1;
	double (*zvel0)[258][258] = (double (*)[258][258]) t_zvel0;
	double (*zvel1)[258][258] = (double (*)[258][258]) t_zvel1;
	double dt = 0.013;

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
						double left_flux =  (xarea[k][j][i]*(xvel0[k][j][i]+xvel0[k][j+1][i]+xvel0[k+1][j][i]+xvel0[k+1][j+1][i] +xvel0[k][j][i]+xvel0[k][j+1][i]+xvel0[k+1][j][i]+xvel0[k+1][j+1][i])) *0.125*dt*0.5;
						double right_flux= (xarea[k][j][i+1]*(xvel0[k][j][i+1]+xvel0[k][j+1][i+1]+xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1] +xvel0[k][j][i+1]+xvel0[k][j+1][i+1]+xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1])) *0.125*dt*0.5;
						double bottom_flux=(yarea[k][j][i]*(yvel0[k][j][i]+yvel0[k][j][i+1]+yvel0[k+1][j][i]+yvel0[k+1][j][i+1]  +yvel0[k][j][i]+yvel0[k][j][i+1]+yvel0[k+1][j][i]+yvel0[k+1][j][i+1])) *0.125*dt*0.5;
						double top_flux=   (yarea[k][j+1][i]*(yvel0[k][j+1][i]+yvel0[k][j+1][i+1]+yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1] + yvel0[k][j+1][i]+yvel0[k][j+1][i+1]+yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1])) *0.125*dt*0.5;
						double back_flux=  (zarea[k][j][i]*(zvel1[k][j][i]+zvel0[k][j][i+1]+zvel0[k][j+1][i]+zvel0[k][j+1][i+1]  +zvel0[k][j][i]+zvel0[k][j][i+1]+zvel0[k][j+1][i]+zvel0[k][j+1][i+1])) *0.125*dt*0.5;
						double front_flux= (zarea[k+1][j][i]*(zvel0[k+1][j][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i]+zvel0[k+1][j+1][i+1]  +zvel0[k+1][j][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i]+zvel0[k+1][j+1][i+1])) *0.125*dt*0.5;
						double total_flux=right_flux-left_flux+top_flux-bottom_flux+front_flux-back_flux;
						double vol_change_s=volume[k][j][i]/(volume[k][j][i]+total_flux);
						double min_cell_volume=(volume[k][j][i]+right_flux-left_flux+top_flux-bottom_flux+front_flux-back_flux + volume[k][j][i]+right_flux-left_flux+top_flux-bottom_flux + volume[k][j][i]+right_flux-left_flux + volume[k][j][i]+top_flux-bottom_flux);
						double recip_volume=1.0/volume[k][j][i];
						double energy_change=(pressure[k][j][i]/density0[k][j][i]+viscosity[k][j][i]/density0[k][j][i])*total_flux*recip_volume;
						energy1[k][j][i]=0.1*energy0[k][j][i]-energy_change;
						density1[k][j][i]=0.1*density0[k][j][i]*vol_change_s;
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
			for (k = 1; k < N-1; k+=2) {
				for (j = 1; j < N-1; j++) {
					for (i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll k=2,j=1,i=1 print-intrinsics true acc-size 1
						left_flux =  (xarea[k][j][i]*(xvel0[k][j][i]+xvel0[k][j+1][i]+xvel0[k+1][j][i]+xvel0[k+1][j+1][i] +xvel0[k][j][i]+xvel0[k][j+1][i]+xvel0[k+1][j][i]+xvel0[k+1][j+1][i])) *0.125*dt*0.5;
						right_flux= (xarea[k][j][i+1]*(xvel0[k][j][i+1]+xvel0[k][j+1][i+1]+xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1] +xvel0[k][j][i+1]+xvel0[k][j+1][i+1]+xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1])) *0.125*dt*0.5;
						bottom_flux=(yarea[k][j][i]*(yvel0[k][j][i]+yvel0[k][j][i+1]+yvel0[k+1][j][i]+yvel0[k+1][j][i+1]  +yvel0[k][j][i]+yvel0[k][j][i+1]+yvel0[k+1][j][i]+yvel0[k+1][j][i+1])) *0.125*dt*0.5;
						top_flux=   (yarea[k][j+1][i]*(yvel0[k][j+1][i]+yvel0[k][j+1][i+1]+yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1] + yvel0[k][j+1][i]+yvel0[k][j+1][i+1]+yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1])) *0.125*dt*0.5;
						back_flux=  (zarea[k][j][i]*(zvel1[k][j][i]+zvel0[k][j][i+1]+zvel0[k][j+1][i]+zvel0[k][j+1][i+1]  +zvel0[k][j][i]+zvel0[k][j][i+1]+zvel0[k][j+1][i]+zvel0[k][j+1][i+1])) *0.125*dt*0.5;
						front_flux= (zarea[k+1][j][i]*(zvel0[k+1][j][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i]+zvel0[k+1][j+1][i+1]  +zvel0[k+1][j][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i]+zvel0[k+1][j+1][i+1])) *0.125*dt*0.5;
						total_flux=right_flux-left_flux+top_flux-bottom_flux+front_flux-back_flux;
						vol_change_s=volume[k][j][i]/(volume[k][j][i]+total_flux);
						min_cell_volume=(volume[k][j][i]+right_flux-left_flux+top_flux-bottom_flux+front_flux-back_flux + volume[k][j][i]+right_flux-left_flux+top_flux-bottom_flux + volume[k][j][i]+right_flux-left_flux + volume[k][j][i]+top_flux-bottom_flux);
						recip_volume=1.0/volume[k][j][i];
						energy_change=(pressure[k][j][i]/density0[k][j][i]+viscosity[k][j][i]/density0[k][j][i])*total_flux*recip_volume;
						energy1[k][j][i]=energy0[k][j][i]-energy_change;
						density1[k][j][i]=density0[k][j][i]*vol_change_s;
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)256*256*256*98*5/(end_time - start_time)/1e9);
}
