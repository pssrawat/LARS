#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void pdv_opt (double* t_density1, double* t_energy1, double* t_density0, double* t_energy0, double* t_volume, double *t_pressure, double *t_viscosity, double* t_xarea, double* t_yarea, double* t_xvel0, double* t_xvel1, double* t_yvel0, double* t_yvel1, int N); 

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

	double (*density1_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*energy1_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*density1)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*energy1)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*xarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*volume)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*pressure)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*viscosity)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*density0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*energy0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xvel1)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yvel1)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
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
					energy1[j][i]=energy0[j][i]-energy_change;
					density1[j][i]=density0[j][i]*volume_change_s;
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)4096*4096*51*10/(end_time - start_time)/1e9);

	pdv_opt ((double*)density1_opt, (double*)energy1_opt, (double*)density0, (double*)energy0, (double*)volume, (double*)pressure, (double*)viscosity, (double*)xarea, (double*)yarea,  (double*)xvel0, (double*)xvel1, (double*)yvel0, (double*)yvel1, N);

	double error = checkError2D (N, 0, (double*)density1_opt, (double*)density1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
	printf("error %e\n", error);
	error = checkError2D (N, 0, (double*)energy1_opt, (double*)energy1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
	printf("error %e\n", error);
}
