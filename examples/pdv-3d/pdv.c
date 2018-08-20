#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void pdv_opt (double* t_density1, double* t_energy1, double* t_density0, double* t_energy0, double* t_volume, double *t_pressure, double *t_viscosity, double* t_xarea, double* t_yarea, double* t_zarea, double* t_xvel0, double* t_xvel1, double* t_yvel0, double* t_yvel1, double* t_zvel0, double* t_zvel1, int N); 

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

	double (*density1_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*energy1_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*density1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*energy1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*xarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*volume)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*pressure)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*viscosity)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*density0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*energy0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xvel1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
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
						energy1[k][j][i]=energy0[k][j][i]-energy_change;
						density1[k][j][i]=density0[k][j][i]*vol_change_s;

					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*98*5/(end_time - start_time)/1e9);

	pdv_opt ((double*)density1_opt, (double*)energy1_opt, (double*)density0, (double*)energy0, (double*)volume, (double*)pressure, (double*)viscosity, (double*)xarea, (double*)yarea, (double*)zarea, (double*)xvel0, (double*)xvel1, (double*)yvel0, (double*)yvel1, (double*)zvel0, (double*)zvel1, N);

	double error = checkError3D (N, N, 0, (double*)density1_opt, (double*)density1, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n", error);
	error = checkError3D (N, N, 0, (double*)energy1_opt, (double*)energy1, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n", error);
}
