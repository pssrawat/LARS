#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void accelerate_opt (double *t_xvel1, double *t_yvel1, double *t_zvel1, double *t_xarea, double *t_yarea, double *t_zarea, double *t_xvel0, double *t_yvel0, double *t_zvel0, double *t_pressure, double *t_viscosity, double *t_density0, double *t_volume, double dt, int N) {
	double (*xvel1)[4098] = (double (*)[4098]) t_xvel1;
	double (*yvel1)[4098] = (double (*)[4098]) t_yvel1;
	double (*zvel1)[4098] = (double (*)[4098]) t_zvel1;
	double (*xarea)[4098] = (double (*)[4098]) t_xarea;
	double (*yarea)[4098] = (double (*)[4098]) t_yarea;
	double (*zarea)[4098] = (double (*)[4098]) t_zarea;
	double (*xvel0)[4098] = (double (*)[4098]) t_xvel0;
	double (*yvel0)[4098] = (double (*)[4098]) t_yvel0;
	double (*zvel0)[4098] = (double (*)[4098]) t_zvel0;
	double (*pressure)[4098] = (double (*)[4098]) t_pressure;
	double (*viscosity)[4098] = (double (*)[4098]) t_viscosity;
	double (*density0)[4098] = (double (*)[4098]) t_density0;
	double (*volume)[4098] = (double (*)[4098]) t_volume;

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
					double stepbymass_s = 0.5*dt /((density0[j-1][i-1]*volume[j-1][i-1] +density0[j-1][i]*volume[j-1][i] +density0[j][i]*volume[j][i] +density0[j][i-1]*volume[j][i-1]) *0.25);
					xvel1[j][i]=0.1*xvel0[j][i]-stepbymass_s*(xarea[j][i]*(pressure[j][i]-pressure[j][i-1]) +xarea[j-1][i]*(pressure[j-1][i]-pressure[j-1][i-1]));
					yvel1[j][i]=0.1*yvel0[j][i]-stepbymass_s*(yarea[j][i]*(pressure[j][i]-pressure[j-1][i]) +yarea[j][i-1]*(pressure[j][i-1]-pressure[j-1][i-1]));
					xvel1[j][i]=0.1*xvel1[j][i]-stepbymass_s*(xarea[j][i]*(viscosity[j][i]-viscosity[j][i-1]) +xarea[j-1][i]*(viscosity[j-1][i]-viscosity[j-1][i-1]));
					yvel1[j][i]=0.1*yvel1[j][i]-stepbymass_s*(yarea[j][i]*(viscosity[j][i]-viscosity[j-1][i]) +yarea[j][i-1]*(viscosity[j][i-1]-viscosity[j-1][i-1]));
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j+=1) {
				for (i = 1; i < N-1; i+=4) {
#pragma begin stencil1 unroll j=1,i=1 print-intrinsics true acc-size 1
					stepbymass_s = 0.5*dt /((density0[j-1][i-1]*volume[j-1][i-1] +density0[j-1][i]*volume[j-1][i] +density0[j][i]*volume[j][i] +density0[j][i-1]*volume[j][i-1]) *0.25);
					xvel1[j][i]=xvel0[j][i]-stepbymass_s*(xarea[j][i]*(pressure[j][i]-pressure[j][i-1]) +xarea[j-1][i]*(pressure[j-1][i]-pressure[j-1][i-1]));
					yvel1[j][i]=yvel0[j][i]-stepbymass_s*(yarea[j][i]*(pressure[j][i]-pressure[j-1][i]) +yarea[j][i-1]*(pressure[j][i-1]-pressure[j-1][i-1]));
					xvel1[j][i]=xvel1[j][i]-stepbymass_s*(xarea[j][i]*(viscosity[j][i]-viscosity[j][i-1]) +xarea[j-1][i]*(viscosity[j-1][i]-viscosity[j-1][i-1]));
					yvel1[j][i]=yvel1[j][i]-stepbymass_s*(yarea[j][i]*(viscosity[j][i]-viscosity[j-1][i]) +yarea[j][i-1]*(viscosity[j][i-1]-viscosity[j-1][i-1]));
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)4096*4096*38*5/(end_time - start_time)/1e9);
}
