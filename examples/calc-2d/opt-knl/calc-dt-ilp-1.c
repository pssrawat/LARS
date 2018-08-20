#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern double rtclock (void);

void calc_dt_opt (double *t_dt_min, double *t_volume, double *t_soundspeed, double *t_density0, double *t_viscosity_a, double *t_xvel0, double *t_yvel0, double *t_zvel0, double *t_xarea, double *t_yarea, double *t_zarea, double *celldx, double *celldy, double *celldz, double dtc_safe, double dtu_safe, double dtv_safe, double dtw_safe, double dtdiv_safe, int N) {
	double (*dt_min)[4098] = (double (*)[4098]) t_dt_min;
	double (*volume)[4098] = (double (*)[4098]) t_volume;
	double (*soundspeed)[4098] = (double (*)[4098]) t_soundspeed;
	double (*density0)[4098] = (double (*)[4098]) t_density0;
	double (*viscosity_a)[4098] = (double (*)[4098]) t_viscosity_a;
	double (*xvel0)[4098] = (double (*)[4098]) t_xvel0;
	double (*yvel0)[4098] = (double (*)[4098]) t_yvel0;
	double (*zvel0)[4098] = (double (*)[4098]) t_zvel0;
	double (*xarea)[4098] = (double (*)[4098]) t_xarea;
	double (*yarea)[4098] = (double (*)[4098]) t_yarea;
	double (*zarea)[4098] = (double (*)[4098]) t_zarea;

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
					double cc=soundspeed[j][i]*soundspeed[j][i];
					cc=cc+2.0*viscosity_a[j][i]/density0[j][i];
					double dtct=dtc_safe*celldx[i]*celldy[j]/sqrt(cc);
					double dv1=(xvel0[j][i]+xvel0[j+1][i])*xarea[j][i];
					double dv2=(xvel0[j][i+1]+xvel0[j+1][i+1])*xarea[j][i+1];
					double div=dv2-dv1;
					double dtut=dtu_safe*2.0*volume[j][i]/(dv1*dv2*volume[j][i]);

					dv1=(yvel0[j][i]+yvel0[j][i+1])*yarea[j][i];
					dv2=(yvel0[j+1][i]+yvel0[j+1][i+1])*yarea[j+1][i];
					div=div+dv2-dv1;
					double dtvt=dtv_safe*2.0*volume[j][i]/(dv1*dv2*volume[j][i]);
					div=div/(2.0*volume[j][i]);
					double dtdivt=dtdiv_safe*(-1.0/div);
					dt_min[j][i]=0.1*dtct*dtut*dtvt*dtdivt;
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
				for (i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll j=1,i=1 print-intrinsics true acc-size 1
					cc=soundspeed[j][i]*soundspeed[j][i];
					cc=cc+2.0*viscosity_a[j][i]/density0[j][i];
					dtct=dtc_safe*celldx[i]*celldy[j]/sqrt(cc);
					dv1=(xvel0[j][i]+xvel0[j+1][i])*xarea[j][i];
					dv2=(xvel0[j][i+1]+xvel0[j+1][i+1])*xarea[j][i+1];
					div=dv2-dv1;
					dtut=dtu_safe*2.0*volume[j][i]/(dv1*dv2*volume[j][i]);

					dv1=(yvel0[j][i]+yvel0[j][i+1])*yarea[j][i];
					dv2=(yvel0[j+1][i]+yvel0[j+1][i+1])*yarea[j+1][i];
					div=div+dv2-dv1;
					dtvt=dtv_safe*2.0*volume[j][i]/(dv1*dv2*volume[j][i]);
					div=div/(2.0*volume[j][i]);
					dtdivt=dtdiv_safe*(-1.0/div);
					dt_min[j][i]=dtct*dtut*dtvt*dtdivt;
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
 	printf ("opt: %6lf\n", (double)4096*4096*36*5/(end_time - start_time)/1e9);
}
