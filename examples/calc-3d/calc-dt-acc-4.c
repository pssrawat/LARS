#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern double rtclock (void);

void calc_dt_opt (double *t_dt_min, double *t_volume, double *t_soundspeed, double *t_density0, double *t_viscosity_a, double *t_xvel0, double *t_yvel0, double *t_zvel0, double *t_xarea, double *t_yarea, double *t_zarea, double *celldx, double *celldy, double *celldz, double dtc_safe, double dtu_safe, double dtv_safe, double dtw_safe, double dtdiv_safe, int N) {
	double (*dt_min)[258][258] = (double (*)[258][258]) t_dt_min;
	double (*volume)[258][258] = (double (*)[258][258]) t_volume;
	double (*soundspeed)[258][258] = (double (*)[258][258]) t_soundspeed;
	double (*density0)[258][258] = (double (*)[258][258]) t_density0;
	double (*viscosity_a)[258][258] = (double (*)[258][258]) t_viscosity_a;
	double (*xvel0)[258][258] = (double (*)[258][258]) t_xvel0;
	double (*yvel0)[258][258] = (double (*)[258][258]) t_yvel0;
	double (*zvel0)[258][258] = (double (*)[258][258]) t_zvel0;
	double (*xarea)[258][258] = (double (*)[258][258]) t_xarea;
	double (*yarea)[258][258] = (double (*)[258][258]) t_yarea;
	double (*zarea)[258][258] = (double (*)[258][258]) t_zarea;

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
						double ds=1.0/(celldx[i]*celldy[j]*celldz[k]);
						double cc=soundspeed[k][j][i]*soundspeed[k][j][i];
						cc=cc+2.0*viscosity_a[k][j][i]/density0[k][j][i];
						double dtct=ds*cc;
						dtct=dtc_safe*1.0/sqrt(dtct);
						double du1=(xvel0[k][j][i]+xvel0[k][j+1][i]+xvel0[k+1][j][i]+xvel0[k+1][j+1][i])*xarea[k][j][i];
						double du2=(xvel0[k][j][i+1]+xvel0[k][j+1][i+1]+xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1])*xarea[k][j][i];
						double dtut=dtu_safe*4.0*volume[k][j][i]/(du1*du2*0.00001*volume[k][j][i]);
						double dv1=(yvel0[k][j][i]+yvel0[k][j][i+1]+yvel0[k+1][j][i]+yvel0[k+1][j][i+1])*yarea[k][j][i];
						double dv2=(yvel0[k][j+1][i]+yvel0[k][j+1][i+1]+yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1])*yarea[k][j][i];
						double dtvt=dtv_safe*4.0*volume[k][j][i]/(dv1*dv2*0.00001*volume[k][j][i]);
						double dw1=(zvel0[k][j][i]+zvel0[k][j+1][i]+zvel0[k][j][i+1]+zvel0[k][j+1][i+1])*zarea[k][j][i];
						double dw2=(zvel0[k+1][j][i]+zvel0[k+1][j+1][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i+1])*zarea[k][j][i];
						double dtwt=dtw_safe*4.0*volume[k][j][i]/(dw1*dw2*0.00001*volume[k][j][i]);
						double div=du2-du1+dv2-dv1+dw2-dw1;
						double dtdivt=dtdiv_safe*4.0*volume[k][j][i]/(volume[k][j][i]*0.00001*div);
						dt_min[k][j][i]=0.1*dtct*dtut*dtvt*dtwt*dtdivt;
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
			for (k = 1; k < N-1; k+=4) {
				for (j = 1; j < N-1; j++) {
					for (i = 1; i < N-1; i+=1) {
#pragma begin stencil1 unroll k=1,j=1,i=1 print-intrinsics false acc-size 1
						ds=1.0/(celldx[i]*celldy[j]*celldz[k]);
						cc=soundspeed[k][j][i]*soundspeed[k][j][i];
						cc+=2.0*viscosity_a[k][j][i]/density0[k][j][i];
						dtct=ds*cc;
						dtct=dtc_safe*1.0/sqrt(dtct);
						du1=(xvel0[k][j][i]+xvel0[k][j+1][i]+xvel0[k+1][j][i]+xvel0[k+1][j+1][i])*xarea[k][j][i];
						du2=(xvel0[k][j][i+1]+xvel0[k][j+1][i+1]+xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1])*xarea[k][j][i];
						dtut=dtu_safe*4.0*volume[k][j][i]/(du1*du2*0.00001*volume[k][j][i]);
						dv1=(yvel0[k][j][i]+yvel0[k][j][i+1]+yvel0[k+1][j][i]+yvel0[k+1][j][i+1])*yarea[k][j][i];
						dv2=(yvel0[k][j+1][i]+yvel0[k][j+1][i+1]+yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1])*yarea[k][j][i];
						dtvt=dtv_safe*4.0*volume[k][j][i]/(dv1*dv2*0.00001*volume[k][j][i]);
						dw1=(zvel0[k][j][i]+zvel0[k][j+1][i]+zvel0[k][j][i+1]+zvel0[k][j+1][i+1])*zarea[k][j][i];
						dw2=(zvel0[k+1][j][i]+zvel0[k+1][j+1][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i+1])*zarea[k][j][i];
						dtwt=dtw_safe*4.0*volume[k][j][i]/(dw1*dw2*0.00001*volume[k][j][i]);
						div=du2-du1+dv2-dv1+dw2-dw1;
						dtdivt=dtdiv_safe*4.0*volume[k][j][i]/(volume[k][j][i]*0.00001*div);
						dt_min[k][j][i]=dtct*dtut*dtvt*dtwt*dtdivt;

						ds=1.0/(celldx[i]*celldy[j]*celldz[k+1]);
						cc=soundspeed[k+1][j][i]*soundspeed[k+1][j][i];
						cc+=2.0*viscosity_a[k+1][j][i]/density0[k+1][j][i];
						dtct=ds*cc;
						dtct=dtc_safe*1.0/sqrt(dtct);
						du1=(xvel0[k+1][j][i]+xvel0[k+1][j+1][i]+xvel0[k+1+1][j][i]+xvel0[k+1+1][j+1][i])*xarea[k+1][j][i];
						du2=(xvel0[k+1][j][i+1]+xvel0[k+1][j+1][i+1]+xvel0[k+1+1][j][i+1]+xvel0[k+1+1][j+1][i+1])*xarea[k+1][j][i];
						dtut=dtu_safe*4.0*volume[k+1][j][i]/(du1*du2*0.00001*volume[k+1][j][i]);
						dv1=(yvel0[k+1][j][i]+yvel0[k+1][j][i+1]+yvel0[k+1+1][j][i]+yvel0[k+1+1][j][i+1])*yarea[k+1][j][i];
						dv2=(yvel0[k+1][j+1][i]+yvel0[k+1][j+1][i+1]+yvel0[k+1+1][j+1][i]+yvel0[k+1+1][j+1][i+1])*yarea[k+1][j][i];
						dtvt=dtv_safe*4.0*volume[k+1][j][i]/(dv1*dv2*0.00001*volume[k+1][j][i]);
						dw1=(zvel0[k+1][j][i]+zvel0[k+1][j+1][i]+zvel0[k+1][j][i+1]+zvel0[k+1][j+1][i+1])*zarea[k+1][j][i];
						dw2=(zvel0[k+1+1][j][i]+zvel0[k+1+1][j+1][i]+zvel0[k+1+1][j][i+1]+zvel0[k+1+1][j+1][i+1])*zarea[k+1][j][i];
						dtwt=dtw_safe*4.0*volume[k+1][j][i]/(dw1*dw2*0.00001*volume[k+1][j][i]);
						div=du2-du1+dv2-dv1+dw2-dw1;
						dtdivt=dtdiv_safe*4.0*volume[k+1][j][i]/(volume[k+1][j][i]*0.00001*div);
						dt_min[k+1][j][i]=dtct*dtut*dtvt*dtwt*dtdivt;

						ds=1.0/(celldx[i]*celldy[j]*celldz[k+2]);
						cc=soundspeed[k+2][j][i]*soundspeed[k+2][j][i];
						cc+=2.0*viscosity_a[k+2][j][i]/density0[k+2][j][i];
						dtct=ds*cc;
						dtct=dtc_safe*1.0/sqrt(dtct);
						du1=(xvel0[k+2][j][i]+xvel0[k+2][j+1][i]+xvel0[k+2+1][j][i]+xvel0[k+2+1][j+1][i])*xarea[k+2][j][i];
						du2=(xvel0[k+2][j][i+1]+xvel0[k+2][j+1][i+1]+xvel0[k+2+1][j][i+1]+xvel0[k+2+1][j+1][i+1])*xarea[k+2][j][i];
						dtut=dtu_safe*4.0*volume[k+2][j][i]/(du1*du2*0.00001*volume[k+2][j][i]);
						dv1=(yvel0[k+2][j][i]+yvel0[k+2][j][i+1]+yvel0[k+2+1][j][i]+yvel0[k+2+1][j][i+1])*yarea[k+2][j][i];
						dv2=(yvel0[k+2][j+1][i]+yvel0[k+2][j+1][i+1]+yvel0[k+2+1][j+1][i]+yvel0[k+2+1][j+1][i+1])*yarea[k+2][j][i];
						dtvt=dtv_safe*4.0*volume[k+2][j][i]/(dv1*dv2*0.00001*volume[k+2][j][i]);
						dw1=(zvel0[k+2][j][i]+zvel0[k+2][j+1][i]+zvel0[k+2][j][i+1]+zvel0[k+2][j+1][i+1])*zarea[k+2][j][i];
						dw2=(zvel0[k+2+1][j][i]+zvel0[k+2+1][j+1][i]+zvel0[k+2+1][j][i+1]+zvel0[k+2+1][j+1][i+1])*zarea[k+2][j][i];
						dtwt=dtw_safe*4.0*volume[k+2][j][i]/(dw1*dw2*0.00001*volume[k+2][j][i]);
						div=du2-du1+dv2-dv1+dw2-dw1;
						dtdivt=dtdiv_safe*4.0*volume[k+2][j][i]/(volume[k+2][j][i]*0.00001*div);
						dt_min[k+2][j][i]=dtct*dtut*dtvt*dtwt*dtdivt;

						ds=1.0/(celldx[i]*celldy[j]*celldz[k+3]);
						cc=soundspeed[k+3][j][i]*soundspeed[k+3][j][i];
						cc+=2.0*viscosity_a[k+3][j][i]/density0[k+3][j][i];
						dtct=ds*cc;
						dtct=dtc_safe*1.0/sqrt(dtct);
						du1=(xvel0[k+3][j][i]+xvel0[k+3][j+1][i]+xvel0[k+3+1][j][i]+xvel0[k+3+1][j+1][i])*xarea[k+3][j][i];
						du2=(xvel0[k+3][j][i+1]+xvel0[k+3][j+1][i+1]+xvel0[k+3+1][j][i+1]+xvel0[k+3+1][j+1][i+1])*xarea[k+3][j][i];
						dtut=dtu_safe*4.0*volume[k+3][j][i]/(du1*du2*0.00001*volume[k+3][j][i]);
						dv1=(yvel0[k+3][j][i]+yvel0[k+3][j][i+1]+yvel0[k+3+1][j][i]+yvel0[k+3+1][j][i+1])*yarea[k+3][j][i];
						dv2=(yvel0[k+3][j+1][i]+yvel0[k+3][j+1][i+1]+yvel0[k+3+1][j+1][i]+yvel0[k+3+1][j+1][i+1])*yarea[k+3][j][i];
						dtvt=dtv_safe*4.0*volume[k+3][j][i]/(dv1*dv2*0.00001*volume[k+3][j][i]);
						dw1=(zvel0[k+3][j][i]+zvel0[k+3][j+1][i]+zvel0[k+3][j][i+1]+zvel0[k+3][j+1][i+1])*zarea[k+3][j][i];
						dw2=(zvel0[k+3+1][j][i]+zvel0[k+3+1][j+1][i]+zvel0[k+3+1][j][i+1]+zvel0[k+3+1][j+1][i+1])*zarea[k+3][j][i];
						dtwt=dtw_safe*4.0*volume[k+3][j][i]/(dw1*dw2*0.00001*volume[k+3][j][i]);
						div=du2-du1+dv2-dv1+dw2-dw1;
						dtdivt=dtdiv_safe*4.0*volume[k+3][j][i]/(volume[k+3][j][i]*0.00001*div);
						dt_min[k+3][j][i]=dtct*dtut*dtvt*dtwt*dtdivt;
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)256*256*256*67*5/(end_time - start_time)/1e9);
}
