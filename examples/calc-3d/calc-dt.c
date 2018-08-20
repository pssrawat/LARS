#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void calc_dt_opt (double *t_dt_min, double *t_volume, double *t_soundspeed, double *t_density0, double *t_viscosity_a, double *t_xvel0, double *t_yvel0, double *t_zvel0, double *t_xarea, double *t_yarea, double *t_zarea, double *celldx, double *celldy, double *celldz, double dtc_safe, double dtu_safe, double dtv_safe, double dtw_safe, double dtdiv_safe, int N); 

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
	double (*dt_min_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*dt_min)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*soundspeed)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*volume)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*density0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*viscosity_a)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double *celldx = (double *) getRandom1DArray (258);
	double *celldy = (double *) getRandom1DArray (258);
	double *celldz = (double *) getRandom1DArray (258);

	double dtc_safe =   0.00000311;
	double dtu_safe =   0.00000411;
	double dtv_safe =   0.00000511;
	double dtw_safe =   0.00000611;
	double dtdiv_safe = 0.00000711;

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
						dt_min[k][j][i]=dtct*dtut*dtvt*dtwt*dtdivt;
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*67*5/(end_time - start_time)/1e9);

	calc_dt_opt ((double*)dt_min_opt, (double*)volume, (double*)soundspeed, (double*)density0, (double*)viscosity_a, (double*)xvel0, (double*)yvel0, (double*)zvel0, (double*)xarea, (double*)yarea, (double*)zarea, (double*)celldx, (double*)celldy, (double*)celldz, dtc_safe, dtu_safe, dtv_safe, dtw_safe, dtdiv_safe, N);

	double error = checkError3D (N, N, 0, (double*)dt_min_opt, (double*)dt_min, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n", error);
}
