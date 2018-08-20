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
	int N = 4098;
	double (*dt_min_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*dt_min)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*soundspeed)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*volume)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*density0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*viscosity_a)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*zvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*zarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double *celldx = (double *) getRandom1DArray (4098);
	double *celldy = (double *) getRandom1DArray (4098);
	double *celldz = (double *) getRandom1DArray (4098);

	double dtc_safe =   0.0000311;
	double dtu_safe =   0.0000411;
	double dtv_safe =   0.0000511;
	double dtw_safe =   0.0000611;
	double dtdiv_safe = 0.0000711;

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
					dt_min[j][i]=dtct*dtut*dtvt*dtdivt;
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)4096*4096*36*5/(end_time - start_time)/1e9);

	calc_dt_opt ((double*)dt_min_opt, (double*)volume, (double*)soundspeed, (double*)density0, (double*)viscosity_a, (double*)xvel0, (double*)yvel0, (double*)zvel0, (double*)xarea, (double*)yarea, (double*)zarea, (double*)celldx, (double*)celldy, (double*)celldz, dtc_safe, dtu_safe, dtv_safe, dtw_safe, dtdiv_safe, N);

	double error = checkError2D (N, 0, (double*)dt_min_opt, (double*)dt_min, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n", error);
}
