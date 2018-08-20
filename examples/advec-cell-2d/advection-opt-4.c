#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void advection_opt (double *t_volume, double *t_density1, double *t_energy1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_mass_flux_x, double *t_mass_flux_y, double *t_pre_vol, double *t_post_vol, double *t_pre_mass, double *t_post_mass, double *t_advec_vol, double *t_post_ener, double *t_ener_flux, double *vertexdx, double *vertexdy, int N) {
	double (*volume)[4098] = (double (*)[4098]) t_volume;
	double (*density1)[4098] = (double (*)[4098]) t_density1;
	double (*energy1)[4098] = (double (*)[4098]) t_energy1;
	double (*vol_flux_x)[4098] = (double (*)[4098]) t_vol_flux_x;
	double (*vol_flux_y)[4098] = (double (*)[4098]) t_vol_flux_y;
	double (*mass_flux_x)[4098] = (double (*)[4098]) t_mass_flux_x;
	double (*mass_flux_y)[4098] = (double (*)[4098]) t_mass_flux_y;
	double (*pre_vol)[4098] = (double (*)[4098]) t_pre_vol;
	double (*post_vol)[4098] = (double (*)[4098]) t_post_vol;
	double (*pre_mass)[4098] = (double (*)[4098]) t_pre_mass;
	double (*post_mass)[4098] = (double (*)[4098]) t_post_mass;
	double (*advec_vol)[4098] = (double (*)[4098]) t_advec_vol;
	double (*post_ener)[4098] = (double (*)[4098]) t_post_ener;
	double (*ener_flux)[4098] = (double (*)[4098]) t_ener_flux;

	int t, i, j, k;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 2; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					pre_vol[j][i]=volume[j][i]+(vol_flux_y[j+1][i]-vol_flux_y[j][i]+vol_flux_x[j][i+1]-vol_flux_x[j][i]);
					post_vol[j][i]=pre_vol[j][i]-(vol_flux_y[j+1][i]-vol_flux_y[j][i]);
				}
			}

#pragma omp for private(i)
			for (j = 2; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					double sigmat=(vol_flux_y[j][i])/pre_vol[j-1][i];
					double sigma3=(1.0+sigmat)*(vertexdy[j]/vertexdy[j-1]);
					double sigma4=2.0-sigmat;

					double diffuw=density1[j-1][i]-density1[j-2][i];
					double diffdw=density1[j][i]-density1[j-1][i];
					double limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j][i]=vol_flux_y[j][i]*(density1[j-1][i]+limiter);

					double sigmam=(mass_flux_y[j][i])/(density1[j-1][i]*pre_vol[j-1][i]);
					diffuw=energy1[j-1][i]-energy1[j-2][i];
					diffdw=energy1[j][i]-energy1[j-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j][i]=mass_flux_y[j][i]*(energy1[j-1][i]+limiter);
				}
			}

#pragma omp for private(i)
			for (j = 2; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					double pre_mass_s=density1[j][i]*pre_vol[j][i];
					double post_mass_s=pre_mass_s+mass_flux_y[j][i]-mass_flux_y[j+1][i];
					double post_ener_s=(energy1[j][i]*pre_mass_s+ener_flux[j][i]-ener_flux[j+1][i])/post_mass_s;
					double advec_vol_s=pre_vol[j][i]+vol_flux_y[j][i]-vol_flux_y[j+1][i];
					density1[j][i]=0.1*post_mass_s/advec_vol_s;
					energy1[j][i]=0.1*post_ener_s;
				}
			}
		}
	}

	double start_time, end_time;
	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 2; j < N-1; j+=4) {
				for (i = 1; i < N-1; i+=4) {
#pragma begin stencil1 unroll j=4,i=1 print-intrinsics true acc-size 1
					pre_vol[j][i]=volume[j][i]+(vol_flux_y[j+1][i]-vol_flux_y[j][i]+vol_flux_x[j][i+1]-vol_flux_x[j][i]);
					post_vol[j][i]=pre_vol[j][i]-(vol_flux_y[j+1][i]-vol_flux_y[j][i]);
#pragma end stencil1
				}
			}

#pragma omp for private(i)
			for (j = 2; j < N-1; j+=4) {
				for (i = 1; i < N-1; i+=4) {
#pragma begin stencil2 unroll j=4,i=1 print-intrinsics true acc-size 1
					sigmat=(vol_flux_y[j][i])/pre_vol[j-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j]/vertexdy[j-1]);
					sigma4=2.0-sigmat;

					diffuw=density1[j-1][i]-density1[j-2][i];
					diffdw=density1[j][i]-density1[j-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j][i]=vol_flux_y[j][i]*(density1[j-1][i]+limiter);

					sigmam=(mass_flux_y[j][i])/(density1[j-1][i]*pre_vol[j-1][i]);
					diffuw=energy1[j-1][i]-energy1[j-2][i];
					diffdw=energy1[j][i]-energy1[j-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j][i]=mass_flux_y[j][i]*(energy1[j-1][i]+limiter);

#pragma end stencil2
				}
			}

#pragma omp for private(i)
			for (j = 2; j < N-1; j+=4) {
				for (i = 1; i < N-1; i+=4) {
#pragma begin stencil3 unroll j=4,i=1 print-intrinsics true acc-size 1
					pre_mass_s=density1[j][i]*pre_vol[j][i];
					post_mass_s=pre_mass_s+mass_flux_y[j][i]-mass_flux_y[j+1][i];
					post_ener_s=(energy1[j][i]*pre_mass_s+ener_flux[j][i]-ener_flux[j+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j][i]+vol_flux_y[j][i]-vol_flux_y[j+1][i];
					density1[j][i]=0.1*post_mass_s/advec_vol_s;
					energy1[j][i]=0.1*post_ener_s;
#pragma end stencil3
				}
			}
		}
	}
        end_time = rtclock ();
	printf ("opt: %6lf\n", (double)4096*4096*47*5/(end_time - start_time)/1e9);
}
