#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock(void);

void advection_opt(double *t_volume, double *t_density1, double *t_vel1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_mass_flux_x, double *t_mass_flux_y, double *t_pre_vol, double *t_post_vol, double *t_node_flux, double *t_node_mass_pre, double *t_node_mass_post, double *t_mom_flux, double *celldx, int N) {
	double(*volume)[4098] =(double(*)[4098]) t_volume;
	double(*density1)[4098] =(double(*)[4098]) t_density1;
	double(*vel1)[4098] =(double(*)[4098]) t_vel1;
	double(*vol_flux_x)[4098] =(double(*)[4098]) t_vol_flux_x;
	double(*vol_flux_y)[4098] =(double(*)[4098]) t_vol_flux_y;
	double(*mass_flux_x)[4098] =(double(*)[4098]) t_mass_flux_x;
	double(*mass_flux_y)[4098] =(double(*)[4098]) t_mass_flux_y;
	double(*pre_vol)[4098] =(double(*)[4098]) t_pre_vol;
	double(*post_vol)[4098] =(double(*)[4098]) t_post_vol;
	double(*node_flux)[4098] =(double(*)[4098]) t_node_flux;
	double(*mom_flux)[4098] =(double(*)[4098]) t_mom_flux;
	double(*node_mass_pre)[4098] =(double(*)[4098]) t_node_mass_pre;
	double(*node_mass_post)[4098] =(double(*)[4098]) t_node_mass_post;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for(t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
				for(i = 1; i < N-1; i++) {
					post_vol[j][i]= volume[j][i]+vol_flux_y[j+1][i]-vol_flux_y[j][i];
					pre_vol[j][i]=post_vol[j][i]+vol_flux_x[j][i+1]-vol_flux_x[j][i];
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
				for(i = 1; i < N-1; i++) {
					node_flux[j][i]=0.25*(mass_flux_x[j-1][i]+mass_flux_x[j][i]  +mass_flux_x[j-1][i+1]+mass_flux_x[j][i+1]);
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
				for(i = 1; i < N-1; i++) {
					node_mass_post[j][i]=0.25*(density1[j-1][i]*post_vol[j-1][i] +density1[j][i]*post_vol[j][i]  +density1[j-1][i-1]*post_vol[j-1][i-1] +density1[j][i-1]*post_vol[j][i-1]);
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
				for(i = 1; i < N-1; i++) {
					node_mass_pre[j][i]=node_mass_post[j][i]-node_flux[j][i-1]+node_flux[j][i];
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
				for(i = 1; i < N-1; i++) {
					double sigma=(node_flux[j][i])/(node_mass_pre[j][i+1]);
					double width=celldx[i];
					double vdiffuw=vel1[j][i+1]-vel1[j][i+2];
					double vdiffdw=vel1[j][i]-vel1[j][i+1];
					double limiter=1.0*(width*((2.0-sigma)*vdiffdw/width+(1.0+sigma)*vdiffuw/celldx[i+1])/6.0*vdiffuw*vdiffdw);
					double advec_vel_s=vel1[j][i+1]+(1.0-sigma)*limiter;
					mom_flux[j][i]=advec_vel_s*node_flux[j][i];
				}
			}
		}
	}

	start_time = rtclock();
	for(t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for(j = 1; j < N-1; j+=4) {
				for(i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll j=4,i=1 print-intrinsics true acc-size 1
					post_vol[j][i]= volume[j][i]+vol_flux_y[j+1][i]-vol_flux_y[j][i];
					pre_vol[j][i]=post_vol[j][i]+vol_flux_x[j][i+1]-vol_flux_x[j][i];
					node_flux[j][i]=0.25*(mass_flux_x[j-1][i]+mass_flux_x[j][i]  +mass_flux_x[j-1][i+1]+mass_flux_x[j][i+1]);
#pragma end stencil1
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j+=4) {
				for(i = 1; i < N-1; i+=8) {
#pragma begin stencil2 unroll j=4,i=1 print-intrinsics true acc-size 1
					node_mass_post[j][i]=0.25*(density1[j-1][i]*post_vol[j-1][i] +density1[j][i]*post_vol[j][i]  +density1[j-1][i-1]*post_vol[j-1][i-1] +density1[j][i-1]*post_vol[j][i-1]);
					node_mass_pre[j][i]=node_mass_post[j][i]-node_flux[j][i-1]+node_flux[j][i];
#pragma end stencil2
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j+=4) {
				for(i = 1; i < N-1; i+=8) {
#pragma begin stencil3 unroll j=4,i=1 print-intrinsics true acc-size 1
					sigma=(node_flux[j][i])/(node_mass_pre[j][i+1]);
					width=celldx[i];
					vdiffuw=vel1[j][i+1]-vel1[j][i+2];
					vdiffdw=vel1[j][i]-vel1[j][i+1];
					limiter=1.0*(width*((2.0-sigma)*vdiffdw/width+(1.0+sigma)*vdiffuw/celldx[i+1])/6.0*vdiffuw*vdiffdw);
					advec_vel_s=vel1[j][i+1]+(1.0-sigma)*limiter;
					mom_flux[j][i]=advec_vel_s*node_flux[j][i];
#pragma end stencil3
				}
			}

#pragma omp for private(i)
			for(j = 1; j < N-1; j+=4) {
				for(i = 1; i < N-1; i+=8) {
#pragma begin stencil4 unroll j=4,i=1 print-intrinsics true acc-size 1
					vel1[j][i]=(vel1[j][i]*node_mass_pre[j][i]+mom_flux[j][i-1]-mom_flux[j][i])/node_mass_post[j][i];
#pragma end stencil4
				}
			}
		}
	}
	end_time = rtclock();
	printf("opt: %6lf\n", (double)4096*4096*41*5/(end_time - start_time)/1e9);
}
