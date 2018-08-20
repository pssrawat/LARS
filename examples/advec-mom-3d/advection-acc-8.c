#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock(void);

void advection_opt(double *t_volume, double *t_density1, double *t_vel1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_vol_flux_z, double *t_mass_flux_x, double *t_mass_flux_y, double *t_mass_flux_z, double *t_pre_vol, double *t_post_vol, double *t_node_flux, double *t_node_mass_pre, double *t_node_mass_post, double *t_mom_flux, double *celldx, int N) {
	double(*volume)[258][258] =(double(*)[258][258]) t_volume;
	double(*density1)[258][258] =(double(*)[258][258]) t_density1;
	double(*vel1)[258][258] =(double(*)[258][258]) t_vel1;
	double(*vol_flux_x)[258][258] =(double(*)[258][258]) t_vol_flux_x;
	double(*vol_flux_y)[258][258] =(double(*)[258][258]) t_vol_flux_y;
	double(*vol_flux_z)[258][258] =(double(*)[258][258]) t_vol_flux_z;
	double(*mass_flux_x)[258][258] =(double(*)[258][258]) t_mass_flux_x;
	double(*mass_flux_y)[258][258] =(double(*)[258][258]) t_mass_flux_y;
	double(*mass_flux_z)[258][258] =(double(*)[258][258]) t_mass_flux_z;
	double(*pre_vol)[258][258] =(double(*)[258][258]) t_pre_vol;
	double(*post_vol)[258][258] =(double(*)[258][258]) t_post_vol;
	double(*node_flux)[258][258] =(double(*)[258][258]) t_node_flux;
	double(*mom_flux)[258][258] =(double(*)[258][258]) t_mom_flux;
	double(*node_mass_pre)[258][258] =(double(*)[258][258]) t_node_mass_pre;
	double(*node_mass_post)[258][258] =(double(*)[258][258]) t_node_mass_post;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for(t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=1) {
				for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
					for(i = 1; i < N-1; i++) {
						post_vol[k][j][i]= volume[k][j][i]+vol_flux_y[k][j+1][i]-vol_flux_y[k][j][i] +vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i];
						pre_vol[k][j][i]=post_vol[k][j][i]+vol_flux_x[k][j][i+1]-vol_flux_x[k][j][i];
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=1) {
				for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
					for(i = 1; i < N-1; i++) {
						node_flux[k][j][i] = 0.125*(mass_flux_x[k][j-1][i]+mass_flux_x[k][j][i] +mass_flux_x[k][j-1][i+1]+mass_flux_x[k][j][i+1] +mass_flux_x[k-1][j-1][i]+mass_flux_x[k-1][j][i] +mass_flux_x[k-1][j-1][i+1]+mass_flux_x[k-1][j][i+1]);
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=1) {
				for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
					for(i = 1; i < N-1; i++) {
						node_mass_post[k][j][i] = 0.125*(density1[k][j-1][i]*post_vol[k][j-1][i] +density1[k][j][i]*post_vol[k][j][i] +density1[k][j-1][i-1]*post_vol[k][j-1][i-1] +density1[k][j][i-1]*post_vol[k][j][i-1] +density1[k-1][j-1][i]*post_vol[k-1][j-1][i] +density1[k-1][j][i]*post_vol[k-1][j][i] +density1[k-1][j-1][i-1]*post_vol[k-1][j-1][i-1] +density1[k-1][j][i-1]*post_vol[k-1][j][i-1]);
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=1) {
				for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
					for(i = 1; i < N-1; i++) {
						node_mass_pre[k][j][i] = node_mass_post[k][j][i]-node_flux[k][j][i-1]+node_flux[k][j][i];
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-2; k+=1) {
				for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
					for(i = 1; i < N-1; i++) {
						double sigma=(node_flux[k][j][i])/(node_mass_pre[k][j][i+1]);
						double width=celldx[i];
						double vdiffuw=vel1[k][j][i+1]-vel1[k][j][i+2];
						double vdiffdw=vel1[k][j][i]-vel1[k][j][i+1];
						double limiter=1.0*(width*((2.0-sigma)*vdiffdw/width+(1.0+sigma)*vdiffuw/celldx[i+1])/6.0*vdiffuw*vdiffdw);
						double advec_vel_s=vel1[k][j][i+1]+(1.0-sigma)*limiter;
						mom_flux[k][j][i]=advec_vel_s*node_flux[k][j][i];
					}
				}
			}

		}
	}

	start_time = rtclock();
	for(t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=8) {
				for(j = 1; j < N-1; j++) {
					for(i = 1; i < N-1; i+=1) {
#pragma begin stencil1 unroll k=8,j=1,i=1 print-intrinsics false acc-size 1
						post_vol[k][j][i]= volume[k][j][i]+vol_flux_y[k][j+1][i]-vol_flux_y[k][j][i] +vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i];
						pre_vol[k][j][i]=post_vol[k][j][i]+vol_flux_x[k][j][i+1]-vol_flux_x[k][j][i];
						node_flux[k][j][i] = 0.125*(mass_flux_x[k][j-1][i]+mass_flux_x[k][j][i] +mass_flux_x[k][j-1][i+1]+mass_flux_x[k][j][i+1] +mass_flux_x[k-1][j-1][i]+mass_flux_x[k-1][j][i] +mass_flux_x[k-1][j-1][i+1]+mass_flux_x[k-1][j][i+1]);
#pragma end stencil1
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=8) {
				for(j = 1; j < N-1; j++) {
					for(i = 1; i < N-1; i+=1) {
#pragma begin stencil2 unroll k=8,j=1,i=1 print-intrinsics false acc-size 1
			node_mass_post[k][j][i] = 0.125*(density1[k][j-1][i]*post_vol[k][j-1][i] +density1[k][j][i]*post_vol[k][j][i] +density1[k][j-1][i-1]*post_vol[k][j-1][i-1] +density1[k][j][i-1]*post_vol[k][j][i-1] +density1[k-1][j-1][i]*post_vol[k-1][j-1][i] +density1[k-1][j][i]*post_vol[k-1][j][i] +density1[k-1][j-1][i-1]*post_vol[k-1][j-1][i-1] +density1[k-1][j][i-1]*post_vol[k-1][j][i-1]);
						node_mass_pre[k][j][i] = node_mass_post[k][j][i]-node_flux[k][j][i-1]+node_flux[k][j][i];
#pragma end stencil2
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-2; k+=8) {
				for(j = 1; j < N-1; j++) {
					for(i = 1; i < N-1; i+=1) {
#pragma begin stencil3 unroll k=8,j=1,i=1 print-intrinsics false acc-size 1
						sigma=(node_flux[k][j][i])/(node_mass_pre[k][j][i+1]);
						width=celldx[i];
						vdiffuw=vel1[k][j][i+1]-vel1[k][j][i+2];
						vdiffdw=vel1[k][j][i]-vel1[k][j][i+1];
						limiter=1.0*(width*((2.0-sigma)*vdiffdw/width+(1.0+sigma)*vdiffuw/celldx[i+1])/6.0*vdiffuw*vdiffdw);
						advec_vel_s=vel1[k][j][i+1]+(1.0-sigma)*limiter;
						mom_flux[k][j][i]=advec_vel_s*node_flux[k][j][i];
#pragma end stencil3
					}
				}
			}

#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=8) {
				for(j = 1; j < N-1; j++) {
					for(i = 1; i < N-1; i+=1) {
#pragma begin stencil4 unroll k=8,j=1,i=1 print-intrinsics false acc-size 1
						vel1[k][j][i]=(vel1[k][j][i]*node_mass_pre[k][j][i]+mom_flux[k][j][i-1]-mom_flux[k][j][i])/node_mass_post[k][j][i];
#pragma end stencil4
					}
				}
			}
		}
	}
	end_time = rtclock();
	printf("opt: %6lf\n",(double)256*256*256*55*5/(end_time - start_time)/1e9);
}
