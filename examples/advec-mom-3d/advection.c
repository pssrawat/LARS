#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void advection_opt (double *t_volume, double *t_density1, double *t_vel1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_vol_flux_z, double *t_mass_flux_x, double *t_mass_flux_y, double *t_mass_flux_z, double *t_pre_vol, double *t_post_vol, double *t_node_flux, double *t_node_mass_pre, double *t_node_mass_post, double *t_mom_flux, double *celldx, int N); 

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
	double (*volume)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*vol_flux_x)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*vol_flux_y)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*vol_flux_z)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*density1)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*node_mass_post_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*node_mass_post)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_x)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*mass_flux_y)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*mass_flux_z)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*pre_vol_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_vol_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*pre_vol)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_vol)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*node_flux_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*node_flux)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*node_mass_pre_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*node_mass_pre)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mom_flux_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mom_flux)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*vel1_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*vel1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double * celldx = (double *) getRandom1DArray (258);


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

#pragma omp for private(j,i)
			for(k = 1; k < N-1; k+=1) {
				for(j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
					for(i = 1; i < N-1; i++) {
						vel1[k][j][i]=(vel1[k][j][i]*node_mass_pre[k][j][i]+mom_flux[k][j][i-1]-mom_flux[k][j][i])/node_mass_post[k][j][i];
					}
				}
			}
		}
	}
	end_time = rtclock();
	printf("orig: %6lf\n",(double)256*256*256*55*5/(end_time - start_time)/1e9);

	advection_opt ((double*)volume, (double*)density1, (double*)vel1_opt, (double*)vol_flux_x, (double*)vol_flux_y, (double*)vol_flux_z, (double*)mass_flux_x, (double*)mass_flux_y, (double*)mass_flux_z, (double*)pre_vol_opt, (double*)post_vol_opt, (double*)node_flux_opt, (double*)node_mass_pre_opt, (double*)node_mass_post_opt, (double*) mom_flux_opt, celldx, N);

	double error = checkError3D (N, N, 0, (double*)vel1, (double*)vel1_opt, 2, N-2, 2, N-2, 2, N-2);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
