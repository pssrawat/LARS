#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void advection_opt (double *t_volume, double *t_density1, double *t_energy1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_vol_flux_z, double *t_mass_flux_x, double *t_mass_flux_y, double *t_mass_flux_z, double *t_pre_vol, double *t_post_vol, double *t_pre_mass, double *t_post_mass, double *t_advec_vol, double *t_post_ener, double *t_ener_flux, double *vertexdx, double *vertexdy, double *vertexdz, int N); 

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

	double (*density1_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*density1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*energy1_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*energy1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_x_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_y_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_z_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_x)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_y)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*mass_flux_z)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*pre_vol_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_vol_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*pre_vol)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_vol)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*pre_mass_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_mass_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*pre_mass)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_mass)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*advec_vol_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_ener_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*ener_flux_out)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*advec_vol)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*post_ener)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*ener_flux)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double *vertexdx = (double *) getRandom1DArray (258);
	double *vertexdy = (double *) getRandom1DArray (258);
	double *vertexdz = (double *) getRandom1DArray (258);

	int t, i, j, k;
	double start_time, end_time;

	// Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						pre_vol[k][j][i]=  volume[k][j][i]+(vol_flux_x[k][j][i+1]-vol_flux_x[k][j][i] +vol_flux_y[k][j+1][i]-vol_flux_y[k][j][i] +vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i]);
						post_vol[k][j][i]=pre_vol[k][j][i]-(vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i]);
					}
				}
			}

#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						double sigmat=(vol_flux_z[k][j][i])/pre_vol[k-1][j][i];
						double sigma3=(1.0+sigmat)*(vertexdz[k]/vertexdz[k-1]);
						double sigma4=2.0-sigmat;
						double diffuw=density1[k-1][j][i]-density1[k-2][j][i];
						double diffdw=density1[k][j][i]-density1[k-1][j][i];
						double limiter=(1.0-sigmat)*diffdw*diffuw*0.1667*(sigma3*diffuw+sigma4*diffdw);
						mass_flux_z[k][j][i]=vol_flux_z[k][j][i]*(density1[k-1][j][i]+limiter);
						double sigmam=(mass_flux_z[k][j][i])/(density1[k-1][j][i]*pre_vol[k-1][j][i]);
						diffuw=energy1[k-1][j][i]-energy1[k-2][j][i];
						diffdw=energy1[k][j][i]-energy1[k-1][j][i];
						limiter=(1.0-sigmam)*diffdw*diffuw*0.1667*(sigma3*diffuw+sigma4*diffdw);
						ener_flux[k][j][i]=mass_flux_z[k][j][i]*(energy1[k-1][j][i]+limiter);

					}
				}
			}

#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 2; i < N-1; i++) {
						pre_mass[k][j][i]=0.1 * density1[k][j][i]*pre_vol[k][j][i];
						post_mass[k][j][i]=pre_mass[k][j][i]+mass_flux_z[k][j][i]-mass_flux_z[k+1][j][i];
						post_ener[k][j][i]=(energy1[k][j][i]*pre_mass[k][j][i]+ener_flux[k][j][i]-ener_flux[k+1][j][i])/post_mass[k][j][i];
						advec_vol[k][j][i]=pre_vol[k][j][i]+vol_flux_z[k][j][i]-vol_flux_z[k+1][j][i];
						density1[k][j][i]=post_mass[k][j][i]/advec_vol[k][j][i];
						energy1[k][j][i]=post_ener[k][j][i];
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
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						pre_vol[k][j][i]=  volume[k][j][i]+(vol_flux_x[k][j][i+1]-vol_flux_x[k][j][i] +vol_flux_y[k][j+1][i]-vol_flux_y[k][j][i] +vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i]);
						post_vol[k][j][i]=pre_vol[k][j][i]-(vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i]);
					}
				}
			}

#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						double sigmat=(vol_flux_z[k][j][i])/pre_vol[k-1][j][i];
						double sigma3=(1.0+sigmat)*(vertexdz[k]/vertexdz[k-1]);
						double sigma4=2.0-sigmat;
						double diffuw=density1[k-1][j][i]-density1[k-2][j][i];
						double diffdw=density1[k][j][i]-density1[k-1][j][i];
						double limiter=(1.0-sigmat)*diffdw*diffuw*0.1667*(sigma3*diffuw+sigma4*diffdw);
						mass_flux_z[k][j][i]=vol_flux_z[k][j][i]*(density1[k-1][j][i]+limiter);
						double sigmam=(mass_flux_z[k][j][i])/(density1[k-1][j][i]*pre_vol[k-1][j][i]);
						diffuw=energy1[k-1][j][i]-energy1[k-2][j][i];
						diffdw=energy1[k][j][i]-energy1[k-1][j][i];
						limiter=(1.0-sigmam)*diffdw*diffuw*0.1667*(sigma3*diffuw+sigma4*diffdw);
						ener_flux[k][j][i]=mass_flux_z[k][j][i]*(energy1[k-1][j][i]+limiter);

					}
				}
			}

#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 2; i < N-1; i++) {
						pre_mass[k][j][i]=density1[k][j][i]*pre_vol[k][j][i];
						post_mass[k][j][i]=pre_mass[k][j][i]+mass_flux_z[k][j][i]-mass_flux_z[k+1][j][i];
						post_ener[k][j][i]=(energy1[k][j][i]*pre_mass[k][j][i]+ener_flux[k][j][i]-ener_flux[k+1][j][i])/post_mass[k][j][i];
						advec_vol[k][j][i]=pre_vol[k][j][i]+vol_flux_z[k][j][i]-vol_flux_z[k+1][j][i];
						density1[k][j][i]=post_mass[k][j][i]/advec_vol[k][j][i];
						energy1[k][j][i]=post_ener[k][j][i];
					}
				}
			}
		}	
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*49*5/(end_time - start_time)/1e9);

	advection_opt ((double*)volume, (double*)density1_out, (double*)energy1_out, (double*)vol_flux_x, (double*)vol_flux_y, (double*)vol_flux_z, (double*)mass_flux_x_out, (double*)mass_flux_y_out, (double*)mass_flux_z_out, (double*)pre_vol_out, (double*) post_vol_out, (double*)pre_mass_out, (double*)post_mass_out, (double*)advec_vol_out, (double*)post_ener_out, (double*)ener_flux_out, vertexdx, vertexdy, vertexdz, N);

	double error = checkError3D (N, N, 0, (double*)pre_mass, (double*)pre_mass_out, 2, N-2, 2, N-2, 2, N-2);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
