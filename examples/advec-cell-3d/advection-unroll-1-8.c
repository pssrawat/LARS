#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void advection_opt (double *t_volume, double *t_density1, double *t_energy1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_vol_flux_z, double *t_mass_flux_x, double *t_mass_flux_y, double *t_mass_flux_z, double *t_pre_vol, double *t_post_vol, double *t_pre_mass, double *t_post_mass, double *t_advec_vol, double *t_post_ener, double *t_ener_flux, double *vertexdx, double *vertexdy, double *vertexdz, int N) {
	double (*volume)[258][258] = (double (*)[258][258]) t_volume;
	double (*density1)[258][258] = (double (*)[258][258]) t_density1;
	double (*energy1)[258][258] = (double (*)[258][258]) t_energy1;
	double (*vol_flux_x)[258][258] = (double (*)[258][258]) t_vol_flux_x;
	double (*vol_flux_y)[258][258] = (double (*)[258][258]) t_vol_flux_y;
	double (*vol_flux_z)[258][258] = (double (*)[258][258]) t_vol_flux_z;
	double (*mass_flux_x)[258][258] = (double (*)[258][258]) t_mass_flux_x;
	double (*mass_flux_y)[258][258] = (double (*)[258][258]) t_mass_flux_y;
	double (*mass_flux_z)[258][258] = (double (*)[258][258]) t_mass_flux_z;
	double (*pre_vol)[258][258] = (double (*)[258][258]) t_pre_vol;
	double (*post_vol)[258][258] = (double (*)[258][258]) t_post_vol;
	double (*pre_mass)[258][258] = (double (*)[258][258]) t_pre_mass;
	double (*post_mass)[258][258] = (double (*)[258][258]) t_post_mass;
	double (*advec_vol)[258][258] = (double (*)[258][258]) t_advec_vol;
	double (*post_ener)[258][258] = (double (*)[258][258]) t_post_ener;
	double (*ener_flux)[258][258] = (double (*)[258][258]) t_ener_flux;

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
				for (j = 1; j < N-1; j+=8) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						pre_vol[k][j][i]=  volume[k][j][i]+(vol_flux_x[k][j][i+1]-vol_flux_x[k][j][i] +vol_flux_y[k][j+1][i]-vol_flux_y[k][j][i] +vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i]);
						post_vol[k][j][i]=pre_vol[k][j][i]-(vol_flux_z[k+1][j][i]-vol_flux_z[k][j][i]);

						pre_vol[k][j+1][i]=  volume[k][j+1][i]+(vol_flux_x[k][j+1][i+1]-vol_flux_x[k][j+1][i] +vol_flux_y[k][j+1+1][i]-vol_flux_y[k][j+1][i] +vol_flux_z[k+1][j+1][i]-vol_flux_z[k][j+1][i]);
						post_vol[k][j+1][i]=pre_vol[k][j+1][i]-(vol_flux_z[k+1][j+1][i]-vol_flux_z[k][j+1][i]);

						pre_vol[k][j+2][i]=  volume[k][j+2][i]+(vol_flux_x[k][j+2][i+1]-vol_flux_x[k][j+2][i] +vol_flux_y[k][j+2+1][i]-vol_flux_y[k][j+2][i] +vol_flux_z[k+1][j+2][i]-vol_flux_z[k][j+2][i]);
						post_vol[k][j+2][i]=pre_vol[k][j+2][i]-(vol_flux_z[k+1][j+2][i]-vol_flux_z[k][j+2][i]);

						pre_vol[k][j+3][i]=  volume[k][j+3][i]+(vol_flux_x[k][j+3][i+1]-vol_flux_x[k][j+3][i] +vol_flux_y[k][j+3+1][i]-vol_flux_y[k][j+3][i] +vol_flux_z[k+1][j+3][i]-vol_flux_z[k][j+3][i]);
						post_vol[k][j+3][i]=pre_vol[k][j+3][i]-(vol_flux_z[k+1][j+3][i]-vol_flux_z[k][j+3][i]);

						pre_vol[k][j+4][i]=  volume[k][j+4][i]+(vol_flux_x[k][j+4][i+1]-vol_flux_x[k][j+4][i] +vol_flux_y[k][j+4+1][i]-vol_flux_y[k][j+4][i] +vol_flux_z[k+1][j+4][i]-vol_flux_z[k][j+4][i]);
						post_vol[k][j+4][i]=pre_vol[k][j+4][i]-(vol_flux_z[k+1][j+4][i]-vol_flux_z[k][j+4][i]);

						pre_vol[k][j+5][i]=  volume[k][j+5][i]+(vol_flux_x[k][j+5][i+1]-vol_flux_x[k][j+5][i] +vol_flux_y[k][j+5+1][i]-vol_flux_y[k][j+5][i] +vol_flux_z[k+1][j+5][i]-vol_flux_z[k][j+5][i]);
						post_vol[k][j+5][i]=pre_vol[k][j+5][i]-(vol_flux_z[k+1][j+5][i]-vol_flux_z[k][j+5][i]);

						pre_vol[k][j+6][i]=  volume[k][j+6][i]+(vol_flux_x[k][j+6][i+1]-vol_flux_x[k][j+6][i] +vol_flux_y[k][j+6+1][i]-vol_flux_y[k][j+6][i] +vol_flux_z[k+1][j+6][i]-vol_flux_z[k][j+6][i]);
						post_vol[k][j+6][i]=pre_vol[k][j+6][i]-(vol_flux_z[k+1][j+6][i]-vol_flux_z[k][j+6][i]);

						pre_vol[k][j+7][i]=  volume[k][j+7][i]+(vol_flux_x[k][j+7][i+1]-vol_flux_x[k][j+7][i] +vol_flux_y[k][j+7+1][i]-vol_flux_y[k][j+7][i] +vol_flux_z[k+1][j+7][i]-vol_flux_z[k][j+7][i]);
						post_vol[k][j+7][i]=pre_vol[k][j+7][i]-(vol_flux_z[k+1][j+7][i]-vol_flux_z[k][j+7][i]);
					}
				}
			}

#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=2) {
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

						sigmat=(vol_flux_z[k+1][j][i])/pre_vol[k+1-1][j][i];
						sigma3=(1.0+sigmat)*(vertexdz[k+1]/vertexdz[k]);
						sigma4=2.0-sigmat;
						diffuw=density1[k+1-1][j][i]-density1[k+1-2][j][i];
						diffdw=density1[k+1][j][i]-density1[k+1-1][j][i];
						limiter=(1.0-sigmat)*diffdw*diffuw*0.1667*(sigma3*diffuw+sigma4*diffdw);
						mass_flux_z[k+1][j][i]=vol_flux_z[k+1][j][i]*(density1[k+1-1][j][i]+limiter);
						sigmam=(mass_flux_z[k+1][j][i])/(density1[k+1-1][j][i]*pre_vol[k+1-1][j][i]);
						diffuw=energy1[k+1-1][j][i]-energy1[k+1-2][j][i];
						diffdw=energy1[k+1][j][i]-energy1[k+1-1][j][i];
						limiter=(1.0-sigmam)*diffdw*diffuw*0.1667*(sigma3*diffuw+sigma4*diffdw);
						ener_flux[k+1][j][i]=mass_flux_z[k+1][j][i]*(energy1[k+1-1][j][i]+limiter);
					}
				}
			}

#pragma omp for private(j,i)
			for (k = 2; k < N-1; k+=1) {
				for (j = 1; j < N-1; j+=8) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 2; i < N-1; i++) {
						pre_mass[k][j][i]=density1[k][j][i]*pre_vol[k][j][i];
						post_mass[k][j][i]=pre_mass[k][j][i]+mass_flux_z[k][j][i]-mass_flux_z[k+1][j][i];
						post_ener[k][j][i]=(energy1[k][j][i]*pre_mass[k][j][i]+ener_flux[k][j][i]-ener_flux[k+1][j][i])/post_mass[k][j][i];
						advec_vol[k][j][i]=pre_vol[k][j][i]+vol_flux_z[k][j][i]-vol_flux_z[k+1][j][i];
						density1[k][j][i]=post_mass[k][j][i]/advec_vol[k][j][i];
						energy1[k][j][i]=post_ener[k][j][i];

                                                pre_mass[k][j+1][i]=density1[k][j+1][i]*pre_vol[k][j+1][i];
                                                post_mass[k][j+1][i]=pre_mass[k][j+1][i]+mass_flux_z[k][j+1][i]-mass_flux_z[k+1][j+1][i];
                                                post_ener[k][j+1][i]=(energy1[k][j+1][i]*pre_mass[k][j+1][i]+ener_flux[k][j+1][i]-ener_flux[k+1][j+1][i])/post_mass[k][j+1][i];
                                                advec_vol[k][j+1][i]=pre_vol[k][j+1][i]+vol_flux_z[k][j+1][i]-vol_flux_z[k+1][j+1][i];
                                                density1[k][j+1][i]=post_mass[k][j+1][i]/advec_vol[k][j+1][i];
                                                energy1[k][j+1][i]=post_ener[k][j+1][i];

                                                pre_mass[k][j+2][i]=density1[k][j+2][i]*pre_vol[k][j+2][i];
                                                post_mass[k][j+2][i]=pre_mass[k][j+2][i]+mass_flux_z[k][j+2][i]-mass_flux_z[k+1][j+2][i];
                                                post_ener[k][j+2][i]=(energy1[k][j+2][i]*pre_mass[k][j+2][i]+ener_flux[k][j+2][i]-ener_flux[k+1][j+2][i])/post_mass[k][j+2][i];
                                                advec_vol[k][j+2][i]=pre_vol[k][j+2][i]+vol_flux_z[k][j+2][i]-vol_flux_z[k+1][j+2][i];
                                                density1[k][j+2][i]=post_mass[k][j+2][i]/advec_vol[k][j+2][i];
                                                energy1[k][j+2][i]=post_ener[k][j+2][i];

                                                pre_mass[k][j+3][i]=density1[k][j+3][i]*pre_vol[k][j+3][i];
                                                post_mass[k][j+3][i]=pre_mass[k][j+3][i]+mass_flux_z[k][j+3][i]-mass_flux_z[k+1][j+3][i];
                                                post_ener[k][j+3][i]=(energy1[k][j+3][i]*pre_mass[k][j+3][i]+ener_flux[k][j+3][i]-ener_flux[k+1][j+3][i])/post_mass[k][j+3][i];
                                                advec_vol[k][j+3][i]=pre_vol[k][j+3][i]+vol_flux_z[k][j+3][i]-vol_flux_z[k+1][j+3][i];
                                                density1[k][j+3][i]=post_mass[k][j+3][i]/advec_vol[k][j+3][i];
                                                energy1[k][j+3][i]=post_ener[k][j+3][i];

                                                pre_mass[k][j+4][i]=density1[k][j+4][i]*pre_vol[k][j+4][i];
                                                post_mass[k][j+4][i]=pre_mass[k][j+4][i]+mass_flux_z[k][j+4][i]-mass_flux_z[k+1][j+4][i];
                                                post_ener[k][j+4][i]=(energy1[k][j+4][i]*pre_mass[k][j+4][i]+ener_flux[k][j+4][i]-ener_flux[k+1][j+4][i])/post_mass[k][j+4][i];
                                                advec_vol[k][j+4][i]=pre_vol[k][j+4][i]+vol_flux_z[k][j+4][i]-vol_flux_z[k+1][j+4][i];
                                                density1[k][j+4][i]=post_mass[k][j+4][i]/advec_vol[k][j+4][i];
                                                energy1[k][j+4][i]=post_ener[k][j+4][i];

                                                pre_mass[k][j+5][i]=density1[k][j+5][i]*pre_vol[k][j+5][i];
                                                post_mass[k][j+5][i]=pre_mass[k][j+5][i]+mass_flux_z[k][j+5][i]-mass_flux_z[k+1][j+5][i];
                                                post_ener[k][j+5][i]=(energy1[k][j+5][i]*pre_mass[k][j+5][i]+ener_flux[k][j+5][i]-ener_flux[k+1][j+5][i])/post_mass[k][j+5][i];
                                                advec_vol[k][j+5][i]=pre_vol[k][j+5][i]+vol_flux_z[k][j+5][i]-vol_flux_z[k+1][j+5][i];
                                                density1[k][j+5][i]=post_mass[k][j+5][i]/advec_vol[k][j+5][i];
                                                energy1[k][j+5][i]=post_ener[k][j+5][i];

                                                pre_mass[k][j+6][i]=density1[k][j+6][i]*pre_vol[k][j+6][i];
                                                post_mass[k][j+6][i]=pre_mass[k][j+6][i]+mass_flux_z[k][j+6][i]-mass_flux_z[k+1][j+6][i];
                                                post_ener[k][j+6][i]=(energy1[k][j+6][i]*pre_mass[k][j+6][i]+ener_flux[k][j+6][i]-ener_flux[k+1][j+6][i])/post_mass[k][j+6][i];
                                                advec_vol[k][j+6][i]=pre_vol[k][j+6][i]+vol_flux_z[k][j+6][i]-vol_flux_z[k+1][j+6][i];
                                                density1[k][j+6][i]=post_mass[k][j+6][i]/advec_vol[k][j+6][i];
                                                energy1[k][j+6][i]=post_ener[k][j+6][i];

                                                pre_mass[k][j+7][i]=density1[k][j+7][i]*pre_vol[k][j+7][i];
                                                post_mass[k][j+7][i]=pre_mass[k][j+7][i]+mass_flux_z[k][j+7][i]-mass_flux_z[k+1][j+7][i];
                                                post_ener[k][j+7][i]=(energy1[k][j+7][i]*pre_mass[k][j+7][i]+ener_flux[k][j+7][i]-ener_flux[k+1][j+7][i])/post_mass[k][j+7][i];
                                                advec_vol[k][j+7][i]=pre_vol[k][j+7][i]+vol_flux_z[k][j+7][i]-vol_flux_z[k+1][j+7][i];
                                                density1[k][j+7][i]=post_mass[k][j+7][i]/advec_vol[k][j+7][i];
                                                energy1[k][j+7][i]=post_ener[k][j+7][i];
					}
				}
			}
		}	
	}
	end_time = rtclock ();
	printf ("unroll: %6lf\n", (double)256*256*256*49*5/(end_time - start_time)/1e9);
}
