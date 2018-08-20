#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void advection_opt (double *t_volume, double *t_density1, double *t_energy1, double *t_vol_flux_x, double *t_vol_flux_y, double *t_mass_flux_x, double *t_mass_flux_y, double *t_pre_vol, double *t_post_vol, double *t_pre_mass, double *t_post_mass, double *t_advec_vol, double *t_post_ener, double *t_ener_flux, double *vertexdx, double *vertexdy, int N); 

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
	double (*volume)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*vol_flux_x)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*vol_flux_y)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);

	double (*density1_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*density1)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*energy1_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*energy1)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*mass_flux_x_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*mass_flux_y_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*mass_flux_x)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*mass_flux_y)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*pre_vol_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*post_vol_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*pre_vol)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*post_vol)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*pre_mass_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*post_mass_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*pre_mass)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*post_mass)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double (*advec_vol_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*post_ener_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*ener_flux_opt)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*advec_vol)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*post_ener)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*ener_flux)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

	double * vertexdx = (double *) getRandom1DArray (4098);
	double * vertexdy = (double *) getRandom1DArray (4098);


	int t, i, j, k;
	double start_time, end_time;

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

	start_time = rtclock ();
	for (t=0; t<5; t++) {
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
					density1[j][i]=post_mass_s/advec_vol_s;
					energy1[j][i]=post_ener_s;
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)4096*4096*47*5/(end_time - start_time)/1e9);

	advection_opt ((double*)volume, (double*)density1_opt, (double*)energy1_opt, (double*)vol_flux_x, (double*)vol_flux_y, (double*)mass_flux_x_opt, (double*)mass_flux_y_opt, (double*)pre_vol_opt, (double*)post_vol_opt, (double*)pre_mass_opt, (double*)post_mass_opt, (double*)advec_vol_opt, (double*)post_ener_opt, (double*)ener_flux_opt, vertexdx, vertexdy, N);

	double error = checkError2D (N, 0, (double*)energy1, (double*)energy1_opt, 2, N-2, 2, N-2);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
	error = checkError2D (N, 0, (double*)density1, (double*)density1_opt, 2, N-2, 2, N-2);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
