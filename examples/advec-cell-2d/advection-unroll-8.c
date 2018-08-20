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
			for (j = 2; j < N-1; j+=8) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					pre_vol[j][i]=volume[j][i]+(vol_flux_y[j+1][i]-vol_flux_y[j][i]+vol_flux_x[j][i+1]-vol_flux_x[j][i]);
					post_vol[j][i]=pre_vol[j][i]-(vol_flux_y[j+1][i]-vol_flux_y[j][i]);

					pre_vol[j+1][i]=volume[j+1][i]+(vol_flux_y[j+1+1][i]-vol_flux_y[j+1][i]+vol_flux_x[j+1][i+1]-vol_flux_x[j+1][i]);
					post_vol[j+1][i]=pre_vol[j+1][i]-(vol_flux_y[j+1+1][i]-vol_flux_y[j+1][i]);

					pre_vol[j+2][i]=volume[j+2][i]+(vol_flux_y[j+2+1][i]-vol_flux_y[j+2][i]+vol_flux_x[j+2][i+1]-vol_flux_x[j+2][i]);
					post_vol[j+2][i]=pre_vol[j+2][i]-(vol_flux_y[j+2+1][i]-vol_flux_y[j+2][i]);

					pre_vol[j+3][i]=volume[j+3][i]+(vol_flux_y[j+3+1][i]-vol_flux_y[j+3][i]+vol_flux_x[j+3][i+1]-vol_flux_x[j+3][i]);
					post_vol[j+3][i]=pre_vol[j+3][i]-(vol_flux_y[j+3+1][i]-vol_flux_y[j+3][i]);

					pre_vol[j+4][i]=volume[j+4][i]+(vol_flux_y[j+4+1][i]-vol_flux_y[j+4][i]+vol_flux_x[j+4][i+1]-vol_flux_x[j+4][i]);
					post_vol[j+4][i]=pre_vol[j+4][i]-(vol_flux_y[j+4+1][i]-vol_flux_y[j+4][i]);

					pre_vol[j+4+1][i]=volume[j+4+1][i]+(vol_flux_y[j+4+1+1][i]-vol_flux_y[j+4+1][i]+vol_flux_x[j+4+1][i+1]-vol_flux_x[j+4+1][i]);
					post_vol[j+4+1][i]=pre_vol[j+4+1][i]-(vol_flux_y[j+4+1+1][i]-vol_flux_y[j+4+1][i]);

					pre_vol[j+4+2][i]=volume[j+4+2][i]+(vol_flux_y[j+4+2+1][i]-vol_flux_y[j+4+2][i]+vol_flux_x[j+4+2][i+1]-vol_flux_x[j+4+2][i]);
					post_vol[j+4+2][i]=pre_vol[j+4+2][i]-(vol_flux_y[j+4+2+1][i]-vol_flux_y[j+4+2][i]);

					pre_vol[j+4+3][i]=volume[j+4+3][i]+(vol_flux_y[j+4+3+1][i]-vol_flux_y[j+4+3][i]+vol_flux_x[j+4+3][i+1]-vol_flux_x[j+4+3][i]);
					post_vol[j+4+3][i]=pre_vol[j+4+3][i]-(vol_flux_y[j+4+3+1][i]-vol_flux_y[j+4+3][i]);
				}
			}

#pragma omp for private(i)
			for (j = 2; j < N-1; j+=8) {
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

					sigmat=(vol_flux_y[j+1][i])/pre_vol[j+1-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+1]/vertexdy[j+1-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+1-1][i]-density1[j+1-2][i];
					diffdw=density1[j+1][i]-density1[j+1-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+1][i]=vol_flux_y[j+1][i]*(density1[j+1-1][i]+limiter);
					sigmam=(mass_flux_y[j+1][i])/(density1[j+1-1][i]*pre_vol[j+1-1][i]);
					diffuw=energy1[j+1-1][i]-energy1[j+1-2][i];
					diffdw=energy1[j+1][i]-energy1[j+1-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+1][i]=mass_flux_y[j+1][i]*(energy1[j+1-1][i]+limiter);

					sigmat=(vol_flux_y[j+2][i])/pre_vol[j+2-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+2]/vertexdy[j+2-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+2-1][i]-density1[j+2-2][i];
					diffdw=density1[j+2][i]-density1[j+2-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+2][i]=vol_flux_y[j+2][i]*(density1[j+2-1][i]+limiter);
					sigmam=(mass_flux_y[j+2][i])/(density1[j+2-1][i]*pre_vol[j+2-1][i]);
					diffuw=energy1[j+2-1][i]-energy1[j+2-2][i];
					diffdw=energy1[j+2][i]-energy1[j+2-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+2][i]=mass_flux_y[j+2][i]*(energy1[j+2-1][i]+limiter);

					sigmat=(vol_flux_y[j+3][i])/pre_vol[j+3-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+3]/vertexdy[j+3-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+3-1][i]-density1[j+3-2][i];
					diffdw=density1[j+3][i]-density1[j+3-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+3][i]=vol_flux_y[j+3][i]*(density1[j+3-1][i]+limiter);
					sigmam=(mass_flux_y[j+3][i])/(density1[j+3-1][i]*pre_vol[j+3-1][i]);
					diffuw=energy1[j+3-1][i]-energy1[j+3-2][i];
					diffdw=energy1[j+3][i]-energy1[j+3-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+3][i]=mass_flux_y[j+3][i]*(energy1[j+3-1][i]+limiter);

					sigmat=(vol_flux_y[j+4][i])/pre_vol[j+4-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+4]/vertexdy[j+4-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+4-1][i]-density1[j+4-2][i];
					diffdw=density1[j+4][i]-density1[j+4-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+4][i]=vol_flux_y[j+4][i]*(density1[j+4-1][i]+limiter);
					sigmam=(mass_flux_y[j+4][i])/(density1[j+4-1][i]*pre_vol[j+4-1][i]);
					diffuw=energy1[j+4-1][i]-energy1[j+4-2][i];
					diffdw=energy1[j+4][i]-energy1[j+4-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+4][i]=mass_flux_y[j+4][i]*(energy1[j+4-1][i]+limiter);

					sigmat=(vol_flux_y[j+4+1][i])/pre_vol[j+4+1-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+4+1]/vertexdy[j+4+1-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+4+1-1][i]-density1[j+4+1-2][i];
					diffdw=density1[j+4+1][i]-density1[j+4+1-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+4+1][i]=vol_flux_y[j+4+1][i]*(density1[j+4+1-1][i]+limiter);
					sigmam=(mass_flux_y[j+4+1][i])/(density1[j+4+1-1][i]*pre_vol[j+4+1-1][i]);
					diffuw=energy1[j+4+1-1][i]-energy1[j+4+1-2][i];
					diffdw=energy1[j+4+1][i]-energy1[j+4+1-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+4+1][i]=mass_flux_y[j+4+1][i]*(energy1[j+4+1-1][i]+limiter);

					sigmat=(vol_flux_y[j+4+2][i])/pre_vol[j+4+2-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+4+2]/vertexdy[j+4+2-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+4+2-1][i]-density1[j+4+2-2][i];
					diffdw=density1[j+4+2][i]-density1[j+4+2-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+4+2][i]=vol_flux_y[j+4+2][i]*(density1[j+4+2-1][i]+limiter);
					sigmam=(mass_flux_y[j+4+2][i])/(density1[j+4+2-1][i]*pre_vol[j+4+2-1][i]);
					diffuw=energy1[j+4+2-1][i]-energy1[j+4+2-2][i];
					diffdw=energy1[j+4+2][i]-energy1[j+4+2-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+4+2][i]=mass_flux_y[j+4+2][i]*(energy1[j+4+2-1][i]+limiter);

					sigmat=(vol_flux_y[j+4+3][i])/pre_vol[j+4+3-1][i];
					sigma3=(1.0+sigmat)*(vertexdy[j+4+3]/vertexdy[j+4+3-1]);
					sigma4=2.0-sigmat;
					diffuw=density1[j+4+3-1][i]-density1[j+4+3-2][i];
					diffdw=density1[j+4+3][i]-density1[j+4+3-1][i];
					limiter=(1.0-sigmat)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					mass_flux_y[j+4+3][i]=vol_flux_y[j+4+3][i]*(density1[j+4+3-1][i]+limiter);
					sigmam=(mass_flux_y[j+4+3][i])/(density1[j+4+3-1][i]*pre_vol[j+4+3-1][i]);
					diffuw=energy1[j+4+3-1][i]-energy1[j+4+3-2][i];
					diffdw=energy1[j+4+3][i]-energy1[j+4+3-1][i];
					limiter=(1.0-sigmam)*diffuw*diffdw*0.1667*(sigma3*diffuw+sigma4*diffdw);
					ener_flux[j+4+3][i]=mass_flux_y[j+4+3][i]*(energy1[j+4+3-1][i]+limiter);
				}
			}

#pragma omp for private(i)
			for (j = 2; j < N-1; j+=8) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					double pre_mass_s=density1[j][i]*pre_vol[j][i];
					double post_mass_s=pre_mass_s+mass_flux_y[j][i]-mass_flux_y[j+1][i];
					double post_ener_s=(energy1[j][i]*pre_mass_s+ener_flux[j][i]-ener_flux[j+1][i])/post_mass_s;
					double advec_vol_s=pre_vol[j][i]+vol_flux_y[j][i]-vol_flux_y[j+1][i];
					density1[j][i]=post_mass_s/advec_vol_s;
					energy1[j][i]=post_ener_s;

					pre_mass_s=density1[j+1][i]*pre_vol[j+1][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+1][i]-mass_flux_y[j+1+1][i];
					post_ener_s=(energy1[j+1][i]*pre_mass_s+ener_flux[j+1][i]-ener_flux[j+1+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+1][i]+vol_flux_y[j+1][i]-vol_flux_y[j+1+1][i];
					density1[j+1][i]=post_mass_s/advec_vol_s;
					energy1[j+1][i]=post_ener_s;

					pre_mass_s=density1[j+2][i]*pre_vol[j+2][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+2][i]-mass_flux_y[j+2+1][i];
					post_ener_s=(energy1[j+2][i]*pre_mass_s+ener_flux[j+2][i]-ener_flux[j+2+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+2][i]+vol_flux_y[j+2][i]-vol_flux_y[j+2+1][i];
					density1[j+2][i]=post_mass_s/advec_vol_s;
					energy1[j+2][i]=post_ener_s;

					pre_mass_s=density1[j+3][i]*pre_vol[j+3][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+3][i]-mass_flux_y[j+3+1][i];
					post_ener_s=(energy1[j+3][i]*pre_mass_s+ener_flux[j+3][i]-ener_flux[j+3+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+3][i]+vol_flux_y[j+3][i]-vol_flux_y[j+3+1][i];
					density1[j+3][i]=post_mass_s/advec_vol_s;
					energy1[j+3][i]=post_ener_s;

					pre_mass_s=density1[j+4][i]*pre_vol[j+4][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+4][i]-mass_flux_y[j+4+1][i];
					post_ener_s=(energy1[j+4][i]*pre_mass_s+ener_flux[j+4][i]-ener_flux[j+4+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+4][i]+vol_flux_y[j+4][i]-vol_flux_y[j+4+1][i];
					density1[j+4][i]=post_mass_s/advec_vol_s;
					energy1[j+4][i]=post_ener_s;

					pre_mass_s=density1[j+4+1][i]*pre_vol[j+4+1][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+4+1][i]-mass_flux_y[j+4+1+1][i];
					post_ener_s=(energy1[j+4+1][i]*pre_mass_s+ener_flux[j+4+1][i]-ener_flux[j+4+1+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+4+1][i]+vol_flux_y[j+4+1][i]-vol_flux_y[j+4+1+1][i];
					density1[j+4+1][i]=post_mass_s/advec_vol_s;
					energy1[j+4+1][i]=post_ener_s;

					pre_mass_s=density1[j+4+2][i]*pre_vol[j+4+2][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+4+2][i]-mass_flux_y[j+4+2+1][i];
					post_ener_s=(energy1[j+4+2][i]*pre_mass_s+ener_flux[j+4+2][i]-ener_flux[j+4+2+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+4+2][i]+vol_flux_y[j+4+2][i]-vol_flux_y[j+4+2+1][i];
					density1[j+4+2][i]=post_mass_s/advec_vol_s;
					energy1[j+4+2][i]=post_ener_s;

					pre_mass_s=density1[j+4+3][i]*pre_vol[j+4+3][i];
					post_mass_s=pre_mass_s+mass_flux_y[j+4+3][i]-mass_flux_y[j+4+3+1][i];
					post_ener_s=(energy1[j+4+3][i]*pre_mass_s+ener_flux[j+4+3][i]-ener_flux[j+4+3+1][i])/post_mass_s;
					advec_vol_s=pre_vol[j+4+3][i]+vol_flux_y[j+4+3][i]-vol_flux_y[j+4+3+1][i];
					density1[j+4+3][i]=post_mass_s/advec_vol_s;
					energy1[j+4+3][i]=post_ener_s;
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("unroll: %6lf\n", (double)4096*4096*47*5/(end_time - start_time)/1e9);
}
