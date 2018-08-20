#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern double rtclock (void);

void accelerate_opt (double *t_xvel1, double *t_yvel1, double *t_zvel1, double *t_xarea, double *t_yarea, double *t_zarea, double *t_xvel0, double *t_yvel0, double *t_zvel0, double *t_pressure, double *t_viscosity, double *t_density0, double *t_volume, double dt, int N) {
	double (*xvel1)[258][258] = (double (*)[258][258]) t_xvel1;
	double (*yvel1)[258][258] = (double (*)[258][258]) t_yvel1;
	double (*zvel1)[258][258] = (double (*)[258][258]) t_zvel1;
	double (*xarea)[258][258] = (double (*)[258][258]) t_xarea;
	double (*yarea)[258][258] = (double (*)[258][258]) t_yarea;
	double (*zarea)[258][258] = (double (*)[258][258]) t_zarea;
	double (*xvel0)[258][258] = (double (*)[258][258]) t_xvel0;
	double (*yvel0)[258][258] = (double (*)[258][258]) t_yvel0;
	double (*zvel0)[258][258] = (double (*)[258][258]) t_zvel0;
	double (*pressure)[258][258] = (double (*)[258][258]) t_pressure;
	double (*viscosity)[258][258] = (double (*)[258][258]) t_viscosity;
	double (*density0)[258][258] = (double (*)[258][258]) t_density0;
	double (*volume)[258][258] = (double (*)[258][258]) t_volume;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k+=1) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						double nodal_mass =  (density0[k  ][j-1][i-1]*volume[k  ][j-1][i-1]  
								+density0[k  ][j-1][i  ]*volume[k  ][j-1][i  ]  
								+density0[k  ][j  ][i  ]*volume[k  ][j  ][i  ]  
								+density0[k  ][j  ][i-1]*volume[k  ][j  ][i-1]
								+density0[k-1][j-1][i-1]*volume[k-1][j-1][i-1]  
								+density0[k-1][j-1][i  ]*volume[k-1][j-1][i  ]  
								+density0[k-1][j  ][i  ]*volume[k-1][j  ][i  ]  
								+density0[k-1][j  ][i-1]*volume[k-1][j  ][i-1]) 
							*0.125;
						double stepbymass_s =0.25*dt/nodal_mass;


						double xvel1_s = xvel0[k][j][i] + stepbymass_s * (xarea[k][j][i]*(pressure[k][j][i] + pressure[k][j][i-1]) + xarea[k][j-1][i]*(pressure[k][j-1][i] + pressure[k][j-1][i-1]) + xarea[k-1][j][i]*(pressure[k-1][j][i] + pressure[k-1][j][i-1]) + xarea[k-1][j-1][i]*(pressure[k-1][j-1][i] + pressure[k-1][j-1][i-1]));

						double yvel1_s = yvel0[k][j][i] + stepbymass_s * (yarea[k][j][i]*(pressure[k][j][i] + pressure[k][j-1][i]) + yarea[k][j][i-1]*(pressure[k][j][i-1] + pressure[k][j-1][i-1]) + yarea[k-1][j][i]*(pressure[k-1][j][i] + pressure[k-1][j-1][i]) + yarea[k-1][j][i-1]*(pressure[k-1][j][i-1] + pressure[k-1][j-1][i-1]));

						double zvel1_s = zvel0[k][j][i] + stepbymass_s * (zarea[k][j][i]*(pressure[k][j][i] + pressure[k-1][j][i]) + zarea[k][j-1][i]*(pressure[k][j-1][i] + pressure[k-1][j-1][i]) + zarea[k][j][i-1]*(pressure[k][j][i-1] + pressure[k-1][j][i-1]) + zarea[k][j-1][i-1]*(pressure[k][j-1][i-1] + pressure[k-1][j-1][i-1]));

						xvel1[k][j][i] = 0.1 * xvel1_s + stepbymass_s * (xarea[k][j][i]*(viscosity[k][j][i] + viscosity[k][j][i-1]) + xarea[k][j-1][i]*(viscosity[k][j-1][i] + viscosity[k][j-1][i-1]) + xarea[k-1][j][i]*(viscosity[k-1][j][i] + viscosity[k-1][j][i-1]) + xarea[k-1][j-1][i]*(viscosity[k-1][j-1][i] + viscosity[k-1][j-1][i-1]));

						yvel1[k][j][i] = 0.1 * yvel1_s + stepbymass_s * (yarea[k][j][i]*(viscosity[k][j][i] + viscosity[k][j-1][i]) + yarea[k][j][i-1]*(viscosity[k][j][i-1] + viscosity[k][j-1][i-1]) + yarea[k-1][j][i]*(viscosity[k-1][j][i] + viscosity[k-1][j-1][i]) + yarea[k-1][j][i-1]*(viscosity[k-1][j][i-1] + viscosity[k-1][j-1][i-1])); 

						zvel1[k][j][i] = 0.1 * zvel1_s + stepbymass_s * (zarea[k][j][i]*(viscosity[k][j][i] + viscosity[k-1][j][i]) + zarea[k][j-1][i]*(viscosity[k][j-1][i] + viscosity[k-1][j-1][i]) + zarea[k][j][i-1]*(viscosity[k][j][i-1] + viscosity[k-1][j][i-1]) + zarea[k][j-1][i-1]*(viscosity[k][j-1][i-1] + viscosity[k-1][j-1][i-1]));
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
			for (k = 1; k < N-1; k+=4) {
				for (j = 1; j < N-1; j++) {
					for (i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll k=4,j=1,i=1 print-intrinsics true acc-size 1
						nodal_mass =  (density0[k  ][j-1][i-1]*volume[k  ][j-1][i-1]  
                   							 +density0[k  ][j-1][i  ]*volume[k  ][j-1][i  ]  
                   							 +density0[k  ][j  ][i  ]*volume[k  ][j  ][i  ]  
                   							 +density0[k  ][j  ][i-1]*volume[k  ][j  ][i-1]
                   							 +density0[k-1][j-1][i-1]*volume[k-1][j-1][i-1]  
                   							 +density0[k-1][j-1][i  ]*volume[k-1][j-1][i  ]  
                   							 +density0[k-1][j  ][i  ]*volume[k-1][j  ][i  ]  
                   							 +density0[k-1][j  ][i-1]*volume[k-1][j  ][i-1]) 
                   							 *0.125;
			        		stepbymass_s =0.25*dt/nodal_mass;
						xvel1_s = xvel0[k][j][i] + stepbymass_s * (xarea[k][j][i]*(pressure[k][j][i] + pressure[k][j][i-1]) + xarea[k][j-1][i]*(pressure[k][j-1][i] + pressure[k][j-1][i-1]) + xarea[k-1][j][i]*(pressure[k-1][j][i] + pressure[k-1][j][i-1]) + xarea[k-1][j-1][i]*(pressure[k-1][j-1][i] + pressure[k-1][j-1][i-1]));
						yvel1_s = yvel0[k][j][i] + stepbymass_s * (yarea[k][j][i]*(pressure[k][j][i] + pressure[k][j-1][i]) + yarea[k][j][i-1]*(pressure[k][j][i-1] + pressure[k][j-1][i-1]) + yarea[k-1][j][i]*(pressure[k-1][j][i] + pressure[k-1][j-1][i]) + yarea[k-1][j][i-1]*(pressure[k-1][j][i-1] + pressure[k-1][j-1][i-1]));
						zvel1_s = zvel0[k][j][i] + stepbymass_s * (zarea[k][j][i]*(pressure[k][j][i] + pressure[k-1][j][i]) + zarea[k][j-1][i]*(pressure[k][j-1][i] + pressure[k-1][j-1][i]) + zarea[k][j][i-1]*(pressure[k][j][i-1] + pressure[k-1][j][i-1]) + zarea[k][j-1][i-1]*(pressure[k][j-1][i-1] + pressure[k-1][j-1][i-1]));
						xvel1[k][j][i] = xvel1_s + stepbymass_s * (xarea[k][j][i]*(viscosity[k][j][i] + viscosity[k][j][i-1]) + xarea[k][j-1][i]*(viscosity[k][j-1][i] + viscosity[k][j-1][i-1]) + xarea[k-1][j][i]*(viscosity[k-1][j][i] + viscosity[k-1][j][i-1]) + xarea[k-1][j-1][i]*(viscosity[k-1][j-1][i] + viscosity[k-1][j-1][i-1]));
						yvel1[k][j][i] = yvel1_s + stepbymass_s * (yarea[k][j][i]*(viscosity[k][j][i] + viscosity[k][j-1][i]) + yarea[k][j][i-1]*(viscosity[k][j][i-1] + viscosity[k][j-1][i-1]) + yarea[k-1][j][i]*(viscosity[k-1][j][i] + viscosity[k-1][j-1][i]) + yarea[k-1][j][i-1]*(viscosity[k-1][j][i-1] + viscosity[k-1][j-1][i-1])); 
						zvel1[k][j][i] = zvel1_s + stepbymass_s * (zarea[k][j][i]*(viscosity[k][j][i] + viscosity[k-1][j][i]) + zarea[k][j-1][i]*(viscosity[k][j-1][i] + viscosity[k-1][j-1][i]) + zarea[k][j][i-1]*(viscosity[k][j][i-1] + viscosity[k-1][j][i-1]) + zarea[k][j-1][i-1]*(viscosity[k][j-1][i-1] + viscosity[k-1][j-1][i-1]));
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)256*256*256*57*5/(end_time - start_time)/1e9);
}
