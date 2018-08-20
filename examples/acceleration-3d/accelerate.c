#include "../common/common.h" 
#include "immintrin.h"
#include<stdio.h>
#include<stdlib.h>

extern void accelerate_opt (double *t_xvel1, double *t_yvel1, double *t_zvel1, double *t_xarea, double *t_yarea, double *t_zarea, double *t_xvel0, double *t_yvel0, double *t_zvel0, double *t_pressure, double *t_viscosity, double *t_density0, double *t_volume, double dt, int N); 

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
	double (*xvel1_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*yvel1_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*zvel1_opt)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*xvel1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*yvel1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*zvel1)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

	double (*xvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*pressure)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*viscosity)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*density0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*volume)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);

	double dt = 0.856;

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

						xvel1[k][j][i] = xvel1_s + stepbymass_s * (xarea[k][j][i]*(viscosity[k][j][i] + viscosity[k][j][i-1]) + xarea[k][j-1][i]*(viscosity[k][j-1][i] + viscosity[k][j-1][i-1]) + xarea[k-1][j][i]*(viscosity[k-1][j][i] + viscosity[k-1][j][i-1]) + xarea[k-1][j-1][i]*(viscosity[k-1][j-1][i] + viscosity[k-1][j-1][i-1]));

						yvel1[k][j][i] = yvel1_s + stepbymass_s * (yarea[k][j][i]*(viscosity[k][j][i] + viscosity[k][j-1][i]) + yarea[k][j][i-1]*(viscosity[k][j][i-1] + viscosity[k][j-1][i-1]) + yarea[k-1][j][i]*(viscosity[k-1][j][i] + viscosity[k-1][j-1][i]) + yarea[k-1][j][i-1]*(viscosity[k-1][j][i-1] + viscosity[k-1][j-1][i-1])); 

						zvel1[k][j][i] = zvel1_s + stepbymass_s * (zarea[k][j][i]*(viscosity[k][j][i] + viscosity[k-1][j][i]) + zarea[k][j-1][i]*(viscosity[k][j-1][i] + viscosity[k-1][j-1][i]) + zarea[k][j][i-1]*(viscosity[k][j][i-1] + viscosity[k-1][j][i-1]) + zarea[k][j-1][i-1]*(viscosity[k][j-1][i-1] + viscosity[k-1][j-1][i-1]));
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*57*5/(end_time - start_time)/1e9);

	accelerate_opt ((double*)xvel1_opt, (double*)yvel1_opt, (double*)zvel1_opt, (double*)xarea, (double*)yarea, (double*)zarea, (double*)xvel0, (double*)yvel0, (double*)zvel0, (double*)pressure, (double*)viscosity, (double*)density0, (double *)volume, dt, N);

	double error = checkError3D (N, N, 0, (double*)xvel1_opt, (double*)xvel1, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
	error = checkError3D (N, N, 0, (double*)yvel1_opt, (double*)yvel1, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
	error = checkError3D (N, N, 0, (double*)zvel1_opt, (double*)zvel1, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf ("error %e\n", error);
}
