#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>
#include<math.h>

extern double rtclock (void);

void viscosity_opt (double* t_viscosity,double* t_density, double* t_xvel0, double* t_yvel0, double* t_xarea, double* t_yarea, double* t_pressure, double *celldx, double *celldy, int N) {
	double (*viscosity)[4098] = (double (*)[4098]) t_viscosity;
	double (*density)[4098] = (double (*)[4098]) t_density;
	double (*xvel0)[4098] = (double (*)[4098]) t_xvel0;
	double (*yvel0)[4098] = (double (*)[4098]) t_yvel0;
	double (*xarea)[4098] = (double (*)[4098]) t_xarea;
	double (*yarea)[4098] = (double (*)[4098]) t_yarea;
	double (*pressure)[4098] = (double (*)[4098]) t_pressure;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
				for (i = 1; i < N-1; i++) {
					double ugrad=(xvel0[j][i+1]+xvel0[j+1][i+1])-(xvel0[j][i]+xvel0[j+1][i]);

					double vgrad=(yvel0[j+1][i]+yvel0[j+1][i+1])-(yvel0[j][i]+yvel0[j][i+1]);

					double div = (celldx[i]*(ugrad)+  celldy[j]*(vgrad));

					double strain2 = 0.5*(xvel0[j+1][i] + xvel0[j+1][i+1]-xvel0[j][i]-xvel0[j][i+1])/celldy[j] + 0.5*(yvel0[j][i+1] + yvel0[j+1][i+1]-yvel0[j][i]-yvel0[j+1][i])/celldx[i];

					double pgradx=(pressure[j][i+1]-pressure[j][i-1])/(celldx[i]+celldx[i+1]);
					double pgrady=(pressure[j+1][i]-pressure[j-1][1])/(celldy[j]+celldy[j+1]);

					double pgradx2 = pgradx*pgradx;
					double pgrady2 = pgrady*pgrady;

					double limiter = ((0.5*(ugrad)/celldx[i])*pgradx2+(0.5*(vgrad)/celldy[j])*pgrady2+strain2*pgradx*pgrady) / (pgradx2+pgrady2);

					pgradx = 1.08*pgradx*div;
					pgrady = 1.08*pgrady*div;
					double pgrad = sqrt(pgradx*pgradx+pgrady*pgrady);
					double xgrad = celldx[i]*pgrad/pgradx;
					double ygrad = celldy[j]*pgrad/pgrady;
					double grad  = xgrad+ygrad;
					double grad2 = grad*grad;
					viscosity[j][i]=0.1 * 2.0*density[j][i]*grad2*limiter*limiter;
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<5; t++) {
#pragma omp parallel 
		{
#pragma omp for private(i)
			for (j = 1; j < N-1; j+=1) {
				for (i = 1; i < N-1; i+=8) {
#pragma begin stencil1 unroll j=1,i=1 print-intrinsics true acc-size 1
					ugrad=(xvel0[j][i+1]+xvel0[j+1][i+1])-(xvel0[j][i]+xvel0[j+1][i]);

					vgrad=(yvel0[j+1][i]+yvel0[j+1][i+1])-(yvel0[j][i]+yvel0[j][i+1]);

					div = (celldx[i]*(ugrad)+  celldy[j]*(vgrad));

					strain2 = 0.5*(xvel0[j+1][i] + xvel0[j+1][i+1]-xvel0[j][i]-xvel0[j][i+1])/celldy[j] + 0.5*(yvel0[j][i+1] + yvel0[j+1][i+1]-yvel0[j][i]-yvel0[j+1][i])/celldx[i];

					pgradx=(pressure[j][i+1]-pressure[j][i-1])/(celldx[i]+celldx[i+1]);
					pgrady=(pressure[j+1][i]-pressure[j-1][1])/(celldy[j]+celldy[j+1]);

					pgradx2 = pgradx*pgradx;
					pgrady2 = pgrady*pgrady;

					limiter = ((0.5*(ugrad)/celldx[i])*pgradx2+(0.5*(vgrad)/celldy[j])*pgrady2+strain2*pgradx*pgrady) / (pgradx2+pgrady2);

					pgradx = 1.08*pgradx*div;
					pgrady = 1.08*pgrady*div;
					pgrad = sqrt(pgradx*pgradx+pgrady*pgrady);
					xgrad = celldx[i]*pgrad/pgradx;
					ygrad = celldy[j]*pgrad/pgrady;
					grad  = xgrad+ygrad;
					grad2 = grad*grad;
					viscosity[j][i]=2.0*density[j][i]*grad2*limiter*limiter;
#pragma end stencil1
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)4096*4096*58*5/(end_time - start_time)/1e9);
}
