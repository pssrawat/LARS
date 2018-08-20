#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>
#include<math.h>

extern double rtclock (void);

void viscosity_opt (double* t_viscosity,double* t_density, double* t_xvel0, double* t_yvel0, double* t_zvel0, double* t_xarea, double* t_yarea, double* t_zarea, double* t_pressure, double *celldx, double *celldy, double *celldz, int N) {
	double (*viscosity)[258][258] = (double (*)[258][258]) t_viscosity;
	double (*density)[258][258] = (double (*)[258][258]) t_density;
	double (*xvel0)[258][258] = (double (*)[258][258]) t_xvel0;
	double (*yvel0)[258][258] = (double (*)[258][258]) t_yvel0;
	double (*zvel0)[258][258] = (double (*)[258][258]) t_zvel0;
	double (*xarea)[258][258] = (double (*)[258][258]) t_xarea;
	double (*yarea)[258][258] = (double (*)[258][258]) t_yarea;
	double (*zarea)[258][258] = (double (*)[258][258]) t_zarea;
	double (*pressure)[258][258] = (double (*)[258][258]) t_pressure;

	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i)
			for (k = 1; k < N-1; k++) {
				for (j = 1; j < N-1; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 1; i < N-1; i++) {
						double ugradx1 = xvel0[k][j][i] + xvel0[k][j+1][i] + xvel0[k+1][j][i] + xvel0[k+1][j+1][i];
						double ugradx2 = xvel0[k][j][i+1] + xvel0[k][j+1][i+1] + xvel0[k+1][j][i+1] + xvel0[k+1][j+1][i+1];
						double ugrady1 = xvel0[k][j][i] + xvel0[k][j][i+1] + xvel0[k+1][j][i] + xvel0[k+1][j][i+1];
						double ugrady2 = xvel0[k][j+1][i] + xvel0[k][j+1][i+1] + xvel0[k+1][j+1][i] + xvel0[k+1][j+1][i+1];
						double ugradz1 = xvel0[k][j][i] + xvel0[k][j][i+1] + xvel0[k][j+1][i] + xvel0[k][j+1][i+1];
						double ugradz2 = xvel0[k+1][j][i] + xvel0[k+1][j][i+1] + xvel0[k+1][j+1][i] + xvel0[k+1][j+1][i+1];

						double vgradx1 = yvel0[k][j][i] + yvel0[k][j+1][i] + yvel0[k+1][j][i] + yvel0[k+1][j+1][i];
						double vgradx2 = yvel0[k][j][i+1] + yvel0[k][j+1][i+1] + yvel0[k+1][j][i+1] + yvel0[k+1][j+1][i+1];
						double vgrady1 = yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k+1][j][i] + yvel0[k+1][j][i+1];
						double vgrady2 = yvel0[k][j+1][i] + yvel0[k][j+1][i+1] + yvel0[k+1][j+1][i] + yvel0[k+1][j+1][i+1];
						double vgradz1 = yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k][j+1][i] + yvel0[k][j+1][i+1];
						double vgradz2 = yvel0[k+1][j][i] + yvel0[k+1][j][i+1] + yvel0[k+1][j+1][i] + yvel0[k+1][j+1][i+1];

						double wgradx1 = zvel0[k][j][i] + zvel0[k][j+1][i] + zvel0[k+1][j][i] + zvel0[k+1][j+1][i];
						double wgradx2 = zvel0[k][j][i+1] + zvel0[k][j+1][i+1] + zvel0[k+1][j][i+1] + zvel0[k+1][j+1][i+1];
						double wgrady1 = zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k+1][j][i] + zvel0[k+1][j][i+1];
						double wgrady2 = zvel0[k][j+1][i] + zvel0[k][j+1][i+1] + zvel0[k+1][j+1][i] + zvel0[k+1][j+1][i+1];
						double wgradz1 = zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k][j+1][i] + zvel0[k][j+1][i+1];
						double wgradz2 = zvel0[k+1][j][i] + zvel0[k+1][j][i+1] + zvel0[k+1][j+1][i] + zvel0[k+1][j+1][i+1];

						double div = (xarea[k][j][i]*(ugradx2+ugradx1) + yarea[k][j][i]*(vgrady2+vgrady1))+ zarea[k][j][i]*(wgradz2+wgradz1);

						double xx = 0.00025*(ugradx2+ugradx1)/(celldx[i]);
						double yy = 0.00025*(vgrady2+vgrady1)/(celldy[j]);
						double zz = 0.00025*(wgradz2+wgradz1)/(celldz[k]);
						double xy = 0.00025*(ugrady2+ugrady1)/(celldy[j])+0.00025*(vgradx2+vgradx1)/(celldx[i]);
						double xz = 0.00025*(ugradz2+ugradz1)/(celldz[k])+0.00025*(wgradx2+wgradx1)/(celldx[i]);
						double yz = 0.00025*(vgradz2+vgradz1)/(celldz[k])+0.00025*(wgrady2+wgrady1)/(celldy[j]);

						double pgradx = (pressure[k][j][i+1] + pressure[k][j][i-1]) / (celldx[i] + celldx[i+1]);
						double pgrady = (pressure[k][j+1][i] + pressure[k][j-1][i]) / (celldy[j] + celldy[j+1]);
						double pgradz = (pressure[k+1][j][i] + pressure[k-1][j][i]) / (celldz[k] + celldz[k+1]);

						double pgradx2 = pgradx*pgradx;
						double pgrady2 = pgrady*pgrady;
						double pgradz2 = pgradz*pgradz;

						double limiter = (xx*pgradx2+yy*pgrady2+zz*pgradz2 + xy*pgradx*pgrady+xz*pgradx*pgradz+yz*pgrady*pgradz) / (pgradx2+pgrady2+pgradz2);

						double pgrad = (pgradx*pgradx+pgrady*pgrady+pgradz*pgradz);
						double xgrad = (celldx[i]*pgrad/pgradx);
						double ygrad = (celldy[j]*pgrad/pgrady);
						double zgrad = (celldz[k]*pgrad/pgradz);
						double grad  = xgrad+ygrad+zgrad;
						double grad2 = grad*grad;
						viscosity[k][j][i] = 0.1 * 2.0*density[k][j][i]*grad2*limiter*limiter;
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
						ugradx1 = xvel0[k][j][i] + xvel0[k][j+1][i] + xvel0[k+1][j][i] + xvel0[k+1][j+1][i];
						ugradx2 = xvel0[k][j][i+1] + xvel0[k][j+1][i+1] + xvel0[k+1][j][i+1] + xvel0[k+1][j+1][i+1];
						ugrady1 = xvel0[k][j][i] + xvel0[k][j][i+1] + xvel0[k+1][j][i] + xvel0[k+1][j][i+1];
						ugrady2 = xvel0[k][j+1][i] + xvel0[k][j+1][i+1] + xvel0[k+1][j+1][i] + xvel0[k+1][j+1][i+1];
						ugradz1 = xvel0[k][j][i] + xvel0[k][j][i+1] + xvel0[k][j+1][i] + xvel0[k][j+1][i+1];
						ugradz2 = xvel0[k+1][j][i] + xvel0[k+1][j][i+1] + xvel0[k+1][j+1][i] + xvel0[k+1][j+1][i+1];

						vgradx1 = yvel0[k][j][i] + yvel0[k][j+1][i] + yvel0[k+1][j][i] + yvel0[k+1][j+1][i];
						vgradx2 = yvel0[k][j][i+1] + yvel0[k][j+1][i+1] + yvel0[k+1][j][i+1] + yvel0[k+1][j+1][i+1];
						vgrady1 = yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k+1][j][i] + yvel0[k+1][j][i+1];
						vgrady2 = yvel0[k][j+1][i] + yvel0[k][j+1][i+1] + yvel0[k+1][j+1][i] + yvel0[k+1][j+1][i+1];
						vgradz1 = yvel0[k][j][i] + yvel0[k][j][i+1] + yvel0[k][j+1][i] + yvel0[k][j+1][i+1];
						vgradz2 = yvel0[k+1][j][i] + yvel0[k+1][j][i+1] + yvel0[k+1][j+1][i] + yvel0[k+1][j+1][i+1];

						wgradx1 = zvel0[k][j][i] + zvel0[k][j+1][i] + zvel0[k+1][j][i] + zvel0[k+1][j+1][i];
						wgradx2 = zvel0[k][j][i+1] + zvel0[k][j+1][i+1] + zvel0[k+1][j][i+1] + zvel0[k+1][j+1][i+1];
						wgrady1 = zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k+1][j][i] + zvel0[k+1][j][i+1];
						wgrady2 = zvel0[k][j+1][i] + zvel0[k][j+1][i+1] + zvel0[k+1][j+1][i] + zvel0[k+1][j+1][i+1];
						wgradz1 = zvel0[k][j][i] + zvel0[k][j][i+1] + zvel0[k][j+1][i] + zvel0[k][j+1][i+1];
						wgradz2 = zvel0[k+1][j][i] + zvel0[k+1][j][i+1] + zvel0[k+1][j+1][i] + zvel0[k+1][j+1][i+1];

						div = (xarea[k][j][i]*(ugradx2+ugradx1) + yarea[k][j][i]*(vgrady2+vgrady1))+ zarea[k][j][i]*(wgradz2+wgradz1);

						xx = 0.00025*(ugradx2+ugradx1)/(celldx[i]);
						yy = 0.00025*(vgrady2+vgrady1)/(celldy[j]);
						zz = 0.00025*(wgradz2+wgradz1)/(celldz[k]);
						xy = 0.00025*(ugrady2+ugrady1)/(celldy[j])+0.00025*(vgradx2+vgradx1)/(celldx[i]);
						xz = 0.00025*(ugradz2+ugradz1)/(celldz[k])+0.00025*(wgradx2+wgradx1)/(celldx[i]);
						yz = 0.00025*(vgradz2+vgradz1)/(celldz[k])+0.00025*(wgrady2+wgrady1)/(celldy[j]);

						pgradx = (pressure[k][j][i+1] + pressure[k][j][i-1]) / (celldx[i] + celldx[i+1]);
						pgrady = (pressure[k][j+1][i] + pressure[k][j-1][i]) / (celldy[j] + celldy[j+1]);
						pgradz = (pressure[k+1][j][i] + pressure[k-1][j][i]) / (celldz[k] + celldz[k+1]);

						pgradx2 = pgradx*pgradx;
						pgrady2 = pgrady*pgrady;
						pgradz2 = pgradz*pgradz;

						limiter = (xx*pgradx2+yy*pgrady2+zz*pgradz2 + xy*pgradx*pgrady+xz*pgradx*pgradz+yz*pgrady*pgradz) / (pgradx2+pgrady2+pgradz2);

						pgrad = (pgradx*pgradx+pgrady*pgrady+pgradz*pgradz);
						xgrad = (celldx[i]*pgrad/pgradx);
						ygrad = (celldy[j]*pgrad/pgrady);
						zgrad = (celldz[k]*pgrad/pgradz);
						grad  = xgrad+ygrad+zgrad;
						grad2 = grad*grad;
						viscosity[k][j][i] = 2.0*density[k][j][i]*grad2*limiter*limiter;
#pragma end stencil1
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)256*256*256*139*5/(end_time - start_time)/1e9);
}
