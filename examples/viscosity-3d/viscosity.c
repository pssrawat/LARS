#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>

extern void viscosity_opt (double* viscosity,double* density, double* xvel0, double* yvel0, double* zvel0, double* xarea, double* yarea, double* zarea, double* pressure, double *celldx, double *celldy, double *celldz, int N);

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
	double (*density)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zvel0)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*xarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*yarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*zarea)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double (*pressure)[258][258] = (double (*)[258][258]) getRandom3DArray (258, 258, 258);
	double *celldx = (double *) getRandom1DArray (258);
	double *celldy = (double *) getRandom1DArray (258);
	double *celldz = (double *) getRandom1DArray (258);
	double (*viscosity_ref)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);
	double (*viscosity)[258][258] = (double (*)[258][258]) getZero3DArray (258, 258, 258);

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
						viscosity_ref[k][j][i] = 0.1 * 2.0*density[k][j][i]*grad2*limiter*limiter;
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
						viscosity_ref[k][j][i] = 2.0*density[k][j][i]*grad2*limiter*limiter;
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)256*256*256*139*5/(end_time - start_time)/1e9);

	viscosity_opt ((double*)viscosity,(double*)density, (double*)xvel0, (double*)yvel0, (double*)zvel0, (double*)xarea, (double*)yarea, (double*)zarea, (double*)pressure, celldx, celldy, celldz, N);

	double error = checkError3D (N, N, 0, (double*)viscosity_ref, (double*)viscosity, 1, N-1, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
}
