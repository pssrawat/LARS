#include "../common/common.h" 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern void viscosity_opt (double* viscosity,double* density, double* xvel0, double* yvel0, double* xarea, double* yarea, double* pressure, double *celldx, double *celldy, int N);

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
	double (*density)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yvel0)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*xarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*yarea)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double (*pressure)[4098] = (double (*)[4098]) getRandom2DArray (4098, 4098);
	double *celldx = (double *) getRandom1DArray (4098);
	double *celldy = (double *) getRandom1DArray (4098);
	double (*viscosity_out)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);
	double (*viscosity)[4098] = (double (*)[4098]) getZero2DArray (4098, 4098);

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
					viscosity[j][i]=2.0*density[j][i]*grad2*limiter*limiter;
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("orig: %6lf\n", (double)4096*4096*58*5/(end_time - start_time)/1e9);

	viscosity_opt ((double*)viscosity_out,(double*)density, (double*)xvel0, (double*)yvel0, (double*)xarea, (double*)yarea, (double*)pressure, celldx, celldy, N);

	double error = checkError2D (N, 0, (double*)viscosity_out, (double*)viscosity, 1, N-1, 1, N-1);
	if (error > TOLERANCE)
		printf("error %e\n",error);
}
