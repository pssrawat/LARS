#include "immintrin.h"
#include<stdio.h>
#include <stdlib.h>

extern double rtclock (void);

void diffterm_opt (double* t_difflux, double* t_q, double* t_ux, double* t_uy, double* t_uz, double* t_vx, double* t_vy, double* t_vz, double* t_wx, double* t_wy, double* t_wz, int N) {
	double (*difflux)[308+37][308+37][308+37] = (double (*)[308+37][308+37][308+37]) t_difflux;
	double (*q)[308+37][308+37][308+37] = (double (*)[308+37][308+37][308+37]) t_q;
	double (*ux)[308+37][308+37] = (double (*)[308+37][308+37]) t_ux;
	double (*uy)[308+37][308+37] = (double (*)[308+37][308+37]) t_uy;
	double (*uz)[308+37][308+37] = (double (*)[308+37][308+37]) t_uz;
	double (*wx)[308+37][308+37] = (double (*)[308+37][308+37]) t_wx;
	double (*wy)[308+37][308+37] = (double (*)[308+37][308+37]) t_wy;
	double (*wz)[308+37][308+37] = (double (*)[308+37][308+37]) t_wz;
	double (*vx)[308+37][308+37] = (double (*)[308+37][308+37]) t_vx;
	double (*vy)[308+37][308+37] = (double (*)[308+37][308+37]) t_vy;
	double (*vz)[308+37][308+37] = (double (*)[308+37][308+37]) t_vz;

	double dxinv0 = 0.01;
	double dxinv1 = 0.02;
	double dxinv2 = 0.03;
	double tauxx = 0.3131;
	double tauyy = 0.3132;
	double tauzz = 0.3133;
	double tauxy = 0.3134;
	double tauxz = 0.3135;
	double tauyz = 0.3136;
	double divu = 0.121;
	double uxx = 0.3171;
	double uyy = 0.3172;
	double uzz = 0.3173;
	double vxx = 0.3174;
	double vyy = 0.3175;
	double vzz = 0.3176;
	double wxx = 0.3177;
	double wyy = 0.3178;
	double wzz = 0.3179;
	double txx = 0.3181;
	double tyy = 0.3182;
	double tzz = 0.3183;
	double mechwork = -0.001;
	double uxy = 1.1141;
	double uxz = 1.1142;
	double vyz = 1.1143;
	double wzx = 1.1144;
	double wzy = 1.1145;
	double vyx = 1.1146;


	int t, i, j, k;
	double start_time, end_time;

	//Cold run
	for (t=0; t<1; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						ux[k][j][i] = (0.8 * (q[1][k][j][i+1] - q[1][k][j][i-1]) + -0.2 * (q[1][k][j][i+2] - q[1][k][j][i-2]) + 0.038 * (q[1][k][j][i+3] - q[1][k][j][i-3]) + -0.0035 * (q[1][k][j][i+4] - q[1][k][j][i-4])) * dxinv0;
						vx[k][j][i] = (0.8 * (q[2][k][j][i+1] - q[2][k][j][i-1]) + -0.2 * (q[2][k][j][i+2] - q[2][k][j][i-2]) + 0.038 * (q[2][k][j][i+3] - q[2][k][j][i-3]) + -0.0035 * (q[2][k][j][i+4] - q[2][k][j][i-4])) * dxinv0;
						wx[k][j][i] = (0.8 * (q[3][k][j][i+1] - q[3][k][j][i-1]) + -0.2 * (q[3][k][j][i+2] - q[3][k][j][i-2]) + 0.038 * (q[3][k][j][i+3] - q[3][k][j][i-3]) + -0.0035 * (q[3][k][j][i+4] - q[3][k][j][i-4])) * dxinv0;
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						uy[k][j][i] = (0.8 * (q[1][k][j+1][i] - q[1][k][j-1][i]) + -0.2 * (q[1][k][j+2][i] - q[1][k][j-2][i]) + 0.038 * (q[1][k][j+3][i] - q[1][k][j-3][i]) + -0.0035 * (q[1][k][j+4][i] - q[1][k][j-4][i])) * dxinv1;
						vy[k][j][i] = (0.8 * (q[2][k][j+1][i] - q[2][k][j-1][i]) + -0.2 * (q[2][k][j+2][i] - q[2][k][j-2][i]) + 0.038 * (q[2][k][j+3][i] - q[2][k][j-3][i]) + -0.0035 * (q[2][k][j+4][i] - q[2][k][j-4][i])) * dxinv1;
						wy[k][j][i] = (0.8 * (q[3][k][j+1][i] - q[3][k][j-1][i]) + -0.2 * (q[3][k][j+2][i] - q[3][k][j-2][i]) + 0.038 * (q[3][k][j+3][i] - q[3][k][j-3][i]) + -0.0035 * (q[3][k][j+4][i] - q[3][k][j-4][i])) * dxinv1;
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						uz[k][j][i] = (0.8 * (q[1][k+1][j][i] - q[1][k-1][j][i]) + -0.2 * (q[1][k+2][j][i] - q[1][k-2][j][i]) + 0.038 * (q[1][k+3][j][i] - q[1][k-3][j][i]) + -0.0035 * (q[1][k+4][j][i] - q[1][k-4][j][i])) * dxinv2;
						vz[k][j][i] = (0.8 * (q[2][k+1][j][i] - q[2][k-1][j][i]) + -0.2 * (q[2][k+2][j][i] - q[2][k-2][j][i]) + 0.038 * (q[2][k+3][j][i] - q[2][k-3][j][i]) + -0.0035 * (q[2][k+4][j][i] - q[2][k-4][j][i])) * dxinv2;
						wz[k][j][i] = (0.8 * (q[3][k+1][j][i] - q[3][k-1][j][i]) + -0.2 * (q[3][k+2][j][i] - q[3][k-2][j][i]) + 0.038 * (q[3][k+3][j][i] - q[3][k-3][j][i]) + -0.0035 * (q[3][k+4][j][i] - q[3][k-4][j][i])) * dxinv2;
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						uxx = (-2.847 * q[1][k][j][i]+1.6 * (q[1][k][j][i+1] + q[1][k][j][i-1]) + -0.2 * (q[1][k][j][i+2] + q[1][k][j][i-2]) + 0.0253 * (q[1][k][j][i+3] + q[1][k][j][i-3]) + -0.0017 * (q[1][k][j][i+4] + q[1][k][j][i-4])) * (dxinv0*dxinv0);
						uyy = (-2.847 * q[1][k][j][i]+1.6 * (q[1][k][j+1][i] + q[1][k][j-1][i]) + -0.2 * (q[1][k][j+2][i] + q[1][k][j-2][i]) + 0.0253 * (q[1][k][j+3][i] + q[1][k][j-3][i]) + -0.0017 * (q[1][k][j+4][i] + q[1][k][j-4][i])) * (dxinv1*dxinv1);
						uzz = (-2.847 * q[1][k][j][i]+1.6 * (q[1][k+1][j][i] + q[1][k-1][j][i]) + -0.2 * (q[1][k+2][j][i] + q[1][k-2][j][i]) + 0.0253 * (q[1][k+3][j][i] + q[1][k-3][j][i]) + -0.0017 * (q[1][k+4][j][i] + q[1][k-4][j][i])) * (dxinv2*dxinv2);
						vyx = (0.8 * (vy[k][j][i+1] - vy[k][j][i-1]) + -0.2 * (vy[k][j][i+2] - vy[k][j][i-2]) + 0.038 * (vy[k][j][i+3] - vy[k][j][i-3]) + -0.0035 * (vy[k][j][i+4] - vy[k][j][i-4])) * dxinv0;
						wzx = (0.8 * (wz[k][j][i+1] - wz[k][j][i-1]) + -0.2 * (wz[k][j][i+2] - wz[k][j][i-2]) + 0.038 * (wz[k][j][i+3] - wz[k][j][i-3]) + -0.0035 * (wz[k][j][i+4] - wz[k][j][i-4])) * dxinv0;
						difflux[1][k][j][i] =  0.3311 * (1.333 * uxx + uyy + uzz + 0.333 * (vyx + wzx));
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						vxx = (-2.847 * q[2][k][j][i]+1.6 * (q[2][k][j][i+1] + q[2][k][j][i-1]) + -0.2 * (q[2][k][j][i+2] + q[2][k][j][i-2]) + 0.0253 * (q[2][k][j][i+3] + q[2][k][j][i-3]) + -0.0017 * (q[2][k][j][i+4] + q[2][k][j][i-4])) * (dxinv0*dxinv0);
						vyy = (-2.847 * q[2][k][j][i]+1.6 * (q[2][k][j+1][i] + q[2][k][j-1][i]) + -0.2 * (q[2][k][j+2][i] + q[2][k][j-2][i]) + 0.0253 * (q[2][k][j+3][i] + q[2][k][j-3][i]) + -0.0017 * (q[2][k][j+4][i] + q[2][k][j-4][i])) * (dxinv1*dxinv1);
						vzz = (-2.847 * q[2][k][j][i]+1.6 * (q[2][k+1][j][i] + q[2][k-1][j][i]) + -0.2 * (q[2][k+2][j][i] + q[2][k-2][j][i]) + 0.0253 * (q[2][k+3][j][i] + q[2][k-3][j][i]) + -0.0017 * (q[2][k+4][j][i] + q[2][k-4][j][i])) * (dxinv2*dxinv2);
						uxy = (0.8 * (ux[k][j+1][i] - ux[k][j-1][i]) + -0.2 * (ux[k][j+2][i] - ux[k][j-2][i]) + 0.038 * (ux[k][j+3][i] - ux[k][j-3][i]) + -0.0035 * (ux[k][j+4][i] - ux[k][j-4][i])) * dxinv1;
						wzy = (0.8 * (wz[k][j+1][i] - wz[k][j-1][i]) + -0.2 * (wz[k][j+2][i] - wz[k][j-2][i]) + 0.038 * (wz[k][j+3][i] - wz[k][j-3][i]) + -0.0035 * (wz[k][j+4][i] - wz[k][j-4][i])) * dxinv1;
						difflux[2][k][j][i] =  0.3311 * (vxx+1.333 * vyy + vzz + 0.333 * (uxy + wzy));
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						wxx = (-2.847 * q[3][k][j][i]+1.6 * (q[3][k][j][i+1] + q[3][k][j][i-1]) + -0.2 * (q[3][k][j][i+2] + q[3][k][j][i-2]) + 0.0253 * (q[3][k][j][i+3] + q[3][k][j][i-3]) + -0.0017 * (q[3][k][j][i+4] + q[3][k][j][i-4])) * (dxinv0*dxinv0);
						wyy = (-2.847 * q[3][k][j][i]+1.6 * (q[3][k][j+1][i] + q[3][k][j-1][i]) + -0.2 * (q[3][k][j+2][i] + q[3][k][j-2][i]) + 0.0253 * (q[3][k][j+3][i] + q[3][k][j-3][i]) + -0.0017 * (q[3][k][j+4][i] + q[3][k][j-4][i])) * (dxinv1*dxinv1);
						wzz = (-2.847 * q[3][k][j][i]+1.6 * (q[3][k+1][j][i] + q[3][k-1][j][i]) + -0.2 * (q[3][k+2][j][i] + q[3][k-2][j][i]) + 0.0253 * (q[3][k+3][j][i] + q[3][k-3][j][i]) + -0.0017 * (q[3][k+4][j][i] + q[3][k-4][j][i])) * (dxinv2*dxinv2);
						uxz = (0.8 * (ux[k+1][j][i] - ux[k-1][j][i]) + -0.2 * (ux[k+2][j][i] - ux[k-2][j][i]) + 0.038 * (ux[k+3][j][i] - ux[k-3][j][i]) + -0.0035 * (ux[k+4][j][i] - ux[k-4][j][i])) * dxinv2;
						vyz = (0.8 * (vy[k+1][j][i] - vy[k-1][j][i]) + -0.2 * (vy[k+2][j][i] - vy[k-2][j][i]) + 0.038 * (vy[k+3][j][i] - vy[k-3][j][i]) + -0.0035 * (vy[k+4][j][i] - vy[k-4][j][i])) * dxinv2;
						difflux[3][k][j][i] =  0.3311 * (wxx + wyy+1.333 * wzz + 0.333 * (uxz + vyz));
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k++) {
				for (j = 4; j < N-4; j++) {
#pragma GCC ivdep
#pragma clang loop vectorize (enable) interleave(enable)
					for (i = 4; i < N-4; i++) {
						txx = (-2.847 * q[5][k][j][i]+1.6 * (q[5][k][j][i+1] + q[5][k][j][i-1]) + -0.2 * (q[5][k][j][i+2] + q[5][k][j][i-2]) + 0.0253 * (q[5][k][j][i+3] + q[5][k][j][i-3]) + -0.0017 * (q[5][k][j][i+4] + q[5][k][j][i-4])) * (dxinv0*dxinv0);
						tyy = (-2.847 * q[5][k][j][i]+1.6 * (q[5][k][j+1][i] + q[5][k][j-1][i]) + -0.2 * (q[5][k][j+2][i] + q[5][k][j-2][i]) + 0.0253 * (q[5][k][j+3][i] + q[5][k][j-3][i]) + -0.0017 * (q[5][k][j+4][i] + q[5][k][j-4][i])) * (dxinv1*dxinv1);
						tzz = (-2.847 * q[5][k][j][i]+1.6 * (q[5][k+1][j][i] + q[5][k-1][j][i]) + -0.2 * (q[5][k+2][j][i] + q[5][k-2][j][i]) + 0.0253 * (q[5][k+3][j][i] + q[5][k-3][j][i]) + -0.0017 * (q[5][k+4][j][i] + q[5][k-4][j][i])) * (dxinv2*dxinv2);
						divu = 0.666 * (ux[k][j][i] + vy[k][j][i] + wz[k][j][i]);
						tauxx = 2e0 * ux[k][j][i] - divu;
						tauyy = 2e0 * vy[k][j][i] - divu;
						tauzz = 2e0 * wz[k][j][i] - divu;
						tauxy = uy[k][j][i] + vx[k][j][i];
						tauxz = uz[k][j][i] + wx[k][j][i];
						tauyz = vz[k][j][i] + wy[k][j][i];
						mechwork = tauxx * ux[k][j][i] + tauyy * vy[k][j][i] + tauzz * wz[k][j][i] + (tauxy*tauxy) + (tauxz*tauxz) + (tauyz*tauyz);
						mechwork =  0.3311 * mechwork + difflux[1][k][j][i] * q[1][k][j][i] + difflux[2][k][j][i] * q[2][k][j][i] + difflux[3][k][j][i] * q[3][k][j][i];
						difflux[4][k][j][i] =  0.1 * 0.7112 * (txx + tyy + tzz) + mechwork;
					}
				}
			}
		}
	}

	start_time = rtclock ();
	for (t=0; t<2; t++) {
#pragma omp parallel 
		{
#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k+=1) {
				for (j = 4; j < N-4; j++) {
					for (i = 4; i < N-4; i+=8) {
#pragma begin stencil1 unroll k=1,j=1,i=1 print-intrinsics true acc-size 1
						ux[k][j][i] = (0.8 * (q[1][k][j][i+1] - q[1][k][j][i-1]) + -0.2 * (q[1][k][j][i+2] - q[1][k][j][i-2]) + 0.038 * (q[1][k][j][i+3] - q[1][k][j][i-3]) + -0.0035 * (q[1][k][j][i+4] - q[1][k][j][i-4])) * dxinv0;
						vx[k][j][i] = (0.8 * (q[2][k][j][i+1] - q[2][k][j][i-1]) + -0.2 * (q[2][k][j][i+2] - q[2][k][j][i-2]) + 0.038 * (q[2][k][j][i+3] - q[2][k][j][i-3]) + -0.0035 * (q[2][k][j][i+4] - q[2][k][j][i-4])) * dxinv0;
						wx[k][j][i] = (0.8 * (q[3][k][j][i+1] - q[3][k][j][i-1]) + -0.2 * (q[3][k][j][i+2] - q[3][k][j][i-2]) + 0.038 * (q[3][k][j][i+3] - q[3][k][j][i-3]) + -0.0035 * (q[3][k][j][i+4] - q[3][k][j][i-4])) * dxinv0;
						uy[k][j][i] = (0.8 * (q[1][k][j+1][i] - q[1][k][j-1][i]) + -0.2 * (q[1][k][j+2][i] - q[1][k][j-2][i]) + 0.038 * (q[1][k][j+3][i] - q[1][k][j-3][i]) + -0.0035 * (q[1][k][j+4][i] - q[1][k][j-4][i])) * dxinv1;
						vy[k][j][i] = (0.8 * (q[2][k][j+1][i] - q[2][k][j-1][i]) + -0.2 * (q[2][k][j+2][i] - q[2][k][j-2][i]) + 0.038 * (q[2][k][j+3][i] - q[2][k][j-3][i]) + -0.0035 * (q[2][k][j+4][i] - q[2][k][j-4][i])) * dxinv1;
						wy[k][j][i] = (0.8 * (q[3][k][j+1][i] - q[3][k][j-1][i]) + -0.2 * (q[3][k][j+2][i] - q[3][k][j-2][i]) + 0.038 * (q[3][k][j+3][i] - q[3][k][j-3][i]) + -0.0035 * (q[3][k][j+4][i] - q[3][k][j-4][i])) * dxinv1;
						uz[k][j][i] = (0.8 * (q[1][k+1][j][i] - q[1][k-1][j][i]) + -0.2 * (q[1][k+2][j][i] - q[1][k-2][j][i]) + 0.038 * (q[1][k+3][j][i] - q[1][k-3][j][i]) + -0.0035 * (q[1][k+4][j][i] - q[1][k-4][j][i])) * dxinv2;
						vz[k][j][i] = (0.8 * (q[2][k+1][j][i] - q[2][k-1][j][i]) + -0.2 * (q[2][k+2][j][i] - q[2][k-2][j][i]) + 0.038 * (q[2][k+3][j][i] - q[2][k-3][j][i]) + -0.0035 * (q[2][k+4][j][i] - q[2][k-4][j][i])) * dxinv2;
						wz[k][j][i] = (0.8 * (q[3][k+1][j][i] - q[3][k-1][j][i]) + -0.2 * (q[3][k+2][j][i] - q[3][k-2][j][i]) + 0.038 * (q[3][k+3][j][i] - q[3][k-3][j][i]) + -0.0035 * (q[3][k+4][j][i] - q[3][k-4][j][i])) * dxinv2;
#pragma end stencil1
					}
				}
			}

#pragma omp for private(j,i) 
			for (k = 4; k < N-4; k+=1) {
				for (j = 4; j < N-4; j++) {
					for (i = 4; i < N-4; i+=8) {
#pragma begin stencil2 unroll k=1,j=1,i=1 print-intrinsics true acc-size 1	
					uxx = (-2.847 * q[1][k][j][i]+1.6 * (q[1][k][j][i+1] + q[1][k][j][i-1]) + -0.2 * (q[1][k][j][i+2] + q[1][k][j][i-2]) + 0.0253 * (q[1][k][j][i+3] + q[1][k][j][i-3]) + -0.0017 * (q[1][k][j][i+4] + q[1][k][j][i-4])) * (dxinv0*dxinv0);
						uyy = (-2.847 * q[1][k][j][i]+1.6 * (q[1][k][j+1][i] + q[1][k][j-1][i]) + -0.2 * (q[1][k][j+2][i] + q[1][k][j-2][i]) + 0.0253 * (q[1][k][j+3][i] + q[1][k][j-3][i]) + -0.0017 * (q[1][k][j+4][i] + q[1][k][j-4][i])) * (dxinv1*dxinv1);
						uzz = (-2.847 * q[1][k][j][i]+1.6 * (q[1][k+1][j][i] + q[1][k-1][j][i]) + -0.2 * (q[1][k+2][j][i] + q[1][k-2][j][i]) + 0.0253 * (q[1][k+3][j][i] + q[1][k-3][j][i]) + -0.0017 * (q[1][k+4][j][i] + q[1][k-4][j][i])) * (dxinv2*dxinv2);
						vyx = (0.8 * (vy[k][j][i+1] - vy[k][j][i-1]) + -0.2 * (vy[k][j][i+2] - vy[k][j][i-2]) + 0.038 * (vy[k][j][i+3] - vy[k][j][i-3]) + -0.0035 * (vy[k][j][i+4] - vy[k][j][i-4])) * dxinv0;
						wzx = (0.8 * (wz[k][j][i+1] - wz[k][j][i-1]) + -0.2 * (wz[k][j][i+2] - wz[k][j][i-2]) + 0.038 * (wz[k][j][i+3] - wz[k][j][i-3]) + -0.0035 * (wz[k][j][i+4] - wz[k][j][i-4])) * dxinv0;
						difflux[1][k][j][i] =  0.3311 * (1.333 * uxx + uyy + uzz + 0.333 * (vyx + wzx));
	
						vxx = (-2.847 * q[2][k][j][i]+1.6 * (q[2][k][j][i+1] + q[2][k][j][i-1]) + -0.2 * (q[2][k][j][i+2] + q[2][k][j][i-2]) + 0.0253 * (q[2][k][j][i+3] + q[2][k][j][i-3]) + -0.0017 * (q[2][k][j][i+4] + q[2][k][j][i-4])) * (dxinv0*dxinv0);
						vyy = (-2.847 * q[2][k][j][i]+1.6 * (q[2][k][j+1][i] + q[2][k][j-1][i]) + -0.2 * (q[2][k][j+2][i] + q[2][k][j-2][i]) + 0.0253 * (q[2][k][j+3][i] + q[2][k][j-3][i]) + -0.0017 * (q[2][k][j+4][i] + q[2][k][j-4][i])) * (dxinv1*dxinv1);
						vzz = (-2.847 * q[2][k][j][i]+1.6 * (q[2][k+1][j][i] + q[2][k-1][j][i]) + -0.2 * (q[2][k+2][j][i] + q[2][k-2][j][i]) + 0.0253 * (q[2][k+3][j][i] + q[2][k-3][j][i]) + -0.0017 * (q[2][k+4][j][i] + q[2][k-4][j][i])) * (dxinv2*dxinv2);
						uxy = (0.8 * (ux[k][j+1][i] - ux[k][j-1][i]) + -0.2 * (ux[k][j+2][i] - ux[k][j-2][i]) + 0.038 * (ux[k][j+3][i] - ux[k][j-3][i]) + -0.0035 * (ux[k][j+4][i] - ux[k][j-4][i])) * dxinv1;
						wzy = (0.8 * (wz[k][j+1][i] - wz[k][j-1][i]) + -0.2 * (wz[k][j+2][i] - wz[k][j-2][i]) + 0.038 * (wz[k][j+3][i] - wz[k][j-3][i]) + -0.0035 * (wz[k][j+4][i] - wz[k][j-4][i])) * dxinv1;
						difflux[2][k][j][i] =  0.3311 * (vxx+1.333 * vyy + vzz + 0.333 * (uxy + wzy));

						wxx = (-2.847 * q[3][k][j][i]+1.6 * (q[3][k][j][i+1] + q[3][k][j][i-1]) + -0.2 * (q[3][k][j][i+2] + q[3][k][j][i-2]) + 0.0253 * (q[3][k][j][i+3] + q[3][k][j][i-3]) + -0.0017 * (q[3][k][j][i+4] + q[3][k][j][i-4])) * (dxinv0*dxinv0);
						wyy = (-2.847 * q[3][k][j][i]+1.6 * (q[3][k][j+1][i] + q[3][k][j-1][i]) + -0.2 * (q[3][k][j+2][i] + q[3][k][j-2][i]) + 0.0253 * (q[3][k][j+3][i] + q[3][k][j-3][i]) + -0.0017 * (q[3][k][j+4][i] + q[3][k][j-4][i])) * (dxinv1*dxinv1);
						wzz = (-2.847 * q[3][k][j][i]+1.6 * (q[3][k+1][j][i] + q[3][k-1][j][i]) + -0.2 * (q[3][k+2][j][i] + q[3][k-2][j][i]) + 0.0253 * (q[3][k+3][j][i] + q[3][k-3][j][i]) + -0.0017 * (q[3][k+4][j][i] + q[3][k-4][j][i])) * (dxinv2*dxinv2);
						uxz = (0.8 * (ux[k+1][j][i] - ux[k-1][j][i]) + -0.2 * (ux[k+2][j][i] - ux[k-2][j][i]) + 0.038 * (ux[k+3][j][i] - ux[k-3][j][i]) + -0.0035 * (ux[k+4][j][i] - ux[k-4][j][i])) * dxinv2;
						vyz = (0.8 * (vy[k+1][j][i] - vy[k-1][j][i]) + -0.2 * (vy[k+2][j][i] - vy[k-2][j][i]) + 0.038 * (vy[k+3][j][i] - vy[k-3][j][i]) + -0.0035 * (vy[k+4][j][i] - vy[k-4][j][i])) * dxinv2;
						difflux[3][k][j][i] =  0.3311 * (wxx + wyy+1.333 * wzz + 0.333 * (uxz + vyz));

						txx = (-2.847 * q[5][k][j][i]+1.6 * (q[5][k][j][i+1] + q[5][k][j][i-1]) + -0.2 * (q[5][k][j][i+2] + q[5][k][j][i-2]) + 0.0253 * (q[5][k][j][i+3] + q[5][k][j][i-3]) + -0.0017 * (q[5][k][j][i+4] + q[5][k][j][i-4])) * (dxinv0*dxinv0);
						tyy = (-2.847 * q[5][k][j][i]+1.6 * (q[5][k][j+1][i] + q[5][k][j-1][i]) + -0.2 * (q[5][k][j+2][i] + q[5][k][j-2][i]) + 0.0253 * (q[5][k][j+3][i] + q[5][k][j-3][i]) + -0.0017 * (q[5][k][j+4][i] + q[5][k][j-4][i])) * (dxinv1*dxinv1);
						tzz = (-2.847 * q[5][k][j][i]+1.6 * (q[5][k+1][j][i] + q[5][k-1][j][i]) + -0.2 * (q[5][k+2][j][i] + q[5][k-2][j][i]) + 0.0253 * (q[5][k+3][j][i] + q[5][k-3][j][i]) + -0.0017 * (q[5][k+4][j][i] + q[5][k-4][j][i])) * (dxinv2*dxinv2);
						divu = 0.666 * (ux[k][j][i] + vy[k][j][i] + wz[k][j][i]);
						tauxx = 2e0 * ux[k][j][i] - divu;
						tauyy = 2e0 * vy[k][j][i] - divu;
						tauzz = 2e0 * wz[k][j][i] - divu;
						tauxy = uy[k][j][i] + vx[k][j][i];
						tauxz = uz[k][j][i] + wx[k][j][i];
						tauyz = vz[k][j][i] + wy[k][j][i];
						mechwork = tauxx * ux[k][j][i] + tauyy * vy[k][j][i] + tauzz * wz[k][j][i] + (tauxy*tauxy) + (tauxz*tauxz) + (tauyz*tauyz);
						mechwork =  0.3311 * mechwork + difflux[1][k][j][i] * q[1][k][j][i] + difflux[2][k][j][i] * q[2][k][j][i] + difflux[3][k][j][i] * q[3][k][j][i];
						difflux[4][k][j][i] =  0.7112 * (txx + tyy + tzz) + mechwork;
#pragma end stencil2
					}
				}
			}
		}
	}
	end_time = rtclock ();
	printf ("opt: %6lf\n", (double)300*300*300*415*2/(end_time - start_time)/1e9);

}

