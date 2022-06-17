#ifndef _MHD_FD1D_H_
#define _MHD_FD1D_H_

/* [IMPORTANT] Select Riemann solver, spatial and temporal accuracy */
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (2)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (3)		 /* Temporal order (1,2,3). Never set >3 */

/* 1D MHD solver */
void mhd_fd1d(double *ro, double *mx, double *my, double *mz,
	      double *by, double *bz, double *en,
	      double bx, double dt, double dx, double de, double mpme,
	      int nx, int xoff, double gamma);

#endif
