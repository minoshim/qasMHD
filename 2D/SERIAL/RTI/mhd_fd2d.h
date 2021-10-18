#ifndef _MHD_FD2D_H_
#define _MHD_FD2D_H_

/* [IMPORTANT] Select Riemann solver, spatial and temporal accuracy, CUCT flag */
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (4)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (3)		 /* Temporal order (1,2,3). Never set >3 */
#define CTW (1)		 /* Flag for CT 2D upwind weighting (Minoshima+19, ApJS,242,14) */

/* 2D MHD solver */
void mhd_fd2d(double *p[], double dt, double dx, double dy,
	      int nm, int nx, int ny, int xoff, int yoff, double gamma,
	      const double *phi_g);

#endif
