#ifndef _MHD_FD3D_H_
#define _MHD_FD3D_H_

/* [IMPORTANT] Select Riemann solver, spatial and temporal accuracy, CUCT flag */
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (2)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (3)		 /* Temporal order (1,2,3). Never set >3 */
#define CTW (1)		 /* Flag for CT 2D upwind weighting (Minoshima+19, ApJS,242,14) */

/* 3D MHD solver */
void mhd_fd3d(double *p[], double dt, double dx, double dy, double dz,
	      int nm, int nx, int ny, int nz, int xoff, int yoff, int zoff, double gamma,
	      int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz);

#endif
