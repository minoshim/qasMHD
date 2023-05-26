#ifndef _MYMACROS_
#define _MYMACROS_

#define XMESH (64)		// Number of cells in X domain
#define YMESH (64)		// Number of cells in Y domain
#define ZMESH (256)		// Number of cells in Z domain
#define MNP_X (2)		// Number of MPI processes in X
#define MNP_Y (2)		// Number of MPI processes in Y
#define MNP_Z (8)		// Number of MPI processes in Z
#define N_OUT (50)		// Number of output
#define DTREC (1.0)		// Time step for output
#define CFL (0.4)		// CFL value

// Select solvers.
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (4)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (3)		 /* Temporal order (1,2,3). Never set >3 */
#define CTW (1)		 /* Flag for CT 2D upwind weighting (Minoshima+19, ApJS,242,14) */

// Flag for random perturbation (see mhd3d_init_.cpp)
#define RANDOM (0)

#endif
