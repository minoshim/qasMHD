#ifndef _MYMACROS_
#define _MYMACROS_

#define XMESH (100)		// Number of cells in X domain
#define YMESH (100)		// Number of cells in Y domain
#define ZMESH (200)		// Number of cells in Z domain
#define MNP_X (2)		// Number of MPI processes in X
#define MNP_Y (2)		// Number of MPI processes in Y
#define MNP_Z (4)		// Number of MPI processes in Z
#define N_OUT (6)		// Number of output
#define DTREC (0.2)		// Time step for output
#define CFL (0.4)		// CFL value

// Select solvers.
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (4)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (3)		 /* Temporal order (1,2,3). Never set >3 */
#define CTW (1)		 /* Flag for CT 2D upwind weighting (Minoshima+19, ApJS,242,14) */

#endif
