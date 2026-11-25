#ifndef _MYMACROS_
#define _MYMACROS_

#define XMESH (256)		// Number of cells in X domain
#define YMESH (256)		// Number of cells in Y domain
#define N_OUT (50)		// Number of output
#define DTREC (4.0)		// Time step for output
#define CFL (0.4)		// CFL value

// Select solvers.
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (2)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (3)		 /* Temporal order (1,2,3). Never set >3 */

// MUSCL limiter for ODR=2. Other values fall back to Minmod.
#define MUSCL_LIMITER (0) /* 0=Minmod, 1=MC, 2=van Leer, 3=Koren */

#define CTW (1)		 /* Flag for CT 2D upwind weighting (Minoshima+19, ApJS,242,14) */

// Flag for random perturbation (see mhd2d_init_.cpp)
#define RANDOM (0)

#endif
