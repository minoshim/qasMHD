#ifndef _MYMACROS_
#define _MYMACROS_

#define XMESH (256)	   // Number of cells in computational domain
#define N_OUT (4095)	   // Number of output
#define DTREC (0.01)	   // Time step for output
#define CFL (0.4)	   // CFL value

// Select solvers.
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (4)		 /* Spatial order (1,2,3,4) */
#define R_K (3)		 /* Temporal order (1,2,3) */

// MUSCL limiter for ODR=2. Other values fall back to Minmod.
#define MUSCL_LIMITER (0) /* 0=Minmod, 1=MC, 2=van Leer, 3=Koren */

#endif
