#ifndef _MYMACROS_
#define _MYMACROS_

#define XMESH (800)	   // Number of cells in computational domain
#define N_OUT (100)	   // Number of output
#define DTREC (0.01)	   // Time step for output
#define CFL (0.4)	   // CFL value

// Select solvers.
#define RMN (2)		 /* Riemann solver (0=Roe, 1=HLLD, 2=LHLLD, 3=MLAU) */
#define ODR (1)		 /* Spatial order (1,2,3,4). Never set >4 */
#define R_K (2)		 /* Temporal order (1,2,3). Never set >3 */

// Select initial condition. See mhd1d_init_.cpp
#define NUM (1)                 
// NUM = 1: Dai & Woodward 1994 (Miyoshi & Kusano 2005, Fig. 5)
// NUM = 2: Brio & Wu 1988 (Miyoshi & Kusano 2005, Fig. 8)
// NUM = 3: Slow switch-off shock (Miyoshi & Kusano 2005, Fig. 9)
// NUM = 4: Slow switch-off rarefaction (Miyoshi & Kusano 2005, Fig. 10)
// NUM = 5: Super-fast expansion (Miyoshi & Kusano 2005, Fig. 11)

#endif
