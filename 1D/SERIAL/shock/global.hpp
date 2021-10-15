#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (800)		// Number of cells in computational domain
#define N_OUT (100)		// Number of output
#define DTOUT (0.01)		// Time step for output

#define NUM (1)			// Select initial condition. See init.hpp
// NUM = 1: Dai & Woodward 1994 (Miyoshi & Kusano 2005, Fig. 5)
// NUM = 2: Brio & Wu 1988 (Miyoshi & Kusano 2005, Fig. 8)
// NUM = 3: Slow switch-off shock (Miyoshi & Kusano 2005, Fig. 9)
// NUM = 4: Slow switch-off rarefaction (Miyoshi & Kusano 2005, Fig. 10)
// NUM = 5: Super-fast expansion (Miyoshi & Kusano 2005, Fig. 11)

namespace global
{
  // Universal parameters (fixed)
  const int xoff=4;		// Number of ghost cells in each side
  const int nx=XMESH+2*xoff;	// Number of cells in whole domain (including offset)
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTOUT;	// Time step for output
  const double pi=4.0*atan(1.0);

  // Parameters
  const double lx=1.0;		// Spatial domain size
  const double xmin=-0.5;	// leftmost x value
  const double dx=lx/XMESH;	// Spatial width
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Tentative parameters
  double gam=5./3.;	// Specific heat ratio
  double bx=1.0;		// Normal magnetic field strength
  double vmax=1.0;	// Maximum wave velocity
  double dt=cfl*dx/vmax;	// Time step for calculation
  int nrec=(int)(dtrec/dt+0.5);	// Number of iterations for output
  int nmax=nrec*nout;		// Number of maximum iteration 

  // Variables
  double *x,*ro,*mx,*my,*mz,*by,*bz,*en;
}
