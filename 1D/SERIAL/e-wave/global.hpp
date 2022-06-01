#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (256)		// Number of cells in computational domain
#define N_OUT (4095)		// Number of output
#define DTOUT (0.01)		// Time step for output

int dnxs[7]={0,0,0,0,0,0,0};
// Boundary condition flag for ro,mx,my,mz,by,bz,en (be sure of variable order) 
// 0=Periodic, +1=Neumann, -1=Dirichlet

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
  const double xmin=0.0;	// leftmost x value
  const double dx=lx/XMESH;	// Spatial width
  const double gam=5./3.;	// Specific heat ratio
  const double bx=1.0;		// Normal magnetic field strength
  const double beta=10.0;	// Plasma beta
  const double ro0=1.0;		// Ambient density
  const double pr0=0.5*beta*bx*bx; // Ambient pressure
  const double cf0=sqrt((bx*bx+gam*pr0)/ro0); // Magnetosonic speed
  const double vx0=0.0*cf0;	// Ambient normal velocity
  const double cfl=0.4;		// CFL number
  const double de=2.0*dx;	// Electron inertia length
  const char fildir[]="dat/";	// Directory for file output

  // Tentative parameters
  double vmax=fabs(vx0)+cf0;	// Maximum wave velocity
  double dt=cfl*dx/vmax;	// Time step for calculation
  int nrec=(int)(dtrec/dt+0.5);	// Number of iterations for output
  int nmax=nrec*nout;		// Number of maximum iteration 

  // Variables
  double *x,*ro,*mx,*my,*mz,*by,*bz,*en;
}
