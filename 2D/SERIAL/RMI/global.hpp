#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (64)		// Number of cells in X domain
#define YMESH (1024)		// Number of cells in Y domain
#define N_OUT (50)		// Number of output
#define DTOUT (1.0)		// Time step for output

#define MAGNET (1)	// Flag for finite B-field (set 0 for Hydro shock)
#define RANDOM (0)		// Flag for random perturbation (see init.hpp)

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={+1,+1,+1,+1,+1,+1,+1,+1};
// Boundary condition flag for ro,mx,my,mz,bx,by,bz,en (be sure of variable order)
// 0=Periodic, +1=Neumann, -1=Dirichlet
int stxs[8]={0,0,0,0,1,0,0,0};
int stys[8]={0,0,0,0,0,1,0,0};
// Staggered grid flag for ro,mx,my,mz,bx,by,bz,en (be sure of variable order)
// Do NOT change

namespace global
{
  // Universal parameters (fixed)
  const int xoff=4,yoff=xoff;	// Number of ghost cells in each side
  const int nx=XMESH+2*xoff,ny=YMESH+2*yoff; // Number of cells in whole domain (including offset)
  const int nd=nx*ny;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTOUT;	// Time step for output
  const double tend=dtrec*nout;	// Simulation end time
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  const double lx=1.0;	// Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;	// Spatial domain in Y
  const double xmin=0;		// Leftmost x value
  const double ymin=-0.5*ly;		// Leftmost y value
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Initial condition parameters
  //// Shock upstream paramters
  const double ro_1=+1.0;	// Density
  const double vx_1=0.0;	// Velocity in X
  const double vy_1=-1.0;	// Velocity in Y
  const double vz_1=0.0;	// Velocity in Z
  const double ma_u=10.0;	// Sound Mach number
  const double pr_1=(ro_1*vy_1*vy_1)/(gam*ma_u*ma_u); // Pressure
  const double beta=1e8;	// Beta
  const double b0_1=(MAGNET)*sqrt(2.0*pr_1/beta);	      // Magnetic field
  const double bx_1=b0_1;
  const double by_1=0.0;
  const double bz_1=0.0;
  const double vref=-0.6;	// Velocity of reference. 0 = shock rest frame
  //// Contact Discon parameters
  const double ro_3=ro_1*1e1;	// Density jump
  const double lambda=lx;	// Wavelength
  const double psi=0.1*lambda;	// Amplitude
  const double mmax=8;		// Number of modes (Available when RAMDOM=1)
  const double dro3=ro_3*0.1;	// Density perturbation (Available when RAMDOM=1)

  // Tentative parameters
  double cf=sqrt((gam+2.0/beta)*pr_1/ro_1);
  double vmax=fabs(vy_1)+cf;	// Rough estimate of vfast+vbulk
  double dt=cfl*dr/vmax;	// Time step
  int nmax=(int)(tend/dt+0.5);	// Number of maximum iteration
  
  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*bx,*by,*bz,*en;
}
