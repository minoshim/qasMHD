#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (1024)		// Number of cells in X domain
#define YMESH (128)		// Number of cells in Y domain
#define N_OUT (50)		// Number of output
#define DTOUT (5.0)		// Time step for output

#define RANDOM (0)		// Flag for random perturbation (see init.hpp)
#define DIFF (1)		// Flag for diffusion

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={+1,+1,-1,+1,+1,-1,+1,+1};
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
  const double lx=64.0;	// Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;	// Spatial domain in Y
  const double xmin=-0.5*lx;		// Leftmost x value
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
  const double lambda=1.0;	// Current sheet thickness
  const double beta=0.2;	// Plasma beta @ lobe
  const double ro0=1.0;		// Density @ CS
  const double ro1=1.0;		// Density @ lobe
  const double b0=1.0;		// Mag field @ lobe
  const double b1=0.05;		// Mag field perturbation by Zenitani
  const double bg=0.0;		// Guide mag field along Z
  const double dv=0.01;		// Random noize perturbation to Vy (avaiable when RANDOM=1)
  const double eta0=DIFF*1.0*lambda/2e2; // Uniform resistivity
  const double prm=0e0;			 // Magnetic Prandtl number
  const double nu0=prm*eta0;		 // Uniform kinetic viscosity

  // Tentative parameters
  double cf=sqrt((1.0+0.5*gam*beta)*(b0*b0+bg*bg)/ro1);
  double vmax=cf;	// Rough estimate of vfast+vbulk
  double dt=cfl*dr/vmax;	// Time step
  int nmax=(int)(tend/dt+0.5);	// Number of maximum iteration
  
  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*bx,*by,*bz,*en;
  double *nu,*eta;
}
