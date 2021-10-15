#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (256)		// Number of cells in X domain
#define YMESH (128)		// Number of cells in Y domain
#define N_OUT (10)		// Number of output
#define DTOUT (0.5)		// Time step for output

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={0,0,0,0,0,0,0,0};
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
  const double lx=2.0;	// Spatial domain in X
  const double ly=1.0;	// Spatial domain in Y
  const double xmin=-lx/2;		// Leftmost x value
  const double ymin=-ly/2;		// Leftmost y value
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Initial condition parameters
  const double ro0=1.0;		// Ambient density
  const double pr0=1e0;		// Ambient pressure
  const double v0=sqrt(5.0);	// Ambient velocity
  // const double theta=0.463647609000806; // Direction of velocity
  const double theta=0.01;		// Lee13
  const double vx0=v0*cos(theta);
  const double vy0=v0*sin(theta);
  const double vz0=0.0;
  const double a0=1e-3;		// Amplitude of magnetic loop
  const double rad=0.3;		// Radius of magnetic loop

  // Tentative parameters
  double vmax=v0+sqrt(gam*pr0/ro0);		// Rough estimate of vfast+vbulk
  double dt=cfl*dr/vmax;	// Time step
  int nmax=(int)(tend/dt+0.5);	// Number of maximum iteration
  
  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*bx,*by,*bz,*en;
}
