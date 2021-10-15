#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (64)		// Number of cells in X domain
#define YMESH (64)		// Number of cells in Y domain
#define ZMESH (4)		// Number of cells in Z domain
#define MNP_X (2)		// Number of MPI processes in X
#define MNP_Y (2)		// Number of MPI processes in Y
#define MNP_Z (1)		// Number of MPI processes in Z
#define N_OUT (5)		// Number of output
#define DTOUT (0.2*3.1416)	// Time step for output

int dnxs[8]={0,0,0,0,0,0,0,0};
int dnys[8]={0,0,0,0,0,0,0,0};
int dnzs[8]={0,0,0,0,0,0,0,0};
// Boundary condition flag for ro,mx,my,mz,bx,by,bz,en (be sure of variable order)
// 0=Periodic, +1=Neumann, -1=Dirichlet
int stxs[8]={0,0,0,0,1,0,0,0};
int stys[8]={0,0,0,0,0,1,0,0};
int stzs[8]={0,0,0,0,0,0,1,0};
// Staggered grid flag for ro,mx,my,mz,bx,by,bz,en (be sure of variable order)
// Do NOT change

namespace global
{
  // Universal parameters (fixed)
  const int xoff=4,yoff=xoff,zoff=xoff;	// Number of ghost cells in each side
  const int nx=XMESH/MNP_X+2*xoff,ny=YMESH/MNP_Y+2*yoff,nz=ZMESH/MNP_Z+2*zoff; // Number of cells in MPI domain (including offset)
  const int nd=nx*ny*nz;
  const int mpi_numx=MNP_X;
  const int mpi_numy=MNP_Y;
  const int mpi_numz=MNP_Z;
  const int mnp=MNP_X*MNP_Y*MNP_Z;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTOUT;	// Time step for output
  const double tend=dtrec*nout;	// Simulation end time
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  const double lx=2.0*pi;	// Spatial domain in X
  const double ly=2.0*pi;	// Spatial domain in Y
  const double lz=2.0*pi;	// Spatial domain in Z
  const double xmin=0;		// Leftmost x value
  const double ymin=0;		// Leftmost y value
  const double zmin=0;		// Leftmost z value
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double dz=lz/ZMESH;	// Spatial width in Z
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double idz=1.0/dz;
  const double gam=5.0/3.0;	// Specific heat ratio
  const double dr=fmin(fmin(dx,dy),dz);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Tentative parameters
  double vmax=3.0;		// Rough estimate of vfast+vbulk
  double dt=cfl*dr/vmax;	// Time step
  int nmax=(int)(tend/dt+0.5);	// Number of maximum iteration
  
  // Variables
  double *x,*y,*z;
  double *ro,*mx,*my,*mz,*bx,*by,*bz,*en;
}
