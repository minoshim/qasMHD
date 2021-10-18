#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "myfunc.h"

#define XMESH (128)		// Number of cells in X domain
#define YMESH (512)		// Number of cells in Y domain
#define MNP_X (1)		// Number of MPI processes in X
#define MNP_Y (4)		// Number of MPI processes in Y
#define N_OUT (50)		// Number of output
#define DTOUT (0.5)		// Time step for output

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
  const int nx=XMESH/MNP_X+2*xoff,ny=YMESH/MNP_Y+2*yoff; // Number of cells in MPI domain (including offset)
  const int nd=nx*ny;
  const int mpi_numx=MNP_X;
  const int mpi_numy=MNP_Y;
  const int mnp=MNP_X*MNP_Y;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTOUT;	// Time step for output
  const double tend=dtrec*nout;	// Simulation end time
  const double pi=4.0*atan(1.0);
  const double dtor=pi/180.;

  // Parameters
  const double lx=20.0;	// Spatial domain in X
  const double ly=lx*(double)YMESH/XMESH;	// Spatial domain in Y
  const double xmin=0;		// Leftmost x value
  const double ymin=-0.5*ly;		// Leftmost y value
  const double dx=lx/XMESH;	// Spatial width in X
  const double dy=ly/YMESH;	// Spatial width in Y
  const double idx=1.0/dx;
  const double idy=1.0/dy;
  const double gam=1.2;	// Specific heat ratio
  const double g0=1.0;		// Gravitational acceleration
  const double dr=((dx < dy)?dx:dy);
  const double cfl=0.4;		// CFL number
  const char fildir[]="dat/";	// Directory for file output

  // Initial condition parameters
  const double beta=1e3;	// Ambient plasma beta
  const double angle_u=90;	// B field angle in upper domain. 90deg: B=Bz, 0deg: B=Bx
  const double angle_l=angle_u;	// B field angle in lower domain.
  const double vamp=0.0;	// Shear velocity amplitude
  const int nmode=2;		// Number of mode for perturbation
  const double wlen=lx/nmode;
  const double lambda=1.0;	// Shear layer width
  const double s0=0.35*ly;	// Shear position
  const double ro_u=1.0;	// Density in upper domain
  const double ro_l=0.2;	// Density in lower domain
  const double b0=1.0;		// B field strength
  const double dv=0.01;		// Perturbation amplitude
  const double lg=8*dy;		// Width of boundary layer around y=0 (for gravity profile)
  
  // Tentative parameters
  double cf=sqrt((1.0+0.5*gam*beta)*b0*b0/fmin(ro_u,ro_l));
  double vmax=vamp+cf;	// Rough estimate of vfast+vbulk
  double dt=cfl*dr/vmax;	// Time step
  int nmax=(int)(tend/dt+0.5);	// Number of maximum iteration
  
  // Variables
  double *x,*y;
  double *ro,*mx,*my,*mz,*bx,*by,*bz,*en;
  double *phi_g;	// Gravitational potential
}
