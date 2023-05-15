#ifndef _CLASS_MHD2D_
#define _CLASS_MHD2D_

#include "mhd_class.hpp"
#include "mymacros.hpp"

class MHD2D : public MHD{

public:
  const int xoff=4,yoff=xoff;		// Number of ghost cells in each side
  const int nx=XMESH+2*xoff,ny=YMESH+2*yoff;	// Number of cells in whole domain (including offset)
  const int nd=nx*ny;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTREC;	// Time step for output
  const double tmax=dtrec*nout;	// Maximum simulation time
  const double cfl=CFL;		// CFL value
  void setdt(int);		// Set time step dt
  void paras();			// Set parameters
  void init_();			// Set initial condition
  void exec_(int);		// Run simulation
  MHD2D();			// Constructor
  virtual ~MHD2D();		// Destructor
  double getdt()		// Get dt
  {
    return dt;
  }
  double getlx()		// Get domain length in x
  {
    return xmax-xmin;
  }
  double getly()		// Get domain length in y
  {
    return ymax-ymin;
  }
  
protected:
  double xmin,xmax,dx,ymin,ymax,dy,dr,dt; // Left/rightmost x and y values, grid size, time step
  char fildir[100];		// Directory for output
  int cnt=0,n=0;		// Counters
  int nmax=nout;		// Maximum step (This will be initialized in setdt)
  double tim=0.0;		// Simulation time
  double trec=dtrec;		// Time for next record
  void bound(double *val[], int nm,
	     const int stxs[], const int dnxs[], const int stys[], const int dnys[]); // Set boundary condition
  void ideal(double);		// ideal MHD solver
  void dout_(int);		// Output data
};

#endif
