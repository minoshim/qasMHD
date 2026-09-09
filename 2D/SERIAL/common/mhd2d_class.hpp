#ifndef _CLASS_MHD2D_
#define _CLASS_MHD2D_

#include <string>

#include "mhd_class.hpp"
#include "mymacros.hpp"

static_assert(RMN >= 0 && RMN <= 3,
	      "RMN must be between 0 and 3");
static_assert(ODR >= 1 && ODR <= 4,
	      "ODR must be between 1 and 4");
static_assert(R_K >= 1 && R_K <= 3,
	      "R_K must be between 1 and 3");

class MHD2D : public MHD{

public:
  const int xoff=4,yoff=xoff;		// Number of ghost cells in each side
  const int nx=XMESH+2*xoff,ny=YMESH+2*yoff;	// Number of cells in whole domain (including offset)
  const int nd=nx*ny;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTREC;	// Time step for output
  const double tmax=dtrec*nout;	// Maximum simulation time
  const double cfl=CFL;		// CFL value
  virtual void setdt(int);	// Set time step dt
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
  std::string fildir="./dat/";	// Directory for output
  int cnt=0,n=0;		// Counters
  bool dt_initialized=false;	// Whether output intervals have been initialized
  int nrec=1;			// Steps per output for fixed dt (initialized in setdt)
  int nmax=nout;		// Maximum step for fixed dt (initialized in setdt)
  double tim=0.0;		// Simulation time
  double trec=dtrec;		// Time for next record
  // Physical domain bounds exclude ghost cells; shifts are measured in cell widths.
  // Use 0.5 for cell-centered coordinates, or 0.0 for the existing OTvortex grid.
  void setup_grid(double xmin_value, double xmax_value,
                  double ymin_value, double ymax_value,
                  double xshift=0.5, double yshift=0.5);
  virtual void bound(double *val[], int nm,
	     const int stxs[], const int dnxs[], const int stys[], const int dnys[]); // Set boundary condition
  virtual void ideal(double);	// ideal MHD solver
  virtual void dout_(int);	// Output data
};

#endif
