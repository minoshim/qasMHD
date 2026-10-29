#ifndef _CLASS_MHD1D_
#define _CLASS_MHD1D_

#include <string>

#include "mhd_class.hpp"
#include "mymacros.hpp"

static_assert(RMN >= 0 && RMN <= 3,
	      "RMN must be between 0 and 3");
static_assert(ODR >= 1 && ODR <= 4,
	      "ODR must be between 1 and 4");
static_assert(R_K >= 1 && R_K <= 3,
	      "R_K must be between 1 and 3");

class MHD1D : public MHD{

public:
  const int xoff=4;		// Number of ghost cells in each side
  const int nx=XMESH+2*xoff;	// Number of cells in whole domain (including offset)
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTREC;	// Time step for output
  const double tmax=dtrec*nout;	// Maximum simulation time
  const double cfl=CFL;		// CFL value
  void setdt(int);		// Set time step dt
  void paras();			// Set parameters
  void init_();			// Set initial condition
  void exec_(int);		// Run simulation
  MHD1D();			// Constructor
  virtual ~MHD1D();		// Destructor
  double getdt()		// Get dt
  {
    return dt;
  }
  
protected:
  double xmin,xmax,dx,dt;	// Left/rightmost x values, grid size, time step
  std::string fildir="./dat/";	// Directory for output
  int cnt=0,n=0;		// Counters
  bool dt_initialized=false;	// Whether output intervals have been initialized
  int nrec=1;			// Step for output (temporary value. This will be initialized in setdt)
  int nmax=nrec*nout; 		// Maximum step (This will be initialized in setdt)
  double tim=0.0;		// Simulation time
  double trec=dtrec;		// Time for next record
  void setup_grid(double, double); // Set computational domain and grid
  void bound(double *val[], int nm, const int dnxs[]); // Set boundary condition
  void ideal(double);		// ideal MHD solver
  void dout_();			// Output data
};

#endif
