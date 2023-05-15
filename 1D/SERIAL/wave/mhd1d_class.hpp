#ifndef _CLASS_MHD1D_
#define _CLASS_MHD1D_

#include "mhd_class.hpp"
#include "mymacros.hpp"

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
  char fildir[100];		// Directory for output
  int cnt=0,n=0;		// Counters
  int nrec=1;			// Step for output (temporary value. This will be initialized in setdt)
  int nmax=nrec*nout; 		// Maximum step (This will be initialized in setdt)
  double tim=0.0;		// Simulation time
  double trec=dtrec;		// Time for next record
  void bound(double *val[], int nm, const int dnxs[]); // Set boundary condition
  void ideal(double);		// ideal MHD solver
  void dout_();			// Output data
};

#endif
