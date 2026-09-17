#ifndef _CLASS_MHD2D_
#define _CLASS_MHD2D_

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "mhd_class.hpp"
#include "mympi_class.hpp"
#include "mymacros.hpp"

static_assert(RMN >= 0 && RMN <= 3, "RMN must be between 0 and 3");
static_assert(ODR >= 1 && ODR <= 4, "ODR must be between 1 and 4");
static_assert(R_K >= 1 && R_K <= 3, "R_K must be between 1 and 3");

class MHD2D : public MYMPI, public MHD{

public:
  const int xoff=4,yoff=xoff;		// Number of ghost cells in each side
  const int mpi_numx=MNP_X;		// Number of MPI processes in X
  const int mpi_numy=MNP_Y;		// Number of MPI processes in Y
  const int mnp=MNP_X*MNP_Y;			// Number of MPI processes
  const int mpi_ranx=mpi_rank%mpi_numx;
  const int mpi_rany=mpi_rank/mpi_numx;
  const int nx=(XMESH+mpi_ranx)/mpi_numx+2*xoff;
  const int ny=(YMESH+mpi_rany)/mpi_numy+2*yoff;
  const int nd=nx*ny;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTREC;	// Time step for output
  const double tmax=dtrec*nout;	// Maximum simulation time
  const double cfl=CFL;		// CFL value
  virtual void setdt(int);	// Set time step dt
  void paras();			// Set parameters
  void init_();			// Set initial condition
  void exec_(int);		// Run simulation
  MHD2D(int*, char***, int);	// Constructor
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
  bool dt_initialized=false;  // Whether output intervals have been initialized
  int nrec=1;                 // Steps per output for fixed dt
  int nmax=nout;              // Maximum step for fixed dt
  double tim=0.0;		// Simulation time
  double trec=dtrec;		// Time for next record
  // Bounds describe the GLOBAL physical domain, excluding ghost cells.
  // Shifts are in cell widths: 0.5 for cell centers, 0.0 for the OTvortex grid.
  // Local sizes and coordinate offsets retain the existing MPI decomposition.
  void setup_grid(double xmin_value, double xmax_value,
                  double ymin_value, double ymax_value,
                  double xshift=0.5, double yshift=0.5);

  // Call only outside OpenMP regions; allocation failure must stop every rank.
  void resize_work(std::vector<double>& work, std::size_t required_size)
  {
    if (work.size() == required_size) return;
    try{
      work.resize(required_size);
    } catch (...){
      std::fprintf(stderr,"Rank %d: cannot allocate MHD2D workspace (%zu doubles).\n",
                   mpi_rank,required_size);
      MPI_Abort(MPI_COMM_WORLD,MPI_ERR_NO_MEM);
      std::abort();
    }
  }

  virtual void bound(double *val[], int nm,
	     const int stxs[], const int dnxs[], const int stys[], const int dnys[]); // Set boundary condition
  virtual void ideal(double);	// ideal MHD solver
  virtual void dout_(int);	// Output data
  int bkup_(int flg);         // Load (0) or save (1); collect status across ranks
  void initialize_dt();      // Initialize output intervals for a new calculation
  void prepare_run(int flg);  // Load and validate a restart, or initialize a new run
  void check_time_parameters() const;
  void check_step() const;
  void check_restart_state() const;
};

#endif
