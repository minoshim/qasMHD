#ifndef _CLASS_MHD3D_
#define _CLASS_MHD3D_

#include "mhd_class.hpp"
#include "mympi_class.hpp"
#include "mymacros.hpp"

class MHD3D : public MYMPI, public MHD{

public:
  const int xoff=4,yoff=xoff,zoff=xoff;	// Number of ghost cells in each side
  const int mpi_numx=MNP_X;		// Number of MPI processes in X
  const int mpi_numy=MNP_Y;		// Number of MPI processes in Y
  const int mpi_numz=MNP_Z;		// Number of MPI processes in Z
  const int mnp=MNP_X*MNP_Y*MNP_Z;	// Number of MPI processes
  const int m_xy=MNP_X*MNP_Y;
  const int mpi_ranx=(mpi_rank%m_xy)%mpi_numx;
  const int mpi_rany=(mpi_rank%m_xy)/mpi_numx;
  const int mpi_ranz=(mpi_rank/m_xy);
  const int nx=(XMESH+mpi_ranx)/mpi_numx+2*xoff; // Number of X cells in MPI domain (including offset)
  const int ny=(YMESH+mpi_rany)/mpi_numy+2*yoff; // Number of Y cells in MPI domain (including offset)
  const int nz=(ZMESH+mpi_ranz)/mpi_numz+2*zoff; // Number of Z cells in MPI domain (including offset)
  const int nd=nx*ny*nz;
  const int nout=N_OUT;		// Number of output
  const double dtrec=DTREC;	// Time step for output
  const double tmax=dtrec*nout;	// Maximum simulation time
  const double cfl=CFL;		// CFL value
  void setdt(int);		// Set time step dt
  void paras();			// Set parameters
  void init_();			// Set initial condition
  void exec_(int);		// Run simulation
  MHD3D(int*, char***, int);	// Constructor
  virtual ~MHD3D();		// Destructor
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
  double getlz()		// Get domain length in z
  {
    return zmax-zmin;
  }

protected:
  double xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz,dr,dt; // Left/rightmost x,y,z values, grid size, time step
  char fildir[100];		// Directory for output
  int cnt=0,n=0;		// Counters
  int nmax=nout;		// Maximum step (This will be initialized in setdt)
  double tim=0.0;		// Simulation time
  double trec=dtrec;		// Time for next record
  void bound(double *val[], int nm,
	     const int stxs[], const int dnxs[], const int stys[], const int dnys[], const int stzs[], const int dnzs[]); // Set boundary condition
  void ideal(double);		// ideal MHD solver
  void dout_(int);		// Output data
  void bkup_(int flg){		// Load (flg == 0) or save (flg == 1) backup data
    if (flg == 0){
      bkup_load(val,nm,nd,&n,&cnt,&tim,&dt,&trec,mpi_rank,fildir);
    } else{
      bkup_save(val,nm,nd,n,cnt,tim,dt,trec,mpi_rank,fildir);
    }
  }
  
};

#endif
