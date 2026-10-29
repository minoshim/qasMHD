#include "mhd3d_class.hpp"
#include "mhd3d_io.hpp"

#include <ctime>
#include <random>
#include <vector>

void MHD3D::init_()
try
{
  // KH instability
  int i,j,k;
  // Initial condition parameters
  const double beta=1e2;	// Ambient plasma beta
  const double angle_u=90.0;	// B field angle in upper domain. 90deg: B=Bz, 0deg: B=Bx
  const double angle_l=angle_u;	// B field angle in lower domain.
  const double vamp=0.5;	// Shear velocity amplitude
  const int nmode=1;		// Number of mode for perturbation
  const double wlen=getlx()/nmode;
  const double lambda=1.0;	// Shear layer width
  const double s0=ymin+0.5*getly(); // Shear position
  const double ro_u=1.0;	// Density in upper domain
  const double ro_l=1.0;	// Density in lower domain
  const double b0=1.0;		// B field strength
  const double pr0=0.5*beta*b0*b0; // Pressure
  const double dv=0.01;		// Perturbation amplitude
  const double dvparax[2]={0,dv};	// Perturbation parameter for vx
  const double dvparay[2]={0,dv};	// Perturbation paeameter for vy

  std::vector<double> dvy;
  resize_work(dvy,static_cast<std::size_t>(nz)*nx);
  // For reproducible runs, replace std::time(nullptr) with a fixed seed.
  unsigned seed=0;
  if (mpi_rank == 0) seed=static_cast<unsigned>(std::time(nullptr));
  mpi3d_io::check_mpi(MPI_Bcast(&seed,1,MPI_UNSIGNED,0,MPI_COMM_WORLD),
                      "Random seed broadcast failed.");
#if (RANDOM)
  // The same X/Z subdomain uses the same stream, independent of Y rank.
  std::seed_seq seeds_vy{seed,static_cast<unsigned>(mpi_ranx),static_cast<unsigned>(mpi_ranz)};
  std::mt19937 engine_vy(seeds_vy);
#endif
  
  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      dvy[nx*k+i]=0.0;
#if (RANDOM)
      dvy[nx*k+i]+=rand_noise_mt(dvparay,engine_vy); // Multiple mode perturbation
#else
      dvy[nx*k+i]+=dv*sin(2*M_PI*x[i]/wlen); // Single mode perturbation
#endif    
    }
  }

  // A separate Z-only stream preserves the vx perturbation even when RANDOM=0.
  // The same Z subdomain uses the same stream, independent of X/Y rank.
  std::seed_seq seeds_vx{seed,static_cast<unsigned>(mpi_ranz)};
  std::mt19937 engine_vx(seeds_vx);
  
  for (k=0;k<nz;k++){
    double vampx=vamp+rand_noise_mt(dvparax,engine_vx); // Perturbation to vx is independent of x and y
    for (j=0;j<ny;j++){
      double angle=0.5*((angle_u+angle_l)+(angle_u-angle_l)*tanh((y[j]-s0)/lambda));
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;

	ro[ss]=0.5*((ro_u+ro_l)+(ro_u-ro_l)*tanh((y[j]-s0)/lambda));
	vx[ss]=+vampx*tanh((y[j]-s0)/lambda);
	vy[ss]=dvy[nx*k+i]*exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda));
	vz[ss]=0.0;
	pr[ss]=pr0;
	
	if (angle == 90){
	  bx[ss]=0.0;
	  by[ss]=0.0;
	  bz[ss]=b0;
	} else{
	  bx[ss]=b0*cos(angle*dtor);
	  by[ss]=0.0;
	  bz[ss]=b0*sin(angle*dtor);
	}
	cx[ss]=bx[ss];
	cy[ss]=by[ss];
	cz[ss]=bz[ss];
	
	cnsvt(ss);
      }
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
} catch (...){
  mpi3d_io::abort_run("Exception while initializing KHI state or random generator.");
}
