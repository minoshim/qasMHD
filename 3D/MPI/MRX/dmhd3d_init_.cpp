#include "dmhd3d_class.hpp"

void DMHD3D::init_()
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

  double dvy[nz][nx];
  unsigned seed;
  double stim;			// Time at node 0, used for seed of rand_noise
  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //Random seed depending on x and z, used for vy perturbation
  seed=(unsigned)(stim*(1+mpi_ranx+(mpi_numx*mpi_ranz)));
  
  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      dvy[k][i]=0.0;
#if (RANDOM)
      dvy[k][i]+=rand_noise(dvparay,seed); // Multiple mode perturbation
#else
      dvy[k][i]+=dv*sin(2*M_PI*x[i]/wlen); // Single mode perturbation
#endif    
    }
  }

  //Reset random seed depending on z, used for vx perturbation
  seed=(unsigned)(stim*(1+mpi_ranz));
  srandom(seed);
  
  for (k=0;k<nz;k++){
    double vampx=vamp+rand_noise(dvparax,seed); // Perturbation to vx is independent of x and y    
    for (j=0;j<ny;j++){
      double angle=0.5*((angle_u+angle_l)+(angle_u-angle_l)*tanh((y[j]-s0)/lambda));
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;

	ro[ss]=0.5*((ro_u+ro_l)+(ro_u-ro_l)*tanh((y[j]-s0)/lambda));
	vx[ss]=+vampx*tanh((y[j]-s0)/lambda);
	vy[ss]=dvy[k][i]*exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda));
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
}
