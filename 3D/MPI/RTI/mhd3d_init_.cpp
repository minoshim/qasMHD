#include "mhd3d_class.hpp"

inline double sqrwave2(double val_u, double val_l, double dx_u, double dx_l);
inline double g_potential(double z, double g0, double lg, int deriv);
double cal_pressure(double y0, double y1, double pr0, int n,
		    double val_u, double val_l, double s0, double lambda, double g0, double lg);

void MHD3D::init_()
{
  // RT instability
  int i,j,k;
  // Initial condition parameters
  const double beta=1e2;	// Ambient plasma beta
  const double angle_u=90;	// B field angle in upper domain. 90deg: B=By, 0deg: B=Bx
  const double angle_l=angle_u;	// B field angle in lower domain.
  const double vamp=0.0;	// Shear velocity amplitude
  const int nmode=1;		// Number of mode for perturbation
  const double wlen=getlx()/nmode;
  const double lambda=1.0;	// Shear layer width
  const double s0=0.35*getlz(); // Shear position (@ +s0 and -s0)
  const double ro_u=1.0;	// Density in upper domain
  const double ro_l=0.2;	// Density in lower domain
  const double pr0=50.0;	// Base pressure
  const double b0=sqrt(2.0*pr0/beta); // B field strength
  const double prmin=1e-4;	// Minimum pressure threshold
  const double dv=0.01;		// Perturbation amplitude
  // Parameters associated with gravity
  const double g0=1.0;		// Gravitational acceleration
  const double lg=8*dy;		// Width of boundary layer around y=0 (for gravity profile)

  double dvz[ny][nx];
  unsigned seed;
  double stim;			// Time at node 0, used for seed of rand_noise
  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //Random seed depending on x and y, used for vz perturbation
  seed=(unsigned)(stim*(1+(mpi_rank%m_xy)));
  
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      dvz[j][i]=0.0;
#if (RANDOM)
      double dvpara[2]={0,dv};
      dvz[j][i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#else
      // dvz[j][i]+=dv*cos(2*M_PI*x[i]/wlen); // Single mode perturbation
      // dvz[j][i]+=dv*cos(2*M_PI*y[j]/wlen); // Single mode perturbation
      dvz[j][i]+=dv*cos(2*M_PI*x[i]/wlen)*cos(2*M_PI*y[j]/wlen); // Single mode perturbation
#endif
    }
  }

  for (k=0;k<nz;k++){
    double angle=sqrwave2(angle_u,angle_l,(z[k]-s0)/lambda,(z[k]+s0)/lambda);
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	
	ro[ss]=sqrwave2(ro_u,ro_l,(z[k]-s0)/lambda,(z[k]+s0)/lambda);
	vx[ss]=sqrwave2(vamp,-vamp,(z[k]-s0)/lambda,(z[k]+s0)/lambda);
	vy[ss]=0.0;	
	vz[ss]=dvz[j][i]*(exp(-((z[k]-s0)*(z[k]-s0))/(4*lambda*lambda))-exp(-((z[k]+s0)*(z[k]+s0))/(4*lambda*lambda)));
	
	// Gravitational potential, phi_g = g0*lg*log(cosh(z/lg))
	// d(phi_g)/dz = g0*tanh(z/lg)
	phi_g[ss]=g_potential(z[k],g0,lg,0);
      
	// Pressure is obtained by integrating the equilibrium
	pr[ss]=cal_pressure(0,z[k],pr0,1+(int)(fabs(z[k])/dz+0.5),
			    ro_u,ro_l,s0,lambda,g0,lg);
	pr[ss]=max(pr[ss],prmin);

	if (angle == 90){
	  bx[ss]=0.0;
	  by[ss]=b0;
	  bz[ss]=0.0;
	} else{
	  bx[ss]=b0*cos(angle*dtor);
	  by[ss]=b0*sin(angle*dtor);
	  bz[ss]=0.0;
	}
	cx[ss]=bx[ss];
	cy[ss]=by[ss];
	cz[ss]=bz[ss];
	
	cnsvt(ss);
	en[ss]+=ro[ss]*phi_g[ss];	// Add gravitational potential
      }
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
}

inline double sqrwave2(double val_u, double val_l, double dx_u, double dx_l)
// Calculate double square wave profile
/* Return val_u (dx_u*dx_l>0) or val_l (dx_u*dx_l<0) */
{
  return 0.5*(2.0+tanh(dx_u)-tanh(dx_l))*(val_u-val_l)+val_l;
}

inline double g_potential(double z, double g0, double lg, int deriv)
// Calculate gravitational potential
// phi_g = g0*lg*log(cosh(z/lg))
// d(phi_g)/dz = g0*tanh(z/lg)
// Return z-derivative if deriv==1
{
  return (deriv == 1)?(g0*tanh(z/lg)):(g0*lg*log(cosh(z/lg)));
}

double cal_pressure(double y0, double y1, double pr0, int n,
		    double val_u, double val_l, double s0, double lambda, double g0, double lg)
// Calculate pressure profile from hydrostatic equilibrium dP/dz = -rho dPhi_g/dz
{
  int j;
  double dyc=(y1-y0)/n;
  double ans=pr0;
  for (j=0;j<n;j++){
    double yc,roc,dphi_gc;
    yc=y0+(j+0.5)*dyc;
    roc=sqrwave2(val_u,val_l,(yc-s0)/lambda,(yc+s0)/lambda);
    dphi_gc=g_potential(yc,g0,lg,1);
    ans+=-roc*dphi_gc*dyc;
  }
  return ans;
}

