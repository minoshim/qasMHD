#include "mhd2d_class.hpp"

double sqrwave2(double val_u, double val_l, double dx_u, double dx_l);
double g_potential(double z, double g0, double lg, int deriv);
double cal_pressure(double y0, double y1, double pr0, int n,
		    double val_u, double val_l, double s0, double lambda, double g0, double lg);

void MHD2D::init_()
{
  int i,j;
  for (i=0;i<nx;i++) x[i]=(i-xoff+0.5)*dx+xmin;
  for (j=0;j<ny;j++) y[j]=(j-yoff+0.5)*dy+ymin;

  // RT instability

  // Initial condition parameters
  const double beta=1e3;	// Ambient plasma beta
  const double angle_u=90;	// B field angle in upper domain. 90deg: B=Bz, 0deg: B=Bx
  const double angle_l=angle_u;	// B field angle in lower domain.
  const double vamp=0.0;	// Shear velocity amplitude
  const int nmode=2;		// Number of mode for perturbation
  const double wlen=getlx()/nmode;
  const double lambda=1.0;	// Shear layer width
  const double s0=0.35*(ymax-ymin); // Shear position (@ +s0 and -s0)
  const double ro_u=1.0;	// Density in upper domain
  const double ro_l=0.2;	// Density in lower domain
  const double b0=1.0;		// B field strength
  const double pr0=0.5*beta*b0*b0; // Base pressure
  const double prmin=1e-4;	// Minimum pressure threshold
  const double dv=0.01;		// Perturbation amplitude
  // Parameters associated with gravity
  const double g0=1.0;		// Gravitational acceleration
  const double lg=8*dy;		// Width of boundary layer around y=0 (for gravity profile)

  double *dvy=new double[nx];
  unsigned seed=(unsigned)time(NULL);
  for (i=0;i<nx;i++){
    dvy[i]=0.0;
#if (RANDOM)
    double dvpara[2]={0,dv};
    dvy[i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#else
    dvy[i]+=dv*cos(2*M_PI*x[i]/wlen); // Single mode perturbation
#endif    
  }
  
  for (j=0;j<ny;j++){
    double angle=sqrwave2(angle_u,angle_l,(y[j]-s0)/lambda,(y[j]+s0)/lambda);
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      
      ro[ss]=sqrwave2(ro_u,ro_l,(y[j]-s0)/lambda,(y[j]+s0)/lambda);
      vx[ss]=sqrwave2(vamp,-vamp,(y[j]-s0)/lambda,(y[j]+s0)/lambda);
      vy[ss]=dvy[i]*(exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda))-exp(-((y[j]+s0)*(y[j]+s0))/(4*lambda*lambda)));
      vz[ss]=0.0;

      // Gravitational potential, phi_g = g0*lg*log(cosh(y/lg))
      // d(phi_g)/dy = g0*tanh(y/lg)
      phi_g[ss]=g_potential(y[j],g0,lg,0);

      // Pressure is obtained by integrating the equilibrium
      pr[ss]=cal_pressure(0,y[j],pr0,1+(int)(fabs(y[j])/dy+0.5),
			  ro_u,ro_l,s0,lambda,g0,lg);
      pr[ss]=max(pr[ss],prmin);
      
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

      // In 2D, cell center Bz is identical to cell edge Bz.
      // Thus cz is not explicitly initialized here.
      // They should be defined in constructer (see mhd2d_class.cpp).
      cnsvt(ss);
      en[ss]+=ro[ss]*phi_g[ss];	// Add gravitational potential
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys);

  delete[] dvy;
}

double sqrwave2(double val_u, double val_l, double dx_u, double dx_l)
// Calculate double square wave profile
/* Return val_u (dx_u*dx_l>0) or val_l (dx_u*dx_l<0) */
{
  return 0.5*(2.0+tanh(dx_u)-tanh(dx_l))*(val_u-val_l)+val_l;
}

double g_potential(double z, double g0, double lg, int deriv)
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

