#include "mhd2d_class.hpp"

void MHD2D::init_()
{
  // Magnetic loop advection in high beta plasma
  int i,j;
  // Initial condition parameters
  const double ro0=1.0;		// Ambient density
  const double pr0=1e0;		// Ambient pressure
  const double v0=sqrt(5.0);	// Ambient velocity
  // const double theta=0.463647609000806; // Direction of velocity
  const double theta=0.01;		// Lee13
  const double vx0=v0*cos(theta);
  const double vy0=v0*sin(theta);
  const double vz0=0.0;
  const double a0=1e-3;		// Amplitude of magnetic loop
  const double rad=0.3;		// Radius of magnetic loop

  double *azp,*azc;
  azp=new double[nd];
  azc=new double[nd];
  // Vector potential
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double r,x0,y0;
      x0=x[i]-0.5*dx;
      y0=y[j]-0.5*dy;
      r=sqrt(x0*x0+y0*y0);
      azp[ss]=(r <= rad)?a0*(rad-r):0;
      x0=x[i];
      y0=y[j];
      r=sqrt(x0*x0+y0*y0);
      azc[ss]=(r <= rad)?a0*(rad-r):0;
    }
  }
  // MHD variables
  double idx=1.0/dx,idy=1.0/dy;
  for (j=1;j<ny-2;j++){
    for (i=1;i<nx-2;i++){
      int ss=nx*j+i;
      ro[ss]=ro0;
      vx[ss]=vx0;
      vy[ss]=vy0;
      vz[ss]=vz0;
      pr[ss]=pr0;

      double val1[4]={azp[nx*(j-1)+i],azp[nx*(j+0)+i],azp[nx*(j+1)+i],azp[nx*(j+2)+i]};
      double val2[4]={azp[nx*j+(i-1)],azp[nx*j+(i+0)],azp[nx*j+(i+1)],azp[nx*j+(i+2)]};
      bx[ss]=+(df1[ODR-1])(&val1[1])*idy;
      by[ss]=-(df1[ODR-1])(&val2[1])*idx;
      bz[ss]=0.0;
      cx[ss]=+(azc[nx*(j+1)+i]-azc[nx*(j-1)+i])*idy*0.5;
      cy[ss]=-(azc[nx*j+(i+1)]-azc[nx*j+(i-1)])*idx*0.5;

      // In 2D, cell center Bz is identical to cell edge Bz.
      // Thus cz is not explicitly initialized here.
      // They should be defined in constructer (see mhd2d_class.cpp).
      cnsvt(ss);
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys);

  delete[] azp;
  delete[] azc;
}

