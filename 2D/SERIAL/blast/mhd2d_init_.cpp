#include "mhd2d_class.hpp"

void MHD2D::init_()
{
  // Blast wave
  int i,j;
  // Initial condition parameters
  const double r_0=0.125;       // Radius of imposed high-P cylinder
  const double ro0=1e0;         // Ambient density
  const double ro1=1e0;         // Density in cylinder
  const double pr0=1e0;         // Ambient pressure
  const double pr1=1e2;         // Pressure in cylinder
  const double b_0=10.0;        // Ambient |B|
  const double ban=30.0;        // B angle relative to x axis

  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double rr=sqrt(x[i]*x[i]+y[j]*y[j]);
      ro[ss]=(rr <= r_0)?ro1:ro0;
      vx[ss]=0.0;
      vy[ss]=0.0;
      vz[ss]=0.0;
      bx[ss]=b_0*cos(ban*dtor);
      by[ss]=b_0*sin(ban*dtor);
      bz[ss]=0.0;
      pr[ss]=(rr <= r_0)?pr1:pr0;
      cx[ss]=bx[ss];
      cy[ss]=by[ss];

      // In 2D, cell center Bz is identical to cell edge Bz.
      // Thus cz is not explicitly initialized here.
      // They should be defined in constructer (see mhd2d_class.cpp).
      cnsvt(ss);
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys);
}

