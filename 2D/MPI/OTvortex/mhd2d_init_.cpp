#include "mhd2d_class.hpp"

void MHD2D::init_()
{
  // Orszag-Tang vortex
  int i,j;
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      ro[ss]=gam*gam;
      vx[ss]=-sin(y[j]);
      vy[ss]=+sin(x[i]);
      vz[ss]=0.0;
      bx[ss]=-sin(y[j]);
      by[ss]=+sin(2*x[i]);
      bz[ss]=0.0;
      pr[ss]=gam;
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

