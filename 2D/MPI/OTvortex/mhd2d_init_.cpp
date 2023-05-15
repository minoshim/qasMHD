#include "mhd2d_class.hpp"

void MHD2D::init_()
{
  int i,j;
  int isum=0,jsum=0,m;
  for (m=0;m<mpi_ranx;m++){
    isum+=(XMESH+m)/mpi_numx;
  }
  for (m=0;m<mpi_rany;m++){
    jsum+=(YMESH+m)/mpi_numy;
  }
  for (i=0;i<nx;i++) x[i]=(i-xoff+isum)*dx+xmin;
  for (j=0;j<ny;j++) y[j]=(j-yoff+jsum)*dy+ymin;
  
  // Orszag-Tang vortex
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

