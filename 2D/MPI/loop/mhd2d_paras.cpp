#include "mhd2d_class.hpp"

void MHD2D::paras()
{
  // Simulation parameters
  xmin=-1.0;
  xmax=+1.0;
  ymin=-0.5;
  ymax=+0.5;
  sprintf(fildir,"./dat/");
  
  // Boundary condition flag for ro,mx,my,mz,bx,by,bz,en
  // 0: periodic, +1: Neumann, -1: Dirichlet, +2: Open, -2: Zero fixed
  dnxs[0]=+0;			// ro
  dnxs[1]=+0;			// mx
  dnxs[2]=+0;			// my
  dnxs[3]=+0;			// mz
  dnxs[4]=+0;			// bx
  dnxs[5]=+0;			// by
  dnxs[6]=+0;			// bz
  dnxs[7]=+0;			// en

  dnys[0]=+0;			// ro
  dnys[1]=+0;			// mx
  dnys[2]=+0;			// my
  dnys[3]=+0;			// mz
  dnys[4]=+0;			// bx
  dnys[5]=+0;			// by
  dnys[6]=+0;			// bz
  dnys[7]=+0;			// en
  
  dx=(xmax-xmin)/XMESH;
  dy=(ymax-ymin)/YMESH;
  dr=min(dx,dy);
  dt=cfl*dr;			// dt will be re-calculated later
  int isum=0,jsum=0,m;
  for (m=0;m<mpi_ranx;m++){
    isum+=(XMESH+m)/mpi_numx;
  }
  for (m=0;m<mpi_rany;m++){
    jsum+=(YMESH+m)/mpi_numy;
  }
  for (int i=0;i<nx;i++) x[i]=(i-xoff+isum+0.5)*dx+xmin;
  for (int j=0;j<ny;j++) y[j]=(j-yoff+jsum+0.5)*dy+ymin;

}
