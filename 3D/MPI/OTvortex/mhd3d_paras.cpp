#include "mhd3d_class.hpp"

void MHD3D::paras()
{
  // Simulation parameters
  xmin=0.0;
  xmax=2.0*M_PI;
  ymin=0.0;
  ymax=2.0*M_PI;
  zmin=0.0;
  zmax=2.0*M_PI;
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

  dnzs[0]=+0;			// ro
  dnzs[1]=+0;			// mx
  dnzs[2]=+0;			// my
  dnzs[3]=+0;			// mz
  dnzs[4]=+0;			// bx
  dnzs[5]=+0;			// by
  dnzs[6]=+0;			// bz
  dnzs[7]=+0;			// en
  
  dx=(xmax-xmin)/XMESH;
  dy=(ymax-ymin)/YMESH;
  dz=(zmax-zmin)/ZMESH;
  dr=min(min(dx,dy),dz);
  dt=cfl*dr;			// dt will be re-calculated later
  int isum=0,jsum=0,ksum=0,m;
  for (m=0;m<mpi_ranx;m++){
    isum+=(XMESH+m)/mpi_numx;
  }
  for (m=0;m<mpi_rany;m++){
    jsum+=(YMESH+m)/mpi_numy;
  }
  for (m=0;m<mpi_ranz;m++){
    ksum+=(ZMESH+m)/mpi_numz;
  }
  for (int i=0;i<nx;i++) x[i]=(i-xoff+isum)*dx+xmin;
  for (int j=0;j<ny;j++) y[j]=(j-yoff+jsum)*dy+ymin;
  for (int k=0;k<nz;k++) z[k]=(k-zoff+ksum)*dz+zmin;

}
