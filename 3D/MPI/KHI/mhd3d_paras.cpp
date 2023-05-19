#include "mhd3d_class.hpp"

void MHD3D::paras()
{
  // Simulation parameters
  setgam(2.0);
  xmin=0.0;
  xmax=20.0;
  ymin=0.0;
  ymax=ymin+(xmax-xmin)*(double)YMESH/XMESH;
  zmin=0.0;
  zmax=zmin+(xmax-xmin)*(double)ZMESH/XMESH;
  dx=(xmax-xmin)/XMESH;
  dy=(ymax-ymin)/YMESH;
  dz=(zmax-zmin)/ZMESH;
  dr=min(min(dx,dy),dz);
  dt=cfl*dr;
  // dt will be re-calculated later
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

  dnys[0]=+1;			// ro
  dnys[1]=+1;			// mx
  dnys[2]=-1;			// my
  dnys[3]=+1;			// mz
  dnys[4]=+1;			// bx
  dnys[5]=-1;			// by
  dnys[6]=+1;			// bz
  dnys[7]=+1;			// en

  dnzs[0]=+0;			// ro
  dnzs[1]=+0;			// mx
  dnzs[2]=+0;			// my
  dnzs[3]=+0;			// mz
  dnzs[4]=+0;			// bx
  dnzs[5]=+0;			// by
  dnzs[6]=+0;			// bz
  dnzs[7]=+0;			// en
  
}
