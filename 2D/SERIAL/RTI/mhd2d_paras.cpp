#include "mhd2d_class.hpp"

void MHD2D::paras()
{
  // Simulation parameters
  xmin=0.0;
  xmax=20.0;
  ymin=-0.5*(xmax-xmin)*(double)YMESH/XMESH;
  ymax=+0.5*(xmax-xmin)*(double)YMESH/XMESH;
  dx=(xmax-xmin)/XMESH;
  dy=(ymax-ymin)/YMESH;
  dr=min(dx,dy);
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

  dnys[0]=+2;			// ro
  dnys[1]=+2;			// mx
  dnys[2]=+2;			// my
  dnys[3]=+2;			// mz
  dnys[4]=+2;			// bx
  dnys[5]=+2;			// by
  dnys[6]=+2;			// bz
  dnys[7]=+2;			// en
  
}
