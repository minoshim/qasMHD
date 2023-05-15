#include "mhd1d_class.hpp"

void MHD1D::paras()
{
  // Simulation parameters
  xmin=-0.5;
  xmax=+0.5;
  dx=(xmax-xmin)/XMESH;
  dt=cfl*dx;
  // dt will be re-calculated later
  sprintf(fildir,"./dat/");
  
  // Boundary condition flag for ro,mx,my,mz,bx,by,bz,en
  // 0: periodic, +1: Neumann, -1: Dirichlet, +2: Open, -2: Zero fixed
  dnxs[0]=+2;			// ro
  dnxs[1]=+2;			// mx
  dnxs[2]=+2;			// my
  dnxs[3]=+2;			// mz
  dnxs[4]=+2;			// bx
  dnxs[5]=+2;			// by
  dnxs[6]=+2;			// bz
  dnxs[7]=+2;			// en

}
