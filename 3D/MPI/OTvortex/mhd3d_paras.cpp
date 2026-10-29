#include "mhd3d_class.hpp"

void MHD3D::paras()
{
  // Simulation parameters
  const double xmin_value=0.0;
  const double xmax_value=2.0*M_PI;
  const double ymin_value=0.0;
  const double ymax_value=2.0*M_PI;
  const double zmin_value=0.0;
  const double zmax_value=2.0*M_PI;
  // Global domain bounds, excluding ghost cells; shifts are in cell widths.
  setup_grid(xmin_value,xmax_value,ymin_value,ymax_value,zmin_value,zmax_value,0.0,0.0,0.0);
  fildir="./dat/";
  
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
  
}
