#include "mhd3d_class.hpp"

void MHD3D::paras()
{
  // Simulation parameters
  const double xmin_value=0.0;
  const double xmax_value=20.0;
  const double ymin_value=0.0;
  const double ymax_value=20.0;
  const double zmin_value=-0.5*(xmax_value-xmin_value)*(double)ZMESH/XMESH;
  const double zmax_value=+0.5*(xmax_value-xmin_value)*(double)ZMESH/XMESH;
  // Global domain bounds, excluding ghost cells; shifts are in cell widths.
  setup_grid(xmin_value,xmax_value,ymin_value,ymax_value,zmin_value,zmax_value);
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

  dnzs[0]=+2;			// ro
  dnzs[1]=+2;			// mx
  dnzs[2]=+2;			// my
  dnzs[3]=+2;			// mz
  dnzs[4]=+2;			// bx
  dnzs[5]=+2;			// by
  dnzs[6]=+2;			// bz
  dnzs[7]=+2;			// en
  
}
