#include "mhd2d_class.hpp"

void MHD2D::paras()
{
  // Simulation parameters
  setgam(2.0);
  const double xmin_value=0.0;
  const double xmax_value=20.0;
  const double ymin_value=0.0;
  const double ymax_value=ymin_value+(xmax_value-xmin_value)*(double)YMESH/XMESH;
  // Physical domain bounds (xmin, xmax, ymin, ymax), excluding ghost cells.
  // Coordinates use the default half-cell shifts (0.5, 0.5).
  setup_grid(xmin_value,xmax_value,ymin_value,ymax_value);
  
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

}
