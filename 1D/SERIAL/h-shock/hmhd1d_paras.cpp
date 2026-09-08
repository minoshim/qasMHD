#include "hmhd1d_class.hpp"

void HMHD1D::paras()
{
  // Simulation parameters
  setup_grid(-0.5,+0.5);
  
  // Boundary condition flag for ro,mx,my,mz,bx,by,bz,en
  // 0: periodic, +1: Neumann, -1: Dirichlet
  dnxs[0]=+1;			// ro
  dnxs[1]=+1;			// mx
  dnxs[2]=+1;			// my
  dnxs[3]=+1;			// mz
  dnxs[4]=+1;			// bx
  dnxs[5]=+1;			// by
  dnxs[6]=+1;			// bz
  dnxs[7]=+1;			// en

  // Hall parameters
  di=D_ION;
  de=dx*D_ELE;
  idx=1.0/dx;
  double kpeak=M_PI*idx;
  if (de != 0) kpeak=min(kpeak,1.0/de);
  vphix=di*kpeak;
  
  // Caution:
  // The phase velocity of whistler waves peaks @ kpeak=1.0/de when electron inertia effect is included.
  // Use of vphix=di/de reduces computational cost for Hall solver.
  // However, it also reduces the amount of numerical dissipation for high freq. waves.
}
