#include "hmhd1d_class.hpp"

void HMHD1D::paras()
{
  // Simulation parameters
  xmin=+0.0;
  xmax=+1.0;
  sprintf(fildir,"./dat/");
  
  // Boundary condition flag for ro,mx,my,mz,bx,by,bz,en
  // 0: periodic, +1: Neumann, -1: Dirichlet
  dnxs[0]=+0;			// ro
  dnxs[1]=+0;			// mx
  dnxs[2]=+0;			// my
  dnxs[3]=+0;			// mz
  dnxs[4]=+0;			// bx
  dnxs[5]=+0;			// by
  dnxs[6]=+0;			// bz
  dnxs[7]=+0;			// en

  dx=(xmax-xmin)/XMESH;
  dt=cfl*dx;			// dt will be re-calculated later
  for (int i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx+xmin;
  }

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
