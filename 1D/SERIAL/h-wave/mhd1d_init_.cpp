#include "mhd1d_class.hpp"

void MHD1D::init_()
{
  
  // Initial condition parameters
  const double ro0=1.0;
  const double bx0=1.0;
  const double beta=10.0;
  const double pr0=0.5*beta*bx0*bx0;
  const double cf0=sqrt((bx0*bx0+gam*pr0)/ro0);
  const double vx0=0.0*cf0;
  
  unsigned seed=10;
  // seed=(unsigned)time(NULL);
  double para1[]={0,0.01},para2[]={1,0.01};

  for (int i=0;i<nx;i++){

    ro[i]=ro0;
    vx[i]=vx0;
    vy[i]=0.0;
    vz[i]=0.0;
    bx[i]=bx0;
    by[i]=bx0*rand_noise(para1,seed);
    bz[i]=bx0*rand_noise(para1,seed);
    pr[i]=pr0*rand_noise(para2,seed);
    
    // In 1D, cell center B is identical to cell edge B.
    // Thus cx,cy,cz are not explicitly initialized here.
    // They should be defined in constructer (see mhd1d_class.cpp).
    cnsvt(i);
  }

  // Boundary condition
  bound(val,nm,dnxs);
}

