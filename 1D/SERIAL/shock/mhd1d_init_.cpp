#include "mhd1d_class.hpp"

// Select initial condition. See mhd1d_init_.cpp
#define NUM (1)                 
// NUM = 1: Dai & Woodward 1994 (Miyoshi & Kusano 2005, Fig. 5)
// NUM = 2: Brio & Wu 1988 (Miyoshi & Kusano 2005, Fig. 8)
// NUM = 3: Slow switch-off shock (Miyoshi & Kusano 2005, Fig. 9)
// NUM = 4: Slow switch-off rarefaction (Miyoshi & Kusano 2005, Fig. 10)
// NUM = 5: Super-fast expansion (Miyoshi & Kusano 2005, Fig. 11)

void MHD1D::init_()
{
  
#if (NUM == 2)
  setgam(2.0);
#endif

  for (int i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx+xmin;
    char flag=(x[i] <= 0);

#if (NUM == 2)
    // BrioWu
    ro[i]=(flag)?(1.00):(0.125);
    vx[i]=0;
    vy[i]=0;
    vz[i]=0;
    bx[i]=0.75;
    by[i]=(flag)?(1):(-1);
    bz[i]=0;
    pr[i]=(flag)?(1.0):(0.1);
#elif (NUM ==3)
  // Slow switch-off shock
    ro[i]=(flag)?(1.368):(1.0);
    vx[i]=(flag)?(0.269):(0.0);
    vy[i]=(flag)?(1.0):(0.0);
    vz[i]=(flag)?(0.0):(0.0);
    bx[i]=1.0
    by[i]=(flag)?(0.0):(1.0);
    bz[i]=0.0;
    pr[i]=(flag)?(1.769):(1.0);
#elif (NUM == 4)
  // Slow switch-off rarefaction
    ro[i]=(flag)?(1.0):(0.2);
    vx[i]=(flag)?(0.0):(1.186);
    vy[i]=(flag)?(0.0):(2.967);
    vz[i]=(flag)?(0.0):(0.0);
    bx[i]=1.0
    by[i]=(flag)?(0.0):(1.6405);
    bz[i]=0.0;
    pr[i]=(flag)?(2.0):(0.1368);
#elif (NUM == 5)    
  // Super-fast expansion
    ro[i]=1.0;
    vx[i]=(flag)?(-3.1):(3.1);
    vy[i]=0.0;
    vz[i]=0.0;
    bx[i]=0.0
    by[i]=0.5;
    bz[i]=0.0;
    pr[i]=0.45;
#else  // Default
    // DaiWoodward
    ro[i]=(flag)?(1.08):(1.0);
    vx[i]=(flag)?(1.20):(0.0);
    vy[i]=(flag)?(0.01):(0.0);
    vz[i]=(flag)?(0.50):(0.0);
    bx[i]=2.0/sqrt(4.0*M_PI);
    by[i]=(flag)?(3.6/sqrt(4*M_PI)):(4.0/sqrt(4*M_PI));
    bz[i]=2.0/sqrt(4*M_PI);
    pr[i]=(flag)?(0.95):(1.0);
#endif

    // In 1D, cell center B is identical to cell edge B.
    // Thus cx,cy,cz are not explicitly initialized here.
    // They should be defined in constructer (see mhd1d_class.cpp).
    cnsvt(i);
  }

  // Boundary condition
  bound(val,nm,dnxs);
}

