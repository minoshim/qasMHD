#ifndef _CLASS_HMHD1D_
#define _CLASS_HMHD1D_

#include "mhd1d_class.hpp"

// Hall MHD1D

class HMHD1D: public MHD1D{

public:
  void exec_(int);		// Run simulation
  HMHD1D();			// Constructor
  virtual ~HMHD1D();		// Destructor

protected:
  double *hx,*hy,*hz;	       // Hall velocity
  const double di=dx*D_ION;	       // Ion inertia length
  const double de=dx*D_ELE;	       // "Artificial" electron inertia length
  const double idx=1.0/dx;
  void hall_(double);	       // Hall term solver

  double rotc(const double *cx, const double *cy, double idx, double idy, int xoffset, int yoffset)
  {
    // rotC, where C is cell-center quantity
    return 0.5*((cy[+xoffset]-cy[-xoffset])*idx-(cx[+yoffset]-cx[-yoffset])*idy);
  }
  void hallv(int i)
  {
    // Get Hall velocity
    double fac=di/ro[i];
    hx[i]=0.0;
    hy[i]=-fac*rotc(&cz[i],&cx[i],0.0,idx,0,1);
    hz[i]=-fac*rotc(&cx[i],&cy[i],idx,0.0,1,0);
  }
};

#endif
