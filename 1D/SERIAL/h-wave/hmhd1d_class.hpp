#ifndef _CLASS_HMHD1D_
#define _CLASS_HMHD1D_

#include "mhd1d_class.hpp"

// Hall MHD1D

class HMHD1D: public MHD1D{

public:
  double haldt();
  void paras();			// Set parameters
  void exec_(int);		// Run simulation
  HMHD1D();			// Constructor
  virtual ~HMHD1D();		// Destructor

protected:
  double *hx,*hy,*hz;	       // Hall velocity
  double di;		       // Ion inertia length
  double de;		       // "Artificial" electron inertia length
  double idx;
  double vphix;			// Maximum whistler phase velocity
  void hall_(double);	       // Hall term solver
  
  double rotc(const double *cx, const double *cy, double idx, double idy, int xoffset, int yoffset)
  {
    // rotC, where C is cell-center quantity
    return 0.5*((cy[+xoffset]-cy[-xoffset])*idx-(cx[+yoffset]-cx[-yoffset])*idy);
  }
  void hallv(int i)
  {
    // Get Hall velocity -(di/n)*rotB
    double fac=di/ro[i];
    hx[i]=0.0;
    hy[i]=-fac*rotc(&cz[i],&cx[i],0.0,idx,0,1);
    hz[i]=-fac*rotc(&cx[i],&cy[i],idx,0.0,1,0);
  }
};

#endif
