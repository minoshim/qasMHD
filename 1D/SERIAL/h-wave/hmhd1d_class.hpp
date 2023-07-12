#ifndef _CLASS_HMHD1D_
#define _CLASS_HMHD1D_

#include "mhd1d_class.hpp"

// Hall MHD1D

class HMHD1D: public MHD1D{

public:
  void paras();			// Set parameters
  double haldt();		// Return time step for Hall term
  void exec_(int);		// Run simulation
  HMHD1D();			// Constructor
  virtual ~HMHD1D();		// Destructor

protected:
  double *hx,*hy,*hz;		// Hall velocity
  double di,de; 		// Ion and "artificial" electron inertia lengths
  double idx;
  double vphix;			// Maximum whistler phase velocity
  void hall_(double);		// Hall term solver
  
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
