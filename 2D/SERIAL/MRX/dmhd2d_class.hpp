#ifndef _DCLASS_MHD2D_
#define _DCLASS_MHD2D_

#include "mhd2d_class.hpp"

// Dissipative MHD2D

class DMHD2D: public MHD2D{

public:
  double setdc();		// Set dissipation coefficients and return their maximum
  void init_();			// Set initial condition
  void exec_(int);		// Run simulation
  DMHD2D();			// Constructor
  virtual ~DMHD2D();		// Destructor

protected:
  double nu0;			// Kinematic viscosity coefficient
  double eta0;			// Resistivity coefficient
  void dsptv(double);		// Dissipation solver
};

#endif
