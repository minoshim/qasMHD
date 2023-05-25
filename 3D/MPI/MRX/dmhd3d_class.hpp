#ifndef _DCLASS_MHD3D_
#define _DCLASS_MHD3D_

#include "mhd3d_class.hpp"

// Dissipative MHD3D

class DMHD3D: public MHD3D{

public:
  double setdc();		// Set dissipation coefficients and return their maximum
  void init_();			// Set initial condition
  void exec_(int);		// Run simulation
  DMHD3D(int*, char***, int);	// Constructor
  virtual ~DMHD3D();		// Destructor

protected:
  double nu0;			// Kinematic viscosity coefficient
  double eta0;			// Resistivity coefficient
  void dsptv(double);		// Dissipation solver
};

#endif
