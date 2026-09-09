#ifndef QASMHD_RMI2D_CLASS_HPP
#define QASMHD_RMI2D_CLASS_HPP

#include "mhd2d_class.hpp"
#include <vector>

class RMI2D : public MHD2D{

public:
  void init_();			// Set initial condition and inflow parameters
  void exec_(int);		// Run with one inflow realization per time step

protected:
  struct InflowParameters{
    double rho,rho_noise;
    double vx,vy,vz;
    double pressure;
    double bx,by,bz;
  };
  InflowParameters inflow={};
  std::vector<double> inflow_ro; // Upper ghost-cell densities, held through all RK stages

  void bound(double *values[], int nvars,
	     const int stxs[], const int dnxs[], const int stys[], const int dnys[]) override;
  // No default argument: the six-argument override above selects refresh_inflow=false.
  void bound(double *values[], int nvars,
	     const int stxs[], const int dnxs[], const int stys[], const int dnys[],
	     bool refresh_inflow);
  void apply_inflow(bool);
  void apply_inflow_center(); // Only modifies the auxiliary cell-center magnetic field
};

#endif
