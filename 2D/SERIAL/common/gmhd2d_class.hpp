#ifndef QASMHD_GMHD2D_CLASS_HPP
#define QASMHD_GMHD2D_CLASS_HPP

#include "mhd2d_class.hpp"

// MHD with a prescribed, time-independent gravitational potential.
// The conserved total energy includes ro*phi_g; no self-gravity solver is provided.
class GMHD2D : public MHD2D{

public:
  GMHD2D();
  ~GMHD2D() override;
  void init_();                 // Problem-specific initial state and potential
  void setdt(int) override;

protected:
  // Non-virtual: GMHD2D routines use gravity-aware primitive conversion.
  void prmtv(int i)
  {
    prmtv_gravity(i);
  }

  // Problem-specific boundary conditions, implemented alongside init_().
  void bound(double *values[], int nvars,
             const int stxs[], const int dnxs[], const int stys[], const int dnys[]) override;
  void ideal(double) override;
  void dout_(int) override;
};

#endif
