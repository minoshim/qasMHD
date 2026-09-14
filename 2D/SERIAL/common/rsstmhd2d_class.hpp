#ifndef RSSTMHD2D_CLASS_
#define RSSTMHD2D_CLASS_

#include "mhd2d_class.hpp"

#ifndef RSST_FORM
#define RSST_FORM 0 // Preserve PVS for cases without an explicit selection.
#endif
static_assert(RSST_FORM == 0 || RSST_FORM == 1,
              "RSST_FORM must be 0 (PVS) or 1 (PMS)");

// MHD2D with Reduced Speed of Sound Technique (Iijima+19)

class RSSTMHD2D: public MHD2D{

public:
  void paras();		// Set parameters
  void setdt(int) override;	// Set time step
  RSSTMHD2D();		// Constructor
  virtual ~RSSTMHD2D();	// Destructor

protected:
  double* ixi;			// RSST factor 1/\xi
  double cs_bnd;		// Sound speed threshold and lower bound for reduction
  double ma_bnd;		// Maximum Mach number allowed by sound speed reduction
  void check_rsst_parameters() const; // Require finite, positive RSST parameters
  void rsst_(int);		// Set RSST factor (1/\xi)
  void ideal(double) override;
};

#endif
