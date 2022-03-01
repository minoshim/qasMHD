#ifndef _MYFUNC_H_
#define _MYFUNC_H_

#include "common_func.h"
#include "mhd_func.h"
#include "mhd_fd3d.h"
#include "boundary.h"

/* Define routines used in main program */

double sqrwave2(double val_u, double val_l, double dx_u, double dx_l)
// Calculate double square wave profile
/* Return val_u (dx_u*dx_l>0) or val_l (dx_u*dx_l<0) */
{
  return 0.5*(2.0+tanh(dx_u)-tanh(dx_l))*(val_u-val_l)+val_l;
}

double g_potential(double z, double g0, double lg, int deriv)
// Calculate gravitational potential
// phi_g = g0*lg*log(cosh(z/lg))
// d(phi_g)/dz = g0*tanh(z/lg)
// Return z-derivative if deriv==1
{
  return (deriv == 1)?(g0*tanh(z/lg)):(g0*lg*log(cosh(z/lg)));
}

#endif