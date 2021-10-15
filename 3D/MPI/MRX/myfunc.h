#ifndef _MYFUNC_H_
#define _MYFUNC_H_

#include "common_func.h"
#include "mhd_func.h"
#include "mhd_fd3d.h"
#include "boundary.h"

/* Define routines used in main program */

inline double harris_field(double x, const double *params)
{
  /* Harris magnetic field */
  double x0=params[0],width=params[1]+1e-15;
  return( tanh((x-x0)/width) );
}
inline double harris_density(double x, const double *params)
{
  /* Harris distribution of the density */
  double x0=params[0],width=params[1]+1e-15;
  double cosh1=cosh((x-x0)/width);
  return( 1.0/(cosh1*cosh1) );
}
inline double rand_noise(const double *params, unsigned seed)
{
  /* Return uniform random distribution (params[0] +- paramas[1]) */
  static int r_flag=0;
  if (r_flag == 0){
    srandom(seed);
    r_flag=1;
  }
  return(params[0]+params[1]*((double)random()/RAND_MAX-0.5)*2.0);
}

#endif
