#ifndef _MYFUNC_H_
#define _MYFUNC_H_

#include "common_func.h"
#include "mhd_func.h"
#include "mhd_fd3d.h"
#include "boundary.h"

/* Define routines used in main program */

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
