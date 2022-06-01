#ifndef _MYFUNC_H_
#define _MYFUNC_H_

#include "common_func.h"
#include "mhd_func.h"
#include "mhd_fd1d.h"
#include "emhd_func1d.h"
#include "boundary.h"

/* Define routines used in main program */

inline double sqrwave(double x, double x0, double x1, double dx)
{
  /* Return square wave distribution */
  return(0.5*(tanh((x-x0)/dx)-tanh((x-x1)/dx)));
}

#endif
