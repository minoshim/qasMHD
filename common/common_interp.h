#ifndef _COMMON_INTERP_H_
#define _COMMON_INTERP_H_

/* Central interpolation */
/* f = address @ i-1/2, Return = value @ i */
inline double cal_fcen_2nd(const double *f) /* 2nd order */
{
  return(0.5*(f[0]+f[1]));
}
inline double cal_fcen_4th(const double *f) /* 4th order */
{
  return(0.0625*(9.0*(f[0]+f[1])-(f[-1]+f[2])));
}

/* Central difference */
/* f = address @ i-1/2, Return = 1st derivative @ i */
inline double cal_df_2nd(const double *f) /* 2nd order */
{
  return(f[1]-f[0]);
}
inline double cal_df_4th(const double *f) /* 4th order */
{
  return((27.0*(f[1]-f[0])-(f[2]-f[-1]))/24.0);
}

/* 2nd central difference */
/* f = address @ i, Return = 2nd derivative @ i */
inline double cal_d2f_2nd(const double *f) /* 2nd order */
{
  return ( (f[1]+f[-1])-2.0*f[0] );
}
inline double cal_d2f_4th(const double *f) /* 4th order */
{
  return ( (-(f[2]+f[-2])+16.0*(f[1]+f[-1])-30.0*f[0])/12.0 );
}

/* Rotation of cell-center variables with 2nd order accuracy */
/* cx,cy @ i,j. Return = (dcy/dx-dcx/dy) @ i,j  */
inline double rotc(const double *cx, const double *cy, double idx, double idy, int xoffset, int yoffset)
{
  // rotC, where C is cell-center quantity
  return 0.5*((cy[+xoffset]-cy[-xoffset])*idx-(cx[+yoffset]-cx[-yoffset])*idy);
}

/* Linear interpolation */
/* f = address @ i, *fl = left state @ i+1/2, *fr = Right state @i-1/2 */
void cal_flr_1st(const double *f, double *fl, double *fr); /* 1st order */
void cal_flr_2nd(const double *f, double *fl, double *fr); /* 2nd order  */
void cal_flr_3rd(const double *f, double *fl, double *fr); /* 3rd order */
void cal_flr_4th(const double *f, double *fl, double *fr); /* 4th order */

/* MUSCL routines */
/* f = address @ i, *fl = left state @ i+1/2, *fr = Right state @i-1/2 */
void muscl_mm_cal_flr(const double *f, double *fl, double *fr); /* MinMod */
void muscl_mc_cal_flr(const double *f, double *fl, double *fr); /* MC */
void muscl_kr_cal_flr(const double *f, double *fl, double *fr); /* Koren */

/* WCNS routines */
/* f = address @ i, *fl = left state @ i+1/2, *fr = Right state @i-1/2 */
void wcns3_cal_flr(const double *f, double *fl, double *fr); /* 3rd order WCNS (Minoshima+19) */
void wcns4_cal_flr(const double *f, double *fl, double *fr); /* 4th order WCNS (Minoshima+19) */

#endif
