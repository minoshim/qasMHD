#ifndef _COMMON_INTERP_H_
#define _COMMON_INTERP_H_

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

/* Central interpolation */
/* f = address @ i-1/2, Return = value @ i */
double cal_fcen_2nd(const double *f); /* 2nd order */
double cal_fcen_4th(const double *f); /* 4th order */

/* Central difference */
/* f = address @ i-1/2, Return = 1st derivative @ i */
double cal_df_2nd(const double *f); /* 2nd order */
double cal_df_4th(const double *f); /* 4th order */

/* 2nd central difference */
/* f = address @ i, Return = 2nd derivative @ i */
double cal_d2f_2nd(const double *f); /* 2nd order */
double cal_d2f_4th(const double *f); /* 4th order */

#endif
