#ifndef _CLASS_MHD_
#define _CLASS_MHD_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "mhd_common.h"

class MHD{

public:
  static const int nm=8;	// Number of MHD variables
  const double dtor=M_PI/180.;
  
  void setgam(double gam_value)
  {
    // Set gamma value
    gam=gam_value;
  }
  double getgam()
  {
    // Get gamma value
    return gam;
  }
  double v_snd(int i)
  {
    // Get sound velocity
    return sqrt(gam*pr[i]/ro[i]);
  }
  double v_alf(int i)
  {
    // Get Alfven velocity
    return sqrt((cx[i]*cx[i]+cy[i]*cy[i]+cz[i]*cz[i])/ro[i]);
  }
  double vfast(int i)
  {
    // Get fast magnetosonic velocity
    return sqrt((gam*pr[i]+(cx[i]*cx[i]+cy[i]*cy[i]+cz[i]*cz[i]))/ro[i]);
  }
  
protected:
  double gam=5.0/3.0;		// Specific heat ratio
  double *x,*y,*z;
  double *ro,*mx,*my,*mz,*en;
  double *bx,*by,*bz;		// B @ cell edge
  double *vx,*vy,*vz,*pr;
  double *cx,*cy,*cz;		// B @ cell center
  double *nu,*eta;		// Kinematic viscosity and resistivity
  double *phi_g;		// Gravitational potential
  double *val[nm];		// Pointer for MHD variables
  int dnxs[nm]={0},dnys[nm]={0},dnzs[nm]={0}; // Boundary flag
  int stxs[nm]={0},stys[nm]={0},stzs[nm]={0}; // Staggered grid flag

  void vlcty(int i)
  {
    // Momentum => velocity
    double iro=1.0/ro[i];
    vx[i]=mx[i]*iro;
    vy[i]=my[i]*iro;
    vz[i]=mz[i]*iro;
  }
  void mmntm(int i)
  {
    // Velocity => momentum
    mx[i]=ro[i]*vx[i];
    my[i]=ro[i]*vy[i];
    mz[i]=ro[i]*vz[i];
  }
  void cnsvt(int i)
  {
    // Primitive => conservative
    mmntm(i);
    en[i]=pr[i]/(gam-1)+0.5*(ro[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])+(cx[i]*cx[i]+cy[i]*cy[i]+cz[i]*cz[i]));
  }
  void prmtv(int i)
  {
    // Conservative => primitive
    vlcty(i);
    pr[i]=(gam-1)*(en[i]-0.5*(ro[i]*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])+(cx[i]*cx[i]+cy[i]*cy[i]+cz[i]*cz[i])));
  }
  double bcell(const double* b_ct, int offset, double (*fcen)(const double *f))
  {
    // B @ cell center from cell edge
    double val[4]={b_ct[-offset],b_ct[0],b_ct[+offset],b_ct[+2*offset]};
    return fcen(&val[1]);
  }

  // Routines for MHD solver from mhd_common.h
  /* Riemann solvers */
  void (*riemann[4])(double, double, double, double, double, double, double,
		     double, double, double, double, double, double, double,
		     double, double, const double*,
		     double*, double*, double*, double*, double*, double*, double*)={calc_flux_roe,calc_flux_hlld,calc_flux_lhlld,calc_flux_mlau};
  /* Linear interpolation function to face LR states*/
  void (*l_interp[4])(const double *f, double *fl, double *fr)={cal_flr_1st,cal_flr_2nd,cal_flr_3rd,cal_flr_4th};
  /* Interpolation function to face LR states */
  void (*interpol[4])(const double *f, double *fl, double *fr)={cal_flr_1st,muscl_mm_cal_flr,wcns3_cal_flr,wcns4_cal_flr};
  /* Interpolation from face to center */
  double (*fcen[4])(const double *f)={cal_fcen_2nd,cal_fcen_2nd,cal_fcen_4th,cal_fcen_4th};
  /* Spatial difference */
  double (*df1[4])(const double *f)={cal_df_2nd,cal_df_2nd,cal_df_4th,cal_df_4th};
};

#endif
