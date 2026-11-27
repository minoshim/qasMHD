#include "common_interp.h"
#include <math.h>
#include <cmath>
#include "common_func.h"

/* Linear interpolation */
/* f = address @ i, *fl = left state @ i+1/2, *fr = Right state @i-1/2 */
void cal_flr_1st(const double *f, double *fl, double *fr) /* 1st order */
{
  (*fl)=f[0];
  (*fr)=f[0];
}
void cal_flr_2nd(const double *f, double *fl, double *fr) /* 2nd order  */
{
  /* (*fl)=0.25*(-f[-1]+4.0*f[0]+f[+1]); */
  /* (*fr)=0.25*(-f[+1]+4.0*f[0]+f[-1]); */
  (*fl)=(-f[-1]+5.0*f[0]+2.0*f[+1])/6.0;
  (*fr)=(-f[+1]+5.0*f[0]+2.0*f[-1])/6.0;
}
void cal_flr_3rd(const double *f, double *fl, double *fr) /* 3rd order */
{
  (*fl)=0.125*(-f[-1]+6.0*f[0]+3.0*f[+1]);
  (*fr)=0.125*(-f[+1]+6.0*f[0]+3.0*f[-1]);
}
void cal_flr_4th(const double *f, double *fl, double *fr) /* 4th order */
{
  (*fl)=0.003125*(9.0*f[-2]-56.0*f[-1]+234.0*f[0]+144.0*f[+1]-11.0*f[+2]);
  (*fr)=0.003125*(9.0*f[+2]-56.0*f[+1]+234.0*f[0]+144.0*f[-1]-11.0*f[-2]);
}

/* MUSCL routines */
/* f = address @ i, *fl = left state @ i+1/2, *fr = Right state @i-1/2 */
void muscl_mm_cal_flr(const double *f, double *fl, double *fr) /* MinMod */
{
  double df=0.5*minmod(f[0]-f[-1],f[+1]-f[0]);
  (*fl)=f[0]+df;
  (*fr)=f[0]-df;
}
void muscl_mc_cal_flr(const double *f, double *fl, double *fr) /* MC */
{
  double alpha=2.0;
  double df=0.5*minmod(alpha*minmod(f[0]-f[-1],f[+1]-f[0]),0.5*(f[+1]-f[-1]));
  (*fl)=f[0]+df;
  (*fr)=f[0]-df;
}
void muscl_vl_cal_flr(const double *f, double *fl, double *fr) /* van Leer */
{
  double d[3]={f[0]-f[-1],f[+1]-f[0]};
  d[2]=d[0]*d[1];
  double df=(d[2]>0)?(d[2]/(d[0]+d[1])):0;
  (*fl)=f[0]+df;
  (*fr)=f[0]-df;
}
void muscl_kr_cal_flr(const double *f, double *fl, double *fr) /* Koren */
{
  double df1=2.0*(f[0]-f[-1]);
  double df2=2.0*(f[+1]-f[0]);
  double dfl=0.5*minmod3(df1,df2,(0.5*df1+df2)/3.0);
  double dfr=0.5*minmod3(df1,df2,(df1+0.5*df2)/3.0);
  (*fl)=f[0]+dfl;
  (*fr)=f[0]-dfr;
}

void muscl_skr4_cal_flr(const double *f, double *fl, double *fr) /* Smooth koren-type */
{
  const double df1=f[0]-f[-1];
  const double df2=f[+1]-f[0];

  // if (!((df1>0.0 && df2>0.0) || (df1<0.0 && df2<0.0))){
  if (df1*df2 <= 0){
    (*fl)=f[0];
    (*fr)=f[0];
    return;
  }

  const bool sw=(std::abs(df2)<=std::abs(df1));
  const double r=sw ? df2/df1 : df1/df2;
  const double df=sw ? df2 : df1;

  const double num1=1.0+r*(6.0+r*(-4.0+r*3.0));
  const double den1=1.0+r*(6.0+r*(6.0+r*(-4.0+r*3.0)));

  const double num2=3.0+r*(-4.0+r*(6.0+r));
  const double den2=3.0+r*(-4.0+r*(6.0+r*(6.0+r)));

  const double d1=df*(num1/den1);
  const double d2=df*(num2/den2);

  (*fl)=f[0]+(sw ? d1 : d2);
  (*fr)=f[0]-(sw ? d2 : d1);
}

/* WCNS routines */
/* f = address @ i, *fl = left state @ i+1/2, *fr = Right state @i-1/2 */
void wcns3_cal_flr(const double *f, double *fl, double *fr)
/* 3rd order WCNS interpolation (Minoshima+19) */
{
  double is[2],w[2],c[2];
  double eps=1e-40,c1,tau3;

  /* Smoothness indicator */
  c1=f[-1]-f[0];
  is[0]=c1*c1;
  c1=f[+1]-f[0];
  is[1]=c1*c1;

  /* Yamaleev+ 09 */
  tau3=f[+1]-2.0*f[0]+f[-1];
  tau3=tau3*tau3;

  is[0]=(1.0+tau3/(is[0]+eps));
  is[1]=(1.0+tau3/(is[1]+eps));

  /* Linear weight (3rd order) */
  c[0]=0.25;
  c[1]=0.75;

  /* Left-side value at 1.5 */
  w[0]=c[0]*is[0];
  w[1]=c[1]*is[1];
  w[0]/=(w[0]+w[1]);
  w[1]=1.0-w[0];
  *fl=0.5*(w[0]*(-f[-1]+3.0*f[0])+w[1]*(f[0]+f[+1]));

  /* Right-side value at 0.5 */
  w[0]=c[1]*is[0];
  w[1]=c[0]*is[1];
  w[0]/=(w[0]+w[1]);
  w[1]=1.0-w[0];
  *fr=0.5*(w[1]*(-f[+1]+3.0*f[0])+w[0]*(f[0]+f[-1]));
}
void wcns4_cal_flr(const double *f, double *fl, double *fr)
/* 4th order WCNS interpolation (Minoshima+19) */
{
  double is[3],w[3],c[3];
  double eps=1e-40,c1,c2,denom;

  /* Smoothness indicator */
  c1=f[-2]-4.0*f[-1]+3.0*f[0];
  c2=f[-2]-2.0*f[-1]+f[0];
  is[0]=(13./12.)*c2*c2+0.25*c1*c1;
  c1=f[-1]-f[+1];
  c2=f[-1]-2.0*f[0]+f[+1];
  is[1]=(13./12.)*c2*c2+0.25*c1*c1;
  c1=f[+2]-4.0*f[+1]+3.0*f[0];
  c2=f[+2]-2.0*f[+1]+f[0];
  is[2]=(13./12.)*c2*c2+0.25*c1*c1;
  is[0]=1.0/(is[0]*is[0]+eps);
  is[1]=1.0/(is[1]*is[1]+eps);
  is[2]=1.0/(is[2]*is[2]+eps);

  /* Linear weight, achieving 5th order of dfdx=(27*(fi+1/2-fi-1/2)-(fi+3/2-fi-3/2))/24 */
  c[0]=0.075;
  c[1]=0.650;
  c[2]=0.275;

  /* Left-side value at 2.5 */
  w[0]=c[0]*is[0];
  w[1]=c[1]*is[1];
  w[2]=c[2]*is[2];
  denom=1.0/(w[0]+w[1]+w[2]);
  w[0]*=denom;
  w[1]*=denom;
  w[2]=1.0-w[0]-w[1];
  *fl=0.125*(+w[0]*(3.0*f[-2]-10.0*f[-1]+15.0*f[0])
	     +w[1]*(-f[-1]+6.0*f[0]+3.0*f[+1])
	     +w[2]*(+3.0*f[0]+6.0*f[+1]-f[+2]));

  /* Right-side value at 1.5 */
  w[0]=c[2]*is[0];
  w[1]=c[1]*is[1];
  w[2]=c[0]*is[2];
  denom=1.0/(w[0]+w[1]+w[2]);
  w[0]*=denom;
  w[1]*=denom;
  w[2]=1.0-w[0]-w[1];
  *fr=0.125*(+w[2]*(3.0*f[+2]-10.0*f[+1]+15.0*f[0])
	     +w[1]*(-f[+1]+6.0*f[0]+3.0*f[-1])
	     +w[0]*(+3.0*f[0]+6.0*f[-1]-f[-2]));
}
