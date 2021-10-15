#include "mhd_fd1d.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mhd_common.h"
#include "boundary.h"

/* [IMPORTANT] Macros are defined in mhd_fd1d.h to select Riemann solver, spatial and temporal accuracy */

/* Riemann solvers */
void (*riemann[])(double, double, double, double, double, double, double, 
		  double, double, double, double, double, double, double, 
		  double, double, const double*,
		  double*, double*, double*, double*, double*, double*, double*)={calc_flux_roe,calc_flux_hlld,calc_flux_lhlld,calc_flux_mlau};
/* Interpolation function */
void (*interpol[])(const double *f, double *fl, double *fr)={cal_flr_1st,muscl_mm_cal_flr,wcns3_cal_flr,wcns4_cal_flr};
/* Spatial difference */
double (*df1[])(const double *f)={cal_df_2nd,cal_df_2nd,cal_df_4th,cal_df_4th};


void mhd_fd1d(double *ro, double *mx, double *my, double *mz,
	      double *by, double *bz, double *en,
	      double bx, double dt, double dx,
	      int nx, int xoff, double gamma)
/* 1D MHD simulation  */
/* ro: density */
/* m?: moment */
/* en: energy normalized by B^2/(4 pi) */
/* b?: magnetic field */
/* gamma: specific heat ratio */
{
  int i,ss,rk;
  const int nm=7;		/* Number of variables */
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const double dtdx=dt/dx;
  void (*func_flux)(double, double, double, double, double, double, double, 
		    double, double, double, double, double, double, double, 
		    double, double, const double*,
		    double*, double*, double*, double*, double*, double*, double*)=riemann[RMN];
  void (*func_lr)(const double *f, double *fl, double *fr)=interpol[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ut,*ul,*ur,*fx,*vx,*vy,*vz,*pr;
  double *p[]={ro,mx,my,mz,by,bz,en};
  
  ut=(double*)malloc(sizeof(double)*nm*nx);
  ul=(double*)malloc(sizeof(double)*nm*nx);
  ur=(double*)malloc(sizeof(double)*nm*nx);
  fx=(double*)malloc(sizeof(double)*nm*nx);
  vx=(double*)malloc(sizeof(double)*nx);
  vy=(double*)malloc(sizeof(double)*nx);
  vz=(double*)malloc(sizeof(double)*nx);
  pr=(double*)malloc(sizeof(double)*nx);

  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nx],p[i],nx);
  }
  
  /* Runge-kutta stage */
  for (rk=0;rk<R_K;rk++){

#ifdef _OPENMP
#pragma omp parallel private(i,ss)
#endif
    {
      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=0;i<nx;i++){
	ss=i;
	mhd_prmtv(ro[ss],mx[ss],my[ss],mz[ss],bx,by[ss],bz[ss],en[ss],
		  &vx[ss],&vy[ss],&vz[ss],&pr[ss],gamma);
      }
      /* Primitive variable at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=2;i<nx-2;i++){
	ss=i;
	int sl=i+1,sr=i;
	double vl[7],vr[7];
	mhd_lrstate(&ro[ss],&vx[ss],&vy[ss],&vz[ss],&by[ss],&bz[ss],&pr[ss],
		    bx,gamma,1,func_lr,vl,vr);
	/* Left-face @ i+1/2 */
	ul[nm*sl+0]=vl[0];	/* ro */
	ul[nm*sl+1]=vl[1];	/* vx */
	ul[nm*sl+2]=vl[2];	/* vy */
	ul[nm*sl+3]=vl[3];	/* vz */
	ul[nm*sl+4]=vl[4];	/* by */
	ul[nm*sl+5]=vl[5];	/* bz */
	ul[nm*sl+6]=vl[6];	/* pr */
	/* Right-face @ i-1/2 */
	ur[nm*sr+0]=vr[0];	/* ro */
	ur[nm*sr+1]=vr[1];	/* vx */
	ur[nm*sr+2]=vr[2];	/* vy */
	ur[nm*sr+3]=vr[3];	/* vz */
	ur[nm*sr+4]=vr[4];	/* by */
	ur[nm*sr+5]=vr[5];	/* bz */
	ur[nm*sr+6]=vr[6];	/* pr */
      }
      /* Numerical flux at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=3;i<nx-2;i++){
	ss=i;	
	double flux[nm]={0};
	double dvsd[2]={(vx[i]-vx[i-1]),0.0};
	func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+4],ul[nm*ss+5],ul[nm*ss+6],
		  ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+4],ur[nm*ss+5],ur[nm*ss+6],
		  bx,gamma,dvsd,
		  &flux[0],&flux[1],&flux[2],&flux[3],&flux[4],&flux[5],&flux[6]);
	fx[nm*ss+0]=flux[0];	/* ro */
	fx[nm*ss+1]=flux[1];	/* mx */
	fx[nm*ss+2]=flux[2];	/* my */
	fx[nm*ss+3]=flux[3];	/* mz */
	fx[nm*ss+4]=flux[4];	/* by */
	fx[nm*ss+5]=flux[5];	/* bz */
	fx[nm*ss+6]=flux[6];	/* en */
      }
      /* Open boundary condition */
      /* Set here to advect boundary data @ t=n to t=n+1 */
      if (rk == 0){             /* Call only for 1st RK step to set BC @ t=n+1 */
#ifdef _OPENMP
#pragma omp single
#endif
	{
	  boundary(ro,mx,my,mz,by,bz,en,bx,dt,dx,nx,xoff,gamma,0,1);
	}
      }
      
      /* Update */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=xoff;i<nx-xoff;i++){
	ss=i;
	double *val[]={&ro[ss],&mx[ss],&my[ss],&mz[ss],&by[ss],&bz[ss],&en[ss]};
	double val0[]={ut[0*nx+ss],ut[1*nx+ss],ut[2*nx+ss],ut[3*nx+ss],ut[4*nx+ss],ut[5*nx+ss],ut[6*nx+ss]};
	mhd_updt1d(val,val0,&fx[nm*ss+0],dtdx,rk_fac[rk],func_df);
      }
      
    } /* OpenMP */

  }

  free(ut);
  free(ul);
  free(ur);
  free(fx);
  free(vx);
  free(vy);
  free(vz);
  free(pr);
}
