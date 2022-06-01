#include "emhd_func1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "common_interp.h"
#include "mhd_fd1d.h"
#include "boundary.h"

/* [IMPORTANT] Following functions relate electric field Enew and its base value Eorg in extended MHD */
/* Their relation is (ro-de^2 \nabla^2) Enew = ro*Eorg */
/* ro: number density (normalized by ro0) */
/* de: electron inertia length in density of ro0 */
/* Eorg: electric field WITHOUT the electron inertia effect (e.g, Eorg=-vxB) */
/* Enew: electric field WITH the electron inertia effect */
/* Reference: Amano, 2015, JCP */

/* Spatial difference */
double (*df2[])(const double *f)={cal_d2f_2nd,cal_d2f_2nd,cal_d2f_4th,cal_d2f_4th};

void emhd_enew2eorg(double *enew, double *eorg, const double *ro,
		    double de, double dx, int nx, int xoff, int dnx)
/* Convert Enew => Eorg */
/* ro is number density and de is electron inertia length (both are normalized) */
/* Boundary condition included */
{
  const double fac=(de*de)/(dx*dx);
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={eorg,enew};
  int dnxs[]={dnx};

  /* Boundary condition for Enew */
  boundary(&p[1],1,nx,xoff,dnxs);

  for (int i=xoff;i<nx-xoff;i++){
    int ss=i;
    eorg[ss]=enew[ss]-fac*func_d2f(&enew[ss])/ro[ss];
  }
  
  /* Boundary condition for Eorg (output) */
  boundary(&p[0],1,nx,xoff,dnxs);
}

void emhd_eorg2enew(double *eorg, double *enew, const double *ro,
		    double de, double dx, int nx, int xoff, int dnx)
/* Wrapper of emhd_eorg2enew functions */
{
  int cnt;
  cnt=emhd_eorg2enew_cg(eorg,enew,ro,de,dx,nx,xoff,dnx);
  if (!cnt) puts("emhd_eorg2enew: Not converged!!");
}

int emhd_eorg2enew_cg(double *eorg, double *enew, const double *ro,
		      double de, double dx, int nx, int xoff, int dnx)
/* Convert Eorg => Enew with Conjugate Gradient method */
/* ro is number density and de is electron inertia length (both are normalized) */
/* Boundary condition included */
{
  int i,ss;
  int cnt=0;
  const int cmax=nx;
  const double eps=1e-8;
  const double fac=(de*de)/(dx*dx);
  double numer,denom,coef,anorm,anormf=0.0;
  double *rk,*pk,*ap;
  rk=(double*)malloc(sizeof(double)*nx);
  pk=(double*)malloc(sizeof(double)*nx);
  ap=(double*)malloc(sizeof(double)*nx);
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={eorg,enew,pk};
  int dnxs[]={dnx};

  /* Initialize */
  boundary(&p[1],1,nx,xoff,dnxs); /* Boundary condition for Enew */
  numer=0.0;
  for (i=xoff;i<nx-xoff;i++){
    ss=i;
    rk[ss]=ro[ss]*eorg[ss]-(ro[ss]*enew[ss]-fac*func_d2f(&enew[ss]));
    pk[ss]=rk[ss];
    numer+=rk[ss]*rk[ss];
    anormf+=fabs(rk[ss]);
  }
  
  /* Iteration */
  do{

    boundary(&p[2],1,nx,xoff,dnxs); /* Boundary condition for Pk */
    denom=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      ap[ss]=(ro[ss]*pk[ss]-fac*func_d2f(&pk[ss]));
      denom+=pk[ss]*ap[ss];
    }
    coef=(denom == 0)?0:numer/denom;

    anorm=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      enew[ss]+=coef*pk[ss];
      rk[ss]-=coef*ap[ss];
      anorm+=fabs(rk[ss]);
    }
    
    denom=numer;
    numer=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      numer+=rk[ss]*rk[ss];
    }
    coef=(denom == 0)?0:numer/denom;

    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      pk[ss]=rk[ss]+coef*pk[ss];
    }
    
  } while( (anorm > eps*anormf) && (++cnt < cmax) );

  /* Boundary condition for Enew (output) */
  boundary(&p[1],1,nx,xoff,dnxs);

  free(rk);
  free(pk);
  free(ap);

  return cmax-cnt;		/* If zero, iteration does NOT converge */
}
