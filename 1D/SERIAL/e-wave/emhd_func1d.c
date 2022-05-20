#include "emhd_func1d.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "common_interp.h"
#include "mhd_fd1d.h"
#include "boundary.h"

/* [IMPORTANT] Following functions relate magnetic field B and its modification Q in Electron MHD */
/* Their relation is (1-de^2 \nabla^2) B = Q, where de is electron inertia length */

/* Spatial difference */
double (*df2[])(const double *f)={cal_d2f_2nd,cal_d2f_2nd,cal_d2f_4th,cal_d2f_4th};

void emhd_b2q(double *b, double *q, double de, double dx, int nx, int xoff, int dnx)
/* Convert B => Q. de is electron inertia length */
{
  const double fac=(de*de)/(dx*dx);
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={q,b};
  int dnxs[]={dnx};
  
  /* Boundary condition for B */
  boundary(&p[1],1,nx,xoff,dnxs);

  for (int i=xoff;i<nx-xoff;i++){
    int ss=i;
    q[ss]=b[ss]-fac*func_d2f(&b[ss]);
  }

  /* Boundary condition for Q (output) */
  boundary(&p[0],1,nx,xoff,dnxs);
}

void emhd_q2b(double *q, double *b, double de, double dx, int nx, int xoff, int dnx)
/* Wrapper of emhd_q2b functions */
{
  int cnt;
  cnt=emhd_q2b_cg(q,b,de,dx,nx,xoff,dnx);
  /* cnt=emhd_q2b_gs(q,b,de,dx,nx,xoff,dnx); */
}

int emhd_q2b_cg(double *q, double *b, double de, double dx, int nx, int xoff, int dnx)
/* Convert Q = > B. de is electron inertia length */
/* Conjugate Gradient method */
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
  double *p[]={q,b,pk};
  int dnxs[]={dnx};

  /* Initialize */
  boundary(&p[1],1,nx,xoff,dnxs); /* Boundary condition for B */
  numer=0.0;
  for (i=xoff;i<nx-xoff;i++){
    ss=i;
    rk[ss]=q[ss]-(b[ss]-fac*func_d2f(&b[ss]));
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
      ap[ss]=(pk[ss]-fac*func_d2f(&pk[ss]));
      denom+=pk[ss]*ap[ss];
    }
    coef=(denom == 0)?0:numer/denom;

    anorm=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      b[ss] +=coef*pk[ss];
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

  /* Boundary condition for B (output) */
  boundary(&p[1],1,nx,xoff,dnxs);

  free(rk);
  free(pk);
  free(ap);

  return cnt;
}

int emhd_q2b_gs(double *q, double *b, double de, double dx, int nx, int xoff, int dnx)
/* Convert Q = > B. de is electron inertia length */
/* Gauss-Seidel method (use for debug) */
{
  int i,ss,pass;
  int cnt=0;
  const int cmax=nx;
  const int pwidth=(ODR >= 3)?3:2; /* Stencil half width in 2nd derivative */
  const double eps=1e-8;
  const double fac=(de*de)/(dx*dx);
  const double coef=(ODR >= 3)?2.5:2.0; /* Coefficient of B_i in 2nd derivative */
  const double denom=1.0/(1.0+coef*fac);
  double anorm,anormf=0.0;
  double numer,resid;
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={q,b};
  int dnxs[]={dnx};

  for (i=xoff;i<nx-xoff;i++){
    ss=i;
    anormf+=fabs(q[ss]-b[ss]);
  }

  /* Boundary condition for B */
  boundary(&p[1],1,nx,xoff,dnxs);
  
  do{
    anorm=0.0;

    for (pass=0;pass<pwidth;pass++){	/* Red-Black */
/* #ifdef _OPENMP */
/* #pragma omp parallel for private(i,ss,numer,resid) reduction(+:anorm) */
/* #endif */
      for (i=xoff+pass;i<nx-xoff;i+=pwidth){
	ss=i;
	numer=q[ss]+fac*(func_d2f(&b[ss])+coef*b[ss]);
	resid=numer*denom-b[ss];
	anorm+=fabs(resid);
	b[ss]+=resid;
      }

      /* Boundary condition for B */
      boundary(&p[1],1,nx,xoff,dnxs);
    }
  } while( (anorm > eps*anormf) && (++cnt < cmax) );

  return cnt;
}
  
