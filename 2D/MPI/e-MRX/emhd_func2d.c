#include "mpi.h"
#include "emhd_func2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "common_interp.h"
#include "mhd_fd2d.h"
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
		    double de, double dx, double dy,
		    int nx, int ny, int xoff, int yoff,
		    int stx, int dnx, int sty, int dny,
		    int mpi_rank, int mpi_numx, int mpi_numy)
/* Convert Enew => Eorg */
/* ro is number density and de is electron inertia length (both are normalized) */
/* Boundary condition included */
{
  const double fac[2]={(de*de)/(dx*dx),(de*de)/(dy*dy)};
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={eorg,enew};
  int stxs[]={stx};
  int dnxs[]={dnx};
  int stys[]={sty};
  int dnys[]={dny};
  
  /* Boundary condition for Enew */
  boundary(&p[1],1,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);

  for (int j=yoff;j<ny-yoff+sty;j++){
    for (int i=xoff;i<nx-xoff+stx;i++){
      int ss=nx*j+i;
      double val0[]={enew[ss-2   ],enew[ss-1   ],enew[ss],enew[ss+1   ],enew[ss+2   ]};
      double val1[]={enew[ss-2*nx],enew[ss-1*nx],enew[ss],enew[ss+1*nx],enew[ss+2*nx]};
      eorg[ss]=enew[ss]-(fac[0]*func_d2f(&val0[2])+fac[1]*func_d2f(&val1[2]))/ro[ss];
    }
  }

  /* Boundary condition for Eorg (output) */
  boundary(&p[0],1,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);
}

void emhd_eorg2enew(double *eorg, double *enew, const double *ro,
		    double de, double dx, double dy,
		    int nx, int ny, int xoff, int yoff,
		    int stx, int dnx, int sty, int dny,
		    int mpi_rank, int mpi_numx, int mpi_numy)
/* Wrapper of emhd_eorg2enew functions */
{
  int cnt;
  cnt=emhd_eorg2enew_cg(eorg,enew,ro,de,dx,dy,nx,ny,xoff,yoff,stx,dnx,sty,dny,mpi_rank,mpi_numx,mpi_numy);
  if (!mpi_rank && !cnt) puts("emhd_eorg2enew: Not converged!!");
}

int emhd_eorg2enew_cg(double *eorg, double *enew, const double *ro,
		      double de, double dx, double dy,
		      int nx, int ny, int xoff, int yoff,
		      int stx, int dnx, int sty, int dny,
		      int mpi_rank, int mpi_numx, int mpi_numy)
/* Convert Eorg => Enew with Conjugate Gradient method */
/* ro is number density and de is electron inertia length (both are normalized) */
/* Boundary condition included */
{
  int i,j,ss;
  int cnt=0;
  const int cmax=nx*ny;
  const double eps=1e-8;
  const double fac[2]={(de*de)/(dx*dx),(de*de)/(dy*dy)};
  double numer,denom,coef,anorm,anormf=0.0;
  double *rk,*pk,*ap;
  rk=(double*)malloc(sizeof(double)*nx*ny);
  pk=(double*)malloc(sizeof(double)*nx*ny);
  ap=(double*)malloc(sizeof(double)*nx*ny);
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={eorg,enew,pk};
  int stxs[]={stx};
  int dnxs[]={dnx};
  int stys[]={sty};
  int dnys[]={dny};
  double dtmp[2],dtmpa[2];

  /* Initialize */
  boundary(&p[1],1,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy); /* B.C. for Enew */
  numer=0.0;
  for (j=yoff;j<ny-yoff+sty;j++){
    for (i=xoff;i<nx-xoff+stx;i++){
      ss=nx*j+i;
      double val0[]={enew[ss-2   ],enew[ss-1   ],enew[ss],enew[ss+1   ],enew[ss+2   ]};
      double val1[]={enew[ss-2*nx],enew[ss-1*nx],enew[ss],enew[ss+1*nx],enew[ss+2*nx]};
      rk[ss]=ro[ss]*eorg[ss]-( ro[ss]*enew[ss]-(fac[0]*func_d2f(&val0[2])+fac[1]*func_d2f(&val1[2])) );
      pk[ss]=rk[ss];
      numer+=rk[ss]*rk[ss];
      anormf+=fabs(rk[ss]);
    }
  }
  dtmp[0]=numer;
  dtmp[1]=anormf;
  MPI_Allreduce(dtmp,dtmpa,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  numer=dtmpa[0];
  anormf=dtmpa[1];

  /* Iteration */
  do{

    boundary(&p[2],1,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy); /* B.C. for Pk */
    denom=0.0;
    for (j=yoff;j<ny-yoff+sty;j++){
      for (i=xoff;i<nx-xoff+stx;i++){
	ss=nx*j+i;
	double val0[]={pk[ss-2   ],pk[ss-1   ],pk[ss],pk[ss+1   ],pk[ss+2   ]};
	double val1[]={pk[ss-2*nx],pk[ss-1*nx],pk[ss],pk[ss+1*nx],pk[ss+2*nx]};
	ap[ss]=( ro[ss]*pk[ss]-(fac[0]*func_d2f(&val0[2])+fac[1]*func_d2f(&val1[2])) );
	denom+=pk[ss]*ap[ss];
      }
    }
    dtmp[0]=denom;
    MPI_Allreduce(dtmp,dtmpa,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    denom=dtmpa[0];
    coef=(denom == 0)?0:numer/denom;

    anorm=0.0;
    for (j=yoff;j<ny-yoff+sty;j++){
      for (i=xoff;i<nx-xoff+stx;i++){
	ss=nx*j+i;
	enew[ss]+=coef*pk[ss];
	rk[ss]-=coef*ap[ss];
	anorm+=fabs(rk[ss]);
      }
    }
    dtmp[0]=anorm;
    MPI_Allreduce(dtmp,dtmpa,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    anorm=dtmpa[0];

    denom=numer;
    numer=0.0;
    for (j=yoff;j<ny-yoff+sty;j++){
      for (i=xoff;i<nx-xoff+stx;i++){
	ss=nx*j+i;
	numer+=rk[ss]*rk[ss];
      }
    }
    dtmp[0]=numer;
    MPI_Allreduce(dtmp,dtmpa,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    numer=dtmpa[0];
    coef=(denom == 0)?0:numer/denom;
    
    for (j=yoff;j<ny-yoff+sty;j++){
      for (i=xoff;i<nx-xoff+stx;i++){
	ss=nx*j+i;
	pk[ss]=rk[ss]+coef*pk[ss];
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
  } while( (anorm > eps*anormf) && (++cnt < cmax) );

  /* Boundary condition for Enew (output) */
  boundary(&p[1],1,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);

  free(rk);
  free(pk);
  free(ap);

  return cmax-cnt;		/* If zero, iteration does NOT converge */
}
