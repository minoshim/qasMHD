#include "mhd_fd2d.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mhd_common.h"
#include "boundary.h"

/* [IMPORTANT] Macros are defined in mhd_fd2d.h to select Riemann solver, spatial and temporal accuracy, CUCT flag. */

/* Riemann solvers */
void (*riemann[])(double, double, double, double, double, double, double,
		  double, double, double, double, double, double, double,
		  double, double, const double*,
		  double*, double*, double*, double*, double*, double*, double*)={calc_flux_roe,calc_flux_hlld,calc_flux_lhlld,calc_flux_mlau};
/* Linear interpolation function to face LR states*/
void (*l_interp[])(const double *f, double *fl, double *fr)={cal_flr_1st,cal_flr_2nd,cal_flr_3rd,cal_flr_4th};
/* Interpolation function to face LR states */
void (*interpol[])(const double *f, double *fl, double *fr)={cal_flr_1st,muscl_mm_cal_flr,wcns3_cal_flr,wcns4_cal_flr};
/* Interpolation from face to center */
double (*fcen[])(const double *f)={cal_fcen_2nd,cal_fcen_2nd,cal_fcen_4th,cal_fcen_4th};
/* Spatial difference */
double (*df1[])(const double *f)={cal_df_2nd,cal_df_2nd,cal_df_4th,cal_df_4th};

/* Boundary condition flag, defined in global.hpp */
extern int dnxs[8];
extern int dnys[8];
/* Staggered grid flag, defined in global.hpp */
extern int stxs[8];
extern int stys[8];

void mhd_fd2d(double *p[], double dt, double dx, double dy,
	      int nm, int nx, int ny, int xoff, int yoff, double gamma,
	      const double *phi_g,
	      int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D finite-difference code for MHD */
/* 2nd, 3rd, and 4th accuracy in space */
/* ROE/HLLD/LHLLD/MLAU + Central-Upwind-CT */
/* Need four offsets for boundary */

/* p[0..7] = ro,mx,my,mz,bx,by,bz,en */
/* nm = Number of variables (8 in multi-D) */
/* gamma: specific heat ratio */

/* phi_g is gravitational potential */
{
  int i,j,ss,rk;
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double dtdx=dt/dx,dtdy=dt/dy;
  const double cflg[]={1,1,1,1,0,0,1,1}; /* Flag for cell center variable */
  void (*func_flux)(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, const double*,
		    double*, double*, double*, double*, double*, double*, double*)=riemann[RMN];
  void (*lfun_lr)(const double *f, double *fl, double *fr)=l_interp[ODR-1];
  void (*func_lr)(const double *f, double *fl, double *fr)=interpol[ODR-1];
  double (*func_bc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ro=p[0],*mx=p[1],*my=p[2],*mz=p[3],*bx=p[4],*by=p[5],*bz=p[6],*en=p[7];
  double *ut,*ul,*ur,*ql,*qr,*fx,*fy,*fc;
  double *vx,*vy,*vz,*pr,*cx,*cy,*cz=bz,*ez;
  double *ct,*dvx,*dvy;

  ut=(double*)malloc(sizeof(double)*nm*nxy);
  ul=(double*)malloc(sizeof(double)*nm*nxy);
  ur=(double*)malloc(sizeof(double)*nm*nxy);
  ql=(double*)malloc(sizeof(double)*nxy);
  qr=(double*)malloc(sizeof(double)*nxy);
  fx=(double*)malloc(sizeof(double)*nm*nxy);
  fy=(double*)malloc(sizeof(double)*nm*nxy);
  fc=(double*)malloc(sizeof(double)*nxy);
  vx=(double*)malloc(sizeof(double)*nxy);
  vy=(double*)malloc(sizeof(double)*nxy);
  vz=(double*)malloc(sizeof(double)*nxy);
  pr=(double*)malloc(sizeof(double)*nxy);
  cx=(double*)malloc(sizeof(double)*nxy);
  cy=(double*)malloc(sizeof(double)*nxy);
  ez=(double*)malloc(sizeof(double)*nxy);
  ct=(double*)malloc(sizeof(double)*nxy);
  dvx=(double*)malloc(sizeof(double)*nxy);
  dvy=(double*)malloc(sizeof(double)*nxy);

  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nxy],p[i],nxy);
  }

  /* Runge-Kutta stage */
  for (rk=0;rk<R_K;rk++){
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,ss)
#endif
    {

      /* Bx and By at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=1;i<nx-2;i++){
	  double val[]={bx[nx*j+(i-1)],bx[nx*j+(i+0)],bx[nx*j+(i+1)],bx[nx*j+(i+2)]};
	  cx[nx*j+i]=func_bc(&val[1]);
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  double val[]={by[nx*(j-1)+i],by[nx*(j+0)+i],by[nx*(j+1)+i],by[nx*(j+2)+i]};
	  cy[nx*j+i]=func_bc(&val[1]);
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[2]={0,0},sty[2]={0,0};
	double *pc[]={cx,cy};
	boundary(pc,2,nx,ny,xoff,yoff,stx,&dnxs[4],sty,&dnys[4],gamma,phi_g,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxy;ss++){
	mhd_prmtv(ro[ss],mx[ss],my[ss],mz[ss],cx[ss],cy[ss],cz[ss],en[ss]-ro[ss]*phi_g[ss],
		  &vx[ss],&vy[ss],&vz[ss],&pr[ss],gamma); /* Subtract gravitational potential */
	ez[ss]=0.0;		/* Necessary initialize at cell corner */
	ct[ss]=0.5;
      }

      /* CUCT upwind weighting (x-y plane for Ez) */
#if (CTW)
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny;j++){
#pragma simd
	for (i=1;i<nx;i++){
	  ss=nx*j+i;
	  ct[ss]=mhd_cuct_weight(&ro[ss],&vx[ss],&vy[ss],&bx[ss],&by[ss],1,nx);
	}
      }
#endif
      
      /* dvx and dvy at cell center for shock detection */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-1;j++){
#pragma simd
	for (i=1;i<nx-1;i++){
	  ss=nx*j+i;
	  dvx[ss]=min((vx[nx*j+(i+1)]-vx[nx*j+i]),(vx[nx*j+i]-vx[nx*j+(i-1)]));
	  dvy[ss]=min((vy[nx*(j+1)+i]-vy[nx*j+i]),(vy[nx*j+i]-vy[nx*(j-1)+i]));
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[2]={0,0},sty[2]={0,0};
	int dnx[2]={-dnxs[1],dnxs[2]},dny[2]={dnys[1],-dnys[2]};
	double *pv[]={dvx,dvy};
	boundary(pv,2,nx,ny,xoff,yoff,stx,dnx,sty,dny,gamma,phi_g,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Primitive variable at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int sl,sr;
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  double vl[7],vr[7];
	  mhd_lrstate(&ro[ss],&vx[ss],&vy[ss],&vz[ss],&cy[ss],&cz[ss],&pr[ss],
		      cx[ss],gamma,1,func_lr,vl,vr);
	  /* Left-face @ i+1/2 */
	  ul[nm*sl+0]=vl[0];	/* ro */
	  ul[nm*sl+1]=vl[1];	/* vx */
	  ul[nm*sl+2]=vl[2];	/* vy */
	  ul[nm*sl+3]=vl[3];	/* vz */
	  ul[nm*sl+4]=bx[sl];	/* bx */
	  ul[nm*sl+5]=vl[4];	/* by */
	  ul[nm*sl+6]=vl[5];	/* bz */
	  ul[nm*sl+7]=vl[6];	/* pr */
	  /* Right-face @ i-1/2 */
	  ur[nm*sr+0]=vr[0];	/* ro */
	  ur[nm*sr+1]=vr[1];	/* vx */
	  ur[nm*sr+2]=vr[2];	/* vy */
	  ur[nm*sr+3]=vr[3];	/* vz */
	  ur[nm*sr+4]=bx[sr];	/* bx */
	  ur[nm*sr+5]=vr[4];	/* by */
	  ur[nm*sr+6]=vr[5];	/* bz */
	  ur[nm*sr+7]=vr[6];	/* pr */
	}
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  double vl,vr;
	  /* Linear interpolation of numerical flux of By */
	  mhd_lr_fb(&vx[ss],&vy[ss],&bx[ss],&cy[ss],1,lfun_lr,&vl,&vr);
	  ql[sl]=vl;		/* by*vx-bx*vy @ i+1/2 Left */
	  qr[sr]=vr;		/* by*vx-bx*vy @ i-1/2 Right */
	}
      }

      /* Numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  double flux[8]={0},bn=0.5*(ul[nm*ss+4]+ur[nm*ss+4]);
	  double dvsd[2]={(vx[nx*j+i]-vx[nx*j+(i-1)]),min(dvy[nx*j+(i-1)],dvy[nx*j+i])};
	  func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+5],ul[nm*ss+6],ul[nm*ss+7],
		    ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+5],ur[nm*ss+6],ur[nm*ss+7],
		    bn,gamma,dvsd,
		    &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);
	  
	  fx[nm*ss+0]=flux[0];	/* ro */
	  fx[nm*ss+1]=flux[1];	/* mx */
	  fx[nm*ss+2]=flux[2];	/* my */
	  fx[nm*ss+3]=flux[3];	/* mz */
	  fx[nm*ss+4]=0;	/* bx */
	  fx[nm*ss+5]=flux[5];	/* by */
	  fx[nm*ss+6]=flux[6];	/* bz */
	  fx[nm*ss+7]=flux[7];	/* en */
	  fx[nm*ss+7]+=flux[0]*0.5*(phi_g[ss- 1]+phi_g[ss]); /* Work done by gravity */
	  /* Split central and upwind parts in numerical flux of By */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fx[nm*ss+5]-=fc[ss];	      /* Upwind part */
	}
      }

      /* Numerical flux of By at cell corner along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sl,sr;
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  double vl[2],vr[2];
	  mhd_lr_single(&fx[nm*ss+5],nm*nx,lfun_lr,&vl[0],&vr[0]); /* Upwind part */
	  mhd_lr_single(&fc[ss],nx,lfun_lr,&vl[1],&vr[1]); /* Central part */
	  ul[2*sl+0]=vl[0];
	  ur[2*sr+0]=vr[0];
	  ul[2*sl+1]=vl[1];
	  ur[2*sr+1]=vr[1];
	}
      }

      /* E-field at cell corner along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=-0.5*((ul[2*ss+0]+ur[2*ss+0])+(1.0-ct[ss])*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* Primitive variable at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-2;j++){
	int sl,sr;
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  double vl[7],vr[7];
	  mhd_lrstate(&ro[ss],&vy[ss],&vz[ss],&vx[ss],&cz[ss],&cx[ss],&pr[ss],
		      cy[ss],gamma,nx,func_lr,vl,vr);
	  /* Left-face @ j+1/2 */
	  ul[nm*sl+0]=vl[0];	/* ro */
	  ul[nm*sl+1]=vl[1];	/* vy */
	  ul[nm*sl+2]=vl[2];	/* vz */
	  ul[nm*sl+3]=vl[3];	/* vx */
	  ul[nm*sl+4]=by[sl];	/* by */
	  ul[nm*sl+5]=vl[4];	/* bz */
	  ul[nm*sl+6]=vl[5];	/* bx */
	  ul[nm*sl+7]=vl[6];	/* pr */
	  /* Right-face @ j-1/2 */
	  ur[nm*sr+0]=vr[0];	/* ro */
	  ur[nm*sr+1]=vr[1];	/* vy */
	  ur[nm*sr+2]=vr[2];	/* vz */
	  ur[nm*sr+3]=vr[3];	/* vx */
	  ur[nm*sr+4]=by[sr];	/* by */
	  ur[nm*sr+5]=vr[4];	/* bz */
	  ur[nm*sr+6]=vr[5];	/* bx */
	  ur[nm*sr+7]=vr[6];	/* pr */
	}
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  sl=nx*(j+1)+i;
	  sr=nx*j+i;
	  double vl,vr;
	  /* Linear interpolation of numerical flux of Bx */
	  mhd_lr_fb(&vy[ss],&vx[ss],&by[ss],&cx[ss],nx,lfun_lr,&vl,&vr);
	  ql[sl]=vl;		/* bx*vy-by*vx @ j+1/2 Left */
	  qr[sr]=vr;		/* bx*vy-by*vx @ j-1/2 Right */
	}
      }

      /* Numerical flux at cell face along Y */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  double flux[8]={0},bn=0.5*(ul[nm*ss+4]+ur[nm*ss+4]);
	  double dvsd[2]={(vy[nx*j+i]-vy[nx*(j-1)+i]),min(dvx[nx*(j-1)+i],dvx[nx*j+i])};
	  func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+5],ul[nm*ss+6],ul[nm*ss+7],
		    ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+5],ur[nm*ss+6],ur[nm*ss+7],
		    bn,gamma,dvsd,
		    &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	  fy[nm*ss+0]=flux[0];	/* ro */
	  fy[nm*ss+2]=flux[1];	/* my */
	  fy[nm*ss+3]=flux[2];	/* mz */
	  fy[nm*ss+1]=flux[3];	/* mx */
	  fy[nm*ss+5]=0;	/* by */
	  fy[nm*ss+6]=flux[5];	/* bz */
	  fy[nm*ss+4]=flux[6];	/* bx */
	  fy[nm*ss+7]=flux[7];	/* en */
	  fy[nm*ss+7]+=flux[0]*0.5*(phi_g[ss-nx]+phi_g[ss]); /* Work done by gravity */
	  /* Split central and upwind parts in numerical flux of Bx */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fy[nm*ss+4]-=fc[ss];	      /* Upwind part */
	}
      }

      /* Numerical flux of Bx at cell corner along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
	int sl,sr;
#pragma simd
	for (i=2;i<nx-2;i++){
	  ss=nx*j+i;
	  sl=nx*j+(i+1);
	  sr=nx*j+i;
	  double vl[2],vr[2];
	  mhd_lr_single(&fy[nm*ss+4],nm,lfun_lr,&vl[0],&vr[0]); /* Upwind part */
	  mhd_lr_single(&fc[ss],1,lfun_lr,&vl[1],&vr[1]);	/* Central part */
	  ul[2*sl+0]=vl[0];
	  ur[2*sr+0]=vr[0];
	  ul[2*sl+1]=vl[1];
	  ur[2*sr+1]=vr[1];
	}
      }

      /* E-field at cell corner along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
#pragma simd
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=+0.5*((ul[2*ss+0]+ur[2*ss+0])+ct[ss]*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* Update variables at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  /* Note: cx and cy are unchanged (dummy) since cflg[4]=cflg[5]=0 */
	  double *val[]={&ro[ss],&mx[ss],&my[ss],&mz[ss],&cx[ss],&cy[ss],&cz[ss],&en[ss]};
	  double val0[]={ut[0*nxy+ss],ut[1*nxy+ss],ut[2*nxy+ss],ut[3*nxy+ss],ut[4*nxy+ss],ut[5*nxy+ss],ut[6*nxy+ss],ut[7*nxy+ss]};
	  mhd_updt2d(val,val0,&fx[nm*ss+0],&fy[nm*ss+0],dtdx,dtdy,rk_fac[rk],cflg,1,nx,func_df);
	  mx[ss]+=-rk_fac[rk][1]*ro[ss]*0.5*(phi_g[ss+ 1]-phi_g[ss- 1])*dtdx; /* Work done by gravity */
	  my[ss]+=-rk_fac[rk][1]*ro[ss]*0.5*(phi_g[ss+nx]-phi_g[ss-nx])*dtdy; /* Work done by gravity */
	}
      }
      /* Update CT Bx */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
#pragma simd
	for (i=xoff;i<nx-xoff+1;i++){
	  ss=nx*j+i;
	  mhd_updt2d_ctb(&bx[ss],ut[4*nxy+ss],&ez[ss],+dtdy,rk_fac[rk],nx,func_df);
	}
      }
      /* Update CT By */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff+1;j++){
#pragma simd
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  mhd_updt2d_ctb(&by[ss],ut[5*nxy+ss],&ez[ss],-dtdx,rk_fac[rk], 1,func_df);
	}
      }

    } /* OpenMP */

    /* Boundary condition */
    boundary(p,nm,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,gamma,phi_g,mpi_rank,mpi_numx,mpi_numy);
  }
  
  free(ut);
  free(ul);
  free(ur);
  free(ql);
  free(qr);
  free(fx);
  free(fy);
  free(fc);
  free(vx);
  free(vy);
  free(vz);
  free(pr);
  free(cx);
  free(cy);
  free(ez);
  free(ct);
  free(dvx);
  free(dvy);
}

