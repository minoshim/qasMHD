#include "mhd_fd2d.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mhd_common.h"
#include "emhd_func2d.h"
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
void (*interpol[])(const double *f, double *fl, double *fr)={cal_flr_1st,muscl_kr_cal_flr,wcns3_cal_flr,wcns4_cal_flr};
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

void mhd_fd2d(double *p[], double dt, double dx, double dy, double de, double mpme,
	      int nm, int nx, int ny, int xoff, int yoff, double gamma,
	      int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D finite-difference code for MHD */
/* 2nd, 3rd, and 4th accuracy in space */
/* ROE/HLLD/LHLLD/MLAU + Central-Upwind-CT */
/* Need four offsets for boundary */

/* p[0..7] = ro,mx,my,mz,bx,by,bz,en */
/* nm = Number of variables (8 in multi-D) */
/* gamma: specific heat ratio */


/* Include modification for extended MHD. de is electron inertia length */
/* mpme is mass ratio (proton/electron) */
/* 2022/06/17: Hall term and Pressure gradient term are NOT included (valid when mpme=1.0) */
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

  /* Modification for extended MHD */
  const double idx=1.0/dx,idy=1.0/dy;
  const double memt=1.0/(1.0+mpme);
  const double mpmt=1.0-memt;
  const double smemt=sqrt(memt);
  const double smpmt=sqrt(mpmt);
  const double smtme=1.0/smemt;
  double *rtmp,*enew,*vp,*ve,*dnvv;
  rtmp=(double*)malloc(sizeof(double)*nxy);
  enew=(double*)malloc(sizeof(double)*nxy*2);
  vp=(double*)malloc(sizeof(double)*nxy*3);
  ve=(double*)malloc(sizeof(double)*nxy*3);
  dnvv=(double*)malloc(sizeof(double)*nxy*3);
  
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
	boundary(pc,2,nx,ny,xoff,yoff,stx,&dnxs[4],sty,&dnys[4],mpi_rank,mpi_numx,mpi_numy);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxy;ss++){
	mhd_prmtv(ro[ss],mx[ss],my[ss],mz[ss],cx[ss],cy[ss],cz[ss],en[ss],
		  &vx[ss],&vy[ss],&vz[ss],&pr[ss],gamma);
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
	boundary(pv,2,nx,ny,xoff,yoff,stx,dnx,sty,dny,mpi_rank,mpi_numx,mpi_numy);
      }

      /* Modification for extended MHD */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-1;j++){
#pragma simd
	for (i=1;i<nx-1;i++){
	  ss=nx*j+i;
	  double fac=smtme*de/ro[ss];
	  double jj[3];
	  jj[0]=calc_jcen(&cy[ss],&cz[ss],idy,  0,nx, 0); /* jx @ cell center */
	  jj[1]=calc_jcen(&cz[ss],&cx[ss],  0,idx, 0, 1); /* jy */
	  jj[2]=calc_jcen(&cx[ss],&cy[ss],idx,idy, 1,nx); /* jz */
	  vp[0*nxy+ss]=vx[ss]+memt*fac*jj[0]; /* Proton velocity @ cell center */
	  vp[1*nxy+ss]=vy[ss]+memt*fac*jj[1];
	  vp[2*nxy+ss]=vz[ss]+memt*fac*jj[2];
	  ve[0*nxy+ss]=vx[ss]-mpmt*fac*jj[0]; /* Electron velocity @ cell center*/
	  ve[1*nxy+ss]=vy[ss]-mpmt*fac*jj[1];
	  ve[2*nxy+ss]=vz[ss]-mpmt*fac*jj[2];
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	double *pvp[]={&vp[0*nxy],&vp[1*nxy],&vp[2*nxy]};
	double *pve[]={&ve[0*nxy],&ve[1*nxy],&ve[2*nxy]};
	boundary(pvp,3,nx,ny,xoff,yoff,&stxs[1],&dnxs[1],&stys[1],&dnys[1],mpi_rank,mpi_numx,mpi_numy);
	boundary(pve,3,nx,ny,xoff,yoff,&stxs[1],&dnxs[1],&stys[1],&dnys[1],mpi_rank,mpi_numx,mpi_numy);
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-1;j++){
#pragma simd
	for (i=1;i<nx-1;i++){
	  ss=nx*j+i;
	  dnvv[0*nxy+ss]=calc_dnvv(&ro[ss],&vp[0*nxy+ss],&vp[1*nxy+ss],&vp[2*nxy+ss],&ve[0*nxy+ss],&ve[1*nxy+ss],&ve[2*nxy+ss],
				   idx,idy,0,1,nx,0);
	  dnvv[1*nxy+ss]=calc_dnvv(&ro[ss],&vp[1*nxy+ss],&vp[2*nxy+ss],&vp[0*nxy+ss],&ve[1*nxy+ss],&ve[2*nxy+ss],&ve[0*nxy+ss],
				   idy,0,idx,nx,0,1);
	  dnvv[2*nxy+ss]=calc_dnvv(&ro[ss],&vp[2*nxy+ss],&vp[0*nxy+ss],&vp[1*nxy+ss],&ve[2*nxy+ss],&ve[0*nxy+ss],&ve[1*nxy+ss],
				   0,idx,idy,0,1,nx);
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[]={0,0,0},sty[]={0,0,0};
	int dnx[]={+dnxs[6],-dnxs[6],-dnxs[5]};
	int dny[]={-dnys[6],+dnys[6],-dnys[4]};
	double *pv[]={&dnvv[0*nxy],&dnvv[1*nxy],&dnvv[2*nxy]};
	boundary(pv,3,nx,ny,xoff,yoff,stx,dnx,sty,dny,mpi_rank,mpi_numx,mpi_numy);
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

	  /* Linear interpolation of numerical flux of By */
	  mhd_lr_fb(&vx[ss],&vy[ss],&bx[ss],&cy[ss],1,lfun_lr,vl,vr);
	  ql[sl]=vl[0];		/* by*vx-bx*vy @ i+1/2 Left */
	  qr[sr]=vr[0];		/* by*vx-bx*vy @ i-1/2 Right */
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
	  /* Split central and upwind parts in numerical flux of By */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fx[nm*ss+5]-=fc[ss];	      /* Upwind part */
	}
      }

      /* Modification for extended MHD */
      {
	int stx[]={1,1};
	int dnx[]={-dnxs[5],-dnxs[6]};
	int sty[]={0,0};
	int dny[]={+dnys[5],+dnys[6]};
#ifdef _OPENMP
#pragma omp for
#endif
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=3;i<nx-2;i++){
	    ss=nx*j+i;
	    rtmp[ss]=0.5*(ro[ss- 1]+ro[ss]);
	    enew[0*nxy+ss]=fx[nm*ss+5]+fc[ss]; /* -ez */
	    enew[1*nxy+ss]=fx[nm*ss+6];	       /* +ey */
	    double fac=mpmt*smemt*de/rtmp[ss];
	    enew[0*nxy+ss]+=-fac*0.5*(dnvv[2*nxy+ss -1]+dnvv[2*nxy+ss]); /* -ez */
	    enew[1*nxy+ss]+=+fac*0.5*(dnvv[1*nxy+ss -1]+dnvv[1*nxy+ss]); /* +ey */
	  }
	}
#ifdef _OPENMP
#pragma omp single
#endif
	{
	  emhd_eorg2enew(&enew[0*nxy],&enew[0*nxy],rtmp,smpmt*de,dx,dy,
			 nx,ny,xoff,yoff,stx[0],dnx[0],sty[0],dny[0],mpi_rank,mpi_numx,mpi_numy);
	  emhd_eorg2enew(&enew[1*nxy],&enew[1*nxy],rtmp,smpmt*de,dx,dy,
			 nx,ny,xoff,yoff,stx[1],dnx[1],sty[1],dny[1],mpi_rank,mpi_numx,mpi_numy);
	}
#ifdef _OPENMP
#pragma omp for
#endif
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=3;i<nx-2;i++){
	    ss=nx*j+i;
	    double dey=+(enew[1*nxy+ss]-(fx[nm*ss+6]));
	    double dez=-(enew[0*nxy+ss]-(fx[nm*ss+5]+fc[ss]));
	    double bytmp=0.5*(cy[ss- 1]+cy[ss]); /* For simplicity. (func_bc can be used) */
	    double bztmp=0.5*(cz[ss- 1]+cz[ss]);
	    fc[ss]     =enew[0*nxy+ss]-fx[nm*ss+5]; /* -ez, central */
	    fx[nm*ss+6]=enew[1*nxy+ss];		    /* +ey */
	    fx[nm*ss+7]+=dey*bztmp-dez*bytmp;
	  }
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

	  /* Linear interpolation of numerical flux of Bx */
	  mhd_lr_fb(&vy[ss],&vx[ss],&by[ss],&cx[ss],nx,lfun_lr,vl,vr);
	  ql[sl]=vl[0];		/* bx*vy-by*vx @ j+1/2 Left */
	  qr[sr]=vr[0];		/* bx*vy-by*vx @ j-1/2 Right */
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
	  /* Split central and upwind parts in numerical flux of Bx */
	  fc[ss]=0.5*(ql[ss]+qr[ss]); /* Central part */
	  fy[nm*ss+4]-=fc[ss];	      /* Upwind part */
	}
      }

      /* Modification for extended MHD */
      {
	int stx[]={0,0};
	int dnx[]={+dnxs[6],+dnxs[4]};
	int sty[]={1,1};
	int dny[]={-dnys[6],-dnys[4]};
#ifdef _OPENMP
#pragma omp for
#endif
	for (j=3;j<ny-2;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    ss=nx*j+i;
	    rtmp[ss]=0.5*(ro[ss-nx]+ro[ss]);
	    enew[0*nxy+ss]=fy[nm*ss+6];	       /* -ex */
	    enew[1*nxy+ss]=fy[nm*ss+4]+fc[ss]; /* +ez */
	    double fac=mpmt*smemt*de/rtmp[ss];
	    enew[0*nxy+ss]+=-fac*0.5*(dnvv[0*nxy+ss-nx]+dnvv[0*nxy+ss]); /* -ex */
	    enew[1*nxy+ss]+=+fac*0.5*(dnvv[2*nxy+ss-nx]+dnvv[2*nxy+ss]); /* +ez */
	  }
	}
#ifdef _OPENMP
#pragma omp single
#endif
	{
	  emhd_eorg2enew(&enew[0*nxy],&enew[0*nxy],rtmp,smpmt*de,dx,dy,
			 nx,ny,xoff,yoff,stx[0],dnx[0],sty[0],dny[0],mpi_rank,mpi_numx,mpi_numy);
	  emhd_eorg2enew(&enew[1*nxy],&enew[1*nxy],rtmp,smpmt*de,dx,dy,
			 nx,ny,xoff,yoff,stx[1],dnx[1],sty[1],dny[1],mpi_rank,mpi_numx,mpi_numy);
	}
#ifdef _OPENMP
#pragma omp for
#endif
	for (j=3;j<ny-2;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    ss=nx*j+i;
	    double dez=+(enew[1*nxy+ss]-(fy[nm*ss+4]+fc[ss]));
	    double dex=-(enew[0*nxy+ss]-(fy[nm*ss+6]));
	    double bztmp=0.5*(cz[ss-nx]+cz[ss]); /* For simplicity. (func_bc can be used) */
	    double bxtmp=0.5*(cx[ss-nx]+cx[ss]);
	    fc[ss]     =enew[1*nxy+ss]-fy[nm*ss+4]; /* +ez, central */
	    fy[nm*ss+6]=enew[0*nxy+ss];		    /* -ex */
	    fy[nm*ss+7]+=dez*bxtmp-dex*bztmp;
	  }
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
    boundary(p,nm,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);
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
  /* Modification for extended MHD */
  free(rtmp);
  free(enew);
  free(vp);
  free(ve);
  free(dnvv);
}

void mhd_diff2d(double *p[], double dt, double dx, double dy,
		double *nu, double *eta,
		int nm, int nx, int ny, int xoff, int yoff, double gamma,
		int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D visco-resistive code */
/* 2nd/4th accuracy in space, CT grid spacing */
/* Need four offsets for boundary */

/* p[0..7] = ro,mx,my,mz,bx,by,bz,en */
/* nu = kinetic viscosity coefficient @ cell center (i,j) */
/* eta = resistivity coefficient @ cell center (i,j) */
/* nm = Number of variables (8 in multi-D) */
/* gamma: specific heat ratio */
{
  int i,j,ss,rk;
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double idx=1.0/dx,idy=1.0/dy;
  const double dtdx=dt*idx,dtdy=dt*idy;
  const double cflg[]={1,1,1,1,0,0,1,1}; /* Flag for cell center variable */
  double (*func_fc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ro=p[0],*mx=p[1],*my=p[2],*mz=p[3],*bx=p[4],*by=p[5],*bz=p[6],*en=p[7];
  double *ut,*fx,*fy;
  double *vx,*vy,*vz,*ex,*ey,*ez;

  ut=(double*)malloc(sizeof(double)*nm*nxy);
  fx=(double*)malloc(sizeof(double)*nm*nxy);
  fy=(double*)malloc(sizeof(double)*nm*nxy);
  vx=(double*)malloc(sizeof(double)*nxy);
  vy=(double*)malloc(sizeof(double)*nxy);
  vz=(double*)malloc(sizeof(double)*nxy);
  ex=(double*)malloc(sizeof(double)*nxy);
  ey=(double*)malloc(sizeof(double)*nxy);
  ez=(double*)malloc(sizeof(double)*nxy);

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

      /* Velocity at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxy;ss++){
	double iro=1.0/ro[ss];
	vx[ss]=mx[ss]*iro;
	vy[ss]=my[ss]*iro;
	vz[ss]=mz[ss]*iro;
      }
      
      /* Resistive E-field at CT grid */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-1;j++){
#pragma simd
	for (i=2;i<nx-1;i++){
	  ss=nx*j+i;
	  mhd_ct_eres(&bx[ss],&by[ss],&bz[ss],&eta[ss],idx,idy,0.0,1,nx,0,func_df,&ex[ss],&ey[ss],&ez[ss]);
	}
      }

      /* Numerical fluxes at cell faces */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=2;j<ny-1;j++){
#pragma simd
	for (i=2;i<nx-1;i++){
	  ss=nx*j+i;
	  double rnu,flux[2][8]={{0.0},{0.0}};
	  
	  /* i-1/2, X direction */
	  rnu=0.5*(ro[ss- 1]*nu[ss- 1]+ro[ss]*nu[ss]);
	  calc_flux_viscous(&vx[ss],&vy[ss],&vz[ss],rnu,idx,idy,0.0,1,nx,0,func_fc,func_df,&flux[0][1],&flux[0][2],&flux[0][3],&flux[0][7]);
	  calc_flux_resistive(&by[ss],&bz[ss],&ey[ss],&ez[ss],1,nx,0,func_fc,&flux[1][5],&flux[1][6],&flux[1][7]);
	  fx[nm*ss+0]=flux[0][0]+flux[1][0];	/* ro */
	  fx[nm*ss+1]=flux[0][1]+flux[1][1];	/* mx */
	  fx[nm*ss+2]=flux[0][2]+flux[1][2];	/* my */
	  fx[nm*ss+3]=flux[0][3]+flux[1][3];	/* mz */
	  fx[nm*ss+4]=0;	/* bx */
	  fx[nm*ss+5]=flux[0][5]+flux[1][5];	/* by */
	  fx[nm*ss+6]=flux[0][6]+flux[1][6];	/* bz */
	  fx[nm*ss+7]=flux[0][7]+flux[1][7];	/* en */

	  /* j-1/2, Y direction */
	  rnu=0.5*(ro[ss-nx]*nu[ss-nx]+ro[ss]*nu[ss]);
	  calc_flux_viscous(&vy[ss],&vz[ss],&vx[ss],rnu,idy,0.0,idx,nx,0,1,func_fc,func_df,&flux[0][1],&flux[0][2],&flux[0][3],&flux[0][7]);
	  calc_flux_resistive(&bz[ss],&bx[ss],&ez[ss],&ex[ss],nx,0,1,func_fc,&flux[1][5],&flux[1][6],&flux[1][7]);
	  fy[nm*ss+0]=flux[0][0]+flux[1][0];	/* ro */
	  fy[nm*ss+2]=flux[0][1]+flux[1][1];	/* my */
	  fy[nm*ss+3]=flux[0][2]+flux[1][2];	/* mz */
	  fy[nm*ss+1]=flux[0][3]+flux[1][3];	/* mx */
	  fy[nm*ss+5]=0;	/* by */
	  fy[nm*ss+6]=flux[0][5]+flux[1][5];	/* bz */
	  fy[nm*ss+4]=flux[0][6]+flux[1][6];	/* bx */
	  fy[nm*ss+7]=flux[0][7]+flux[1][7];	/* en */
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
	  double dummy;
	  double *val[]={&ro[ss],&mx[ss],&my[ss],&mz[ss],&dummy,&dummy,&bz[ss],&en[ss]};
	  double val0[]={ut[0*nxy+ss],ut[1*nxy+ss],ut[2*nxy+ss],ut[3*nxy+ss],ut[4*nxy+ss],ut[5*nxy+ss],ut[6*nxy+ss],ut[7*nxy+ss]};
	  mhd_updt2d(val,val0,&fx[nm*ss+0],&fy[nm*ss+0],dtdx,dtdy,rk_fac[rk],cflg,1,nx,func_df);
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
    boundary(p,nm,nx,ny,xoff,yoff,stxs,dnxs,stys,dnys,mpi_rank,mpi_numx,mpi_numy);
  }
  
  free(ut);
  free(fx);
  free(fy);
  free(vx);
  free(vy);
  free(vz);
  free(ex);
  free(ey);
  free(ez);
}