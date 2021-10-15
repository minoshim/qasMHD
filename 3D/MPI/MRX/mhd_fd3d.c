#include "mhd_fd3d.h"
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "mhd_common.h"
#include "boundary.h"

/* [IMPORTANT] Macros are defined in mhd_fd3d.h to select Riemann solver, spatial and temporal accuracy, CUCT flag. */

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
extern int dnzs[8];
/* Staggered grid flag, defined in global.hpp */
extern int stxs[8];
extern int stys[8];
extern int stzs[8];

void mhd_fd3d(double *p[], double dt, double dx, double dy, double dz,
	      int nm, int nx, int ny, int nz, int xoff, int yoff, int zoff, double gamma,
	      int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D finite-difference code for MHD */
/* 2nd, 3rd, and 4th accuracy in space */
/* ROE/HLLD/LHLLD/MLAU + Central-Upwind-CT */
/* Need four offsets for boundary */

/* p[0..7] = ro,mx,my,mz,bx,by,bz,en */
/* nm = Number of variables (8 in multi-D) */
/* gamma: specific heat ratio */
{
  int i,j,k,ss,s1,rk;
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny,nxyz=nx*ny*nz;
  const double dtdx=dt/dx,dtdy=dt/dy,dtdz=dt/dz;
  const double cflg[]={1,1,1,1,0,0,0,1}; /* Flag for cell center variable */
  void (*func_flux)(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, const double*,
		    double*, double*, double*, double*, double*, double*, double*)=riemann[RMN];
  void (*lfun_lr)(const double *f, double *fl, double *fr)=l_interp[ODR-1];
  void (*func_lr)(const double *f, double *fl, double *fr)=interpol[ODR-1];
  double (*func_bc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ro=p[0],*mx=p[1],*my=p[2],*mz=p[3],*bx=p[4],*by=p[5],*bz=p[6],*en=p[7];
  double *ut,*fx,*fy,*fz,*fc;
  double *vx,*vy,*vz,*pr,*cx,*cy,*cz,*ex,*ey,*ez;
  double *ct,*dvx,*dvy,*dvz;

  ut=(double*)malloc(sizeof(double)*nm*nxyz);
  fx=(double*)malloc(sizeof(double)*nm*nxyz);
  fy=(double*)malloc(sizeof(double)*nm*nxyz);
  fz=(double*)malloc(sizeof(double)*nm*nxyz);
  fc=(double*)malloc(sizeof(double)*2*nxyz);
  vx=(double*)malloc(sizeof(double)*nxyz);
  vy=(double*)malloc(sizeof(double)*nxyz);
  vz=(double*)malloc(sizeof(double)*nxyz);
  pr=(double*)malloc(sizeof(double)*nxyz);
  cx=(double*)malloc(sizeof(double)*nxyz);
  cy=(double*)malloc(sizeof(double)*nxyz);
  cz=(double*)malloc(sizeof(double)*nxyz);
  ex=(double*)malloc(sizeof(double)*nxyz);
  ey=(double*)malloc(sizeof(double)*nxyz);
  ez=(double*)malloc(sizeof(double)*nxyz);
  ct=(double*)malloc(sizeof(double)*3*nxyz);
  dvx=(double*)malloc(sizeof(double)*nxyz);
  dvy=(double*)malloc(sizeof(double)*nxyz);
  dvz=(double*)malloc(sizeof(double)*nxyz);

  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nxyz],p[i],nxyz);
  }

  /* Runge-Kutta stage */
  for (rk=0;rk<R_K;rk++){
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ss,s1)
#endif
    {

      /* B at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=1;i<nx-2;i++){
	    double val[]={bx[nx*(ny*k+j)+(i-1)],bx[nx*(ny*k+j)+(i+0)],
			  bx[nx*(ny*k+j)+(i+1)],bx[nx*(ny*k+j)+(i+2)]};
	    cx[nx*(ny*k+j)+i]=func_bc(&val[1]);
	  }
	}
	for (j=1;j<ny-2;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    double val[]={by[nx*(ny*k+(j-1))+i],by[nx*(ny*k+(j+0))+i],
			  by[nx*(ny*k+(j+1))+i],by[nx*(ny*k+(j+2))+i]};
	    cy[nx*(ny*k+j)+i]=func_bc(&val[1]);
	  }
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=1;k<nz-2;k++){
	for (j=0;j<ny;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    double val[]={bz[nx*(ny*(k-1)+j)+i],bz[nx*(ny*(k+0)+j)+i],
			  bz[nx*(ny*(k+1)+j)+i],bz[nx*(ny*(k+2)+j)+i]};
	    cz[nx*(ny*k+j)+i]=func_bc(&val[1]);
	  }
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[3]={0,0,0},sty[3]={0,0,0},stz[3]={0,0,0};
	double *pc[]={cx,cy,cz};
	boundary(pc,3,nx,ny,nz,xoff,yoff,zoff,stx,&dnxs[4],sty,&dnys[4],stz,&dnzs[4],mpi_rank,mpi_numx,mpi_numy,mpi_numz);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxyz;ss++){
	mhd_prmtv(ro[ss],mx[ss],my[ss],mz[ss],cx[ss],cy[ss],cz[ss],en[ss],
		  &vx[ss],&vy[ss],&vz[ss],&pr[ss],gamma);
	ex[ss]=ey[ss]=ez[ss]=0.0;		/* Necessary initialize at cell corner */
	ct[nxyz*0+ss]=ct[nxyz*1+ss]=ct[nxyz*2+ss]=0.5;	
      }

      /* CUCT upwind weighting */
#if (CTW)
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=0;k<nz;k++){
	for (j=1;j<ny;j++){
#pragma simd
	  /* X-Y plane for Ez */
	  for (i=1;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    ct[nxyz*0+ss]=mhd_cuct_weight(&ro[ss],&vx[ss],&vy[ss],&bx[ss],&by[ss],1,nx);
	  }
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=1;k<nz;k++){
	for (j=1;j<ny;j++){
#pragma simd
	  /* Y-Z plane for Ex */
	  for (i=0;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    ct[nxyz*1+ss]=mhd_cuct_weight(&ro[ss],&vy[ss],&vz[ss],&by[ss],&bz[ss],nx,nxy);
	  }
	}
	for (j=0;j<ny;j++){
#pragma simd
	  /* Z-X plane for Ey */
	  for (i=1;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    ct[nxyz*2+ss]=mhd_cuct_weight(&ro[ss],&vz[ss],&vx[ss],&bz[ss],&bx[ss],nxy,1);
	  }
	}
      }
#endif	/* CTW */

      /* dvx,dvy,dvz at cell center for shock detection */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=1;k<nz-1;k++){
	for (j=1;j<ny-1;j++){
#pragma simd
	  for (i=1;i<nx-1;i++){
	    ss=nx*(ny*k+j)+i;
	    dvx[ss]=min((vx[nx*(ny*k+j)+(i+1)]-vx[nx*(ny*k+j)+i]),
			(vx[nx*(ny*k+j)+i]-vx[nx*(ny*k+j)+(i-1)]));
	    dvy[ss]=min((vy[nx*(ny*k+(j+1))+i]-vy[nx*(ny*k+j)+i]),
			(vy[nx*(ny*k+j)+i]-vy[nx*(ny*k+(j-1))+i]));
	    dvz[ss]=min((vz[nx*(ny*(k+1)+j)+i]-vz[nx*(ny*k+j)+i]),
			(vz[nx*(ny*k+j)+i]-vz[nx*(ny*(k-1)+j)+i]));
	  }
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[3]={0,0,0},sty[3]={0,0,0},stz[3]={0,0,0};
	int dnx[3]={-dnxs[1],dnxs[2],dnxs[3]},dny[3]={dnys[1],-dnys[2],dnys[3]},dnz[3]={dnzs[1],dnzs[2],-dnzs[3]};
	double *pv[]={dvx,dvy,dvz};
	boundary(pv,3,nx,ny,nz,xoff,yoff,zoff,stx,dnx,sty,dny,stz,dnz,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
      }

      /* Numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur,*ql,*qr;
	ul=(double*)malloc(sizeof(double)*nm*nx);
	ur=(double*)malloc(sizeof(double)*nm*nx);
	ql=(double*)malloc(sizeof(double)*2*nx);
	qr=(double*)malloc(sizeof(double)*2*nx);
	for (j=0;j<ny;j++){
	  /* Interpolation */
	  for (i=2;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1l=i+1;
	    s1r=i;
	    sl=nx*(ny*k+j)+(i+1);
	    sr=nx*(ny*k+j)+i;
	    double vl[7],vr[7];
	    mhd_lrstate(&ro[ss],&vx[ss],&vy[ss],&vz[ss],&cy[ss],&cz[ss],&pr[ss],
			cx[ss],gamma,1,func_lr,vl,vr);
	    /* Left-face @ i+1/2 */
	    ul[nm*s1l+0]=vl[0];	 /* ro */
	    ul[nm*s1l+1]=vl[1];	 /* vx */
	    ul[nm*s1l+2]=vl[2];	 /* vy */
	    ul[nm*s1l+3]=vl[3];	 /* vz */
	    ul[nm*s1l+4]=bx[sl]; /* bx */
	    ul[nm*s1l+5]=vl[4];	 /* by */
	    ul[nm*s1l+6]=vl[5];	 /* bz */
	    ul[nm*s1l+7]=vl[6];	 /* pr */
	    /* Right-face @ i-1/2 */
	    ur[nm*s1r+0]=vr[0];	 /* ro */
	    ur[nm*s1r+1]=vr[1];	 /* vx */
	    ur[nm*s1r+2]=vr[2];	 /* vy */
	    ur[nm*s1r+3]=vr[3];	 /* vz */
	    ur[nm*s1r+4]=bx[sr]; /* bx */
	    ur[nm*s1r+5]=vr[4];	 /* by */
	    ur[nm*s1r+6]=vr[5];	 /* bz */
	    ur[nm*s1r+7]=vr[6];	 /* pr */

	    /* Linear interpolation of numerical flux of By */
	    mhd_lr_fb(&vx[ss],&vy[ss],&bx[ss],&cy[ss],1,lfun_lr,vl,vr);
	    ql[2*s1l+0]=vl[0];	/* by*vx-bx*vy @ i+1/2 Left */
	    qr[2*s1r+0]=vr[0];	/* by*vx-bx*vy @ i-1/2 Right */
	    /* Linear interpolation of numerical flux of Bz */
	    mhd_lr_fb(&vx[ss],&vz[ss],&bx[ss],&cz[ss],1,lfun_lr,vl,vr);
	    ql[2*s1l+1]=vl[0];	/* bz*vx-bx*vz @ i+1/2 Left */
	    qr[2*s1r+1]=vr[0];	/* bz*vx-bx*vz @ i-1/2 Right */
	  }
	  /* Numerical flux */
#pragma simd
	  for (i=3;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1=i;
	    double flux[8]={0},bn=0.5*(ul[nm*s1+4]+ur[nm*s1+4]);
	    double dvtm=min(min(dvy[nx*(ny*k+j)+(i-1)],dvy[nx*(ny*k+j)+i]),
			    min(dvz[nx*(ny*k+j)+(i-1)],dvz[nx*(ny*k+j)+i]));
	    double dvsd[2]={(vx[nx*(ny*k+j)+i]-vx[nx*(ny*k+j)+(i-1)]),dvtm};
	    func_flux(ul[nm*s1+0],ul[nm*s1+1],ul[nm*s1+2],ul[nm*s1+3],ul[nm*s1+5],ul[nm*s1+6],ul[nm*s1+7],
		      ur[nm*s1+0],ur[nm*s1+1],ur[nm*s1+2],ur[nm*s1+3],ur[nm*s1+5],ur[nm*s1+6],ur[nm*s1+7],
		      bn,gamma,dvsd,
		      &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	    fx[nm*ss+0]=flux[0]; /* ro */
	    fx[nm*ss+1]=flux[1]; /* mx */
	    fx[nm*ss+2]=flux[2]; /* my */
	    fx[nm*ss+3]=flux[3]; /* mz */
	    fx[nm*ss+4]=0;	 /* bx */
	    fx[nm*ss+5]=flux[5]; /* by */
	    fx[nm*ss+6]=flux[6]; /* bz */
	    fx[nm*ss+7]=flux[7]; /* en */

	    /* Split central and upwind parts in numerical flux of By */
	    fc[2*ss+0]=0.5*(ql[2*s1+0]+qr[2*s1+0]); /* Central part */
	    fx[nm*ss+5]-=fc[2*ss+0];		     /* Upwind part */
	    /* Split central and upwind parts in numerical flux of Bz */
	    fc[2*ss+1]=0.5*(ql[2*s1+1]+qr[2*s1+1]); /* Central part */
	    fx[nm*ss+6]-=fc[2*ss+1];		     /* Upwind part */
	  }
	}
	free(ul);
	free(ur);
	free(ql);
	free(qr);
      }
      /* Numerical flux of By at cell corner along Y to calculate Ez */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*ny);
	ur=(double*)malloc(sizeof(double)*2*ny);
	for (i=3;i<nx-2;i++){
	  /* Interpolation */
#pragma simd
	  for (j=2;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1l=j+1;
	    s1r=j;
	    sl=nx*(ny*k+(j+1))+i;
	    sr=nx*(ny*k+j)+i;
	    double vl,vr;
	    mhd_lr_single(&fx[nm*ss+5],nm*nx,lfun_lr,&vl,&vr); /* Upwind part */
	    ul[2*s1l+0]=vl;
	    ur[2*s1r+0]=vr;
	    mhd_lr_single(&fc[2*ss+0],2*nx,lfun_lr,&vl,&vr); /* Central part */
	    ul[2*s1l+1]=vl;
	    ur[2*s1r+1]=vr;
	  }
	  /* CT method */
#pragma simd
	  for (j=3;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1=j;
	    ez[ss]+=-0.5*((ul[2*s1+0]+ur[2*s1+0])+(1.0-ct[nxyz*0+ss])*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      /* Numerical flux of Bz at cell corner along Z to calculate Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nz);
	ur=(double*)malloc(sizeof(double)*2*nz);
	for (i=3;i<nx-2;i++){
	  /* Interpolation */
#pragma simd
	  for (k=2;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1l=k+1;
	    s1r=k;
	    sl=nx*(ny*(k+1)+j)+i;
	    sr=nx*(ny*k+j)+i;
	    double vl,vr;
	    mhd_lr_single(&fx[nm*ss+6],nm*nxy,lfun_lr,&vl,&vr); /* Upwind part */
	    ul[2*s1l+0]=vl;
	    ur[2*s1r+0]=vr;
	    mhd_lr_single(&fc[2*ss+1],2*nxy,lfun_lr,&vl,&vr); /* Central part */
	    ul[2*s1l+1]=vl;
	    ur[2*s1r+1]=vr;
	  }
	  /* CT method */
#pragma simd
	  for (k=3;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1=k;
	    ey[ss]+=+0.5*((ul[2*s1+0]+ur[2*s1+0])+ct[nxyz*2+ss]*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }

      /* Numerical flux at cell face along Y */	  
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur,*ql,*qr;
	ul=(double*)malloc(sizeof(double)*nm*ny);
	ur=(double*)malloc(sizeof(double)*nm*ny);
	ql=(double*)malloc(sizeof(double)*2*ny);
	qr=(double*)malloc(sizeof(double)*2*ny);
	for (i=0;i<nx;i++){
	  /* Interpolation */
	  for (j=2;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1l=j+1;
	    s1r=j;
	    sl=nx*(ny*k+(j+1))+i;
	    sr=nx*(ny*k+j)+i;
	    double vl[7],vr[7];
	    mhd_lrstate(&ro[ss],&vy[ss],&vz[ss],&vx[ss],&cz[ss],&cx[ss],&pr[ss],
			cy[ss],gamma,nx,func_lr,vl,vr);
	    /* Left-face @ j+1/2 */
	    ul[nm*s1l+0]=vl[0];  /* ro */
	    ul[nm*s1l+1]=vl[1];	 /* vy */
	    ul[nm*s1l+2]=vl[2];	 /* vz */
	    ul[nm*s1l+3]=vl[3];	 /* vx */
	    ul[nm*s1l+4]=by[sl]; /* by */
	    ul[nm*s1l+5]=vl[4];	 /* bz */
	    ul[nm*s1l+6]=vl[5];	 /* bx */
	    ul[nm*s1l+7]=vl[6];	 /* pr */
	    /* Right-face @ j-1/2 */
	    ur[nm*s1r+0]=vr[0];	 /* ro */
	    ur[nm*s1r+1]=vr[1];	 /* vy */
	    ur[nm*s1r+2]=vr[2];	 /* vz */
	    ur[nm*s1r+3]=vr[3];	 /* vx */
	    ur[nm*s1r+4]=by[sr]; /* by */
	    ur[nm*s1r+5]=vr[4];	 /* bz */
	    ur[nm*s1r+6]=vr[5];	 /* bx */
	    ur[nm*s1r+7]=vr[6];	 /* pr */

	    /* Linear interpolation of numerical flux of Bz */
	    mhd_lr_fb(&vy[ss],&vz[ss],&by[ss],&cz[ss],nx,lfun_lr,vl,vr);
	    ql[2*s1l+0]=vl[0];	/* bz*vy-by*vz @ j+1/2 Left */
	    qr[2*s1r+0]=vr[0];	/* bz*vy-by*vz @ j-1/2 Right */
	    /* Linear interpolation of numerical flux of Bx */
	    mhd_lr_fb(&vy[ss],&vx[ss],&by[ss],&cx[ss],nx,lfun_lr,vl,vr);
	    ql[2*s1l+1]=vl[0];	/* bx*vy-by*vx @ j+1/2 Left */
	    qr[2*s1r+1]=vr[0];	/* bx*vy-by*vx @ j-1/2 Right */
	  }
	  /* Numerical flux */
#pragma simd
	  for (j=3;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1=j;
	    double flux[8]={0},bn=0.5*(ul[nm*s1+4]+ur[nm*s1+4]);
	    double dvtm=min(min(dvz[nx*(ny*k+(j-1))+i],dvz[nx*(ny*k+j)+i]),
			    min(dvx[nx*(ny*k+(j-1))+i],dvx[nx*(ny*k+j)+i]));
	    double dvsd[2]={(vy[nx*(ny*k+j)+i]-vy[nx*(ny*k+(j-1))+i]),dvtm};
	    func_flux(ul[nm*s1+0],ul[nm*s1+1],ul[nm*s1+2],ul[nm*s1+3],ul[nm*s1+5],ul[nm*s1+6],ul[nm*s1+7],
		      ur[nm*s1+0],ur[nm*s1+1],ur[nm*s1+2],ur[nm*s1+3],ur[nm*s1+5],ur[nm*s1+6],ur[nm*s1+7],
		      bn,gamma,dvsd,
		      &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	    fy[nm*ss+0]=flux[0]; /* ro */
	    fy[nm*ss+2]=flux[1]; /* my */
	    fy[nm*ss+3]=flux[2]; /* mz */
	    fy[nm*ss+1]=flux[3]; /* mx */
	    fy[nm*ss+5]=0;	 /* by */
	    fy[nm*ss+6]=flux[5]; /* bz */
	    fy[nm*ss+4]=flux[6]; /* bx */
	    fy[nm*ss+7]=flux[7]; /* en */

	    /* Split central and upwind parts in numerical flux of Bz */
	    fc[2*ss+0]=0.5*(ql[2*s1+0]+qr[2*s1+0]); /* Central part */
	    fy[nm*ss+6]-=fc[2*ss+0];		     /* Upwind part */
	    /* Split central and upwind parts in numerical flux of Bx */
	    fc[2*ss+1]=0.5*(ql[2*s1+1]+qr[2*s1+1]); /* Central part */
	    fy[nm*ss+4]-=fc[2*ss+1];		     /* Upwind part */
	  }
	}
	free(ul);
	free(ur);
	free(ql);
	free(qr);
      }
      /* Numerical flux of Bz at cell corner along Z to calculate Ex */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nz);
	ur=(double*)malloc(sizeof(double)*2*nz);
	for (i=0;i<nx;i++){
	  /* Interpolation */
#pragma simd
	  for (k=2;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1l=k+1;
	    s1r=k;
	    sl=nx*(ny*(k+1)+j)+i;
	    sr=nx*(ny*k+j)+i;
	    double vl,vr;
	    mhd_lr_single(&fy[nm*ss+6],nm*nxy,lfun_lr,&vl,&vr); /* Upwind part */
	    ul[2*s1l+0]=vl;
	    ur[2*s1r+0]=vr;
	    mhd_lr_single(&fc[2*ss+0],2*nxy,lfun_lr,&vl,&vr); /* Central part */
	    ul[2*s1l+1]=vl;
	    ur[2*s1r+1]=vr;
	  }
	  /* CT method */
#pragma simd
	  for (k=3;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1=k;
	    ex[ss]+=-0.5*((ul[2*s1+0]+ur[2*s1+0])+(1.0-ct[nxyz*1+ss])*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      /* Numerical flux of Bx at cell corner along X to calculate Ez */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nx);
	ur=(double*)malloc(sizeof(double)*2*nx);
	for (j=3;j<ny-2;j++){
	  /* Interpolation */
#pragma simd
	  for (i=2;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1l=i+1;
	    s1r=i;
	    sl=nx*(ny*k+j)+(i+1);
	    sr=nx*(ny*k+j)+i;
	    double vl,vr;
	    mhd_lr_single(&fy[nm*ss+4],nm,lfun_lr,&vl,&vr); /* Upwind part */
	    ul[2*s1l+0]=vl;
	    ur[2*s1r+0]=vr;
	    mhd_lr_single(&fc[2*ss+1],2,lfun_lr,&vl,&vr); /* Central part */
	    ul[2*s1l+1]=vl;
	    ur[2*s1r+1]=vr;
	  }
	  /* CT method */
#pragma simd
	  for (i=3;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1=i;
	    ez[ss]+=+0.5*((ul[2*s1+0]+ur[2*s1+0])+ct[nxyz*0+ss]*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }

      /* Numerical flux at cell face along Z */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur,*ql,*qr;
	ul=(double*)malloc(sizeof(double)*nm*nz);
	ur=(double*)malloc(sizeof(double)*nm*nz);
	ql=(double*)malloc(sizeof(double)*2*nz);
	qr=(double*)malloc(sizeof(double)*2*nz);
	for (i=0;i<nx;i++){
	  /* Interpolation */
	  for (k=2;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1l=k+1;
	    s1r=k;
	    sl=nx*(ny*(k+1)+j)+i;
	    sr=nx*(ny*k+j)+i;
	    double vl[7],vr[7];
	    mhd_lrstate(&ro[ss],&vz[ss],&vx[ss],&vy[ss],&cx[ss],&cy[ss],&pr[ss],
			cz[ss],gamma,nxy,func_lr,vl,vr);
	    /* Left-face @ k+1/2 */
	    ul[nm*s1l+0]=vl[0];	 /* ro */
	    ul[nm*s1l+1]=vl[1];	 /* vz */
	    ul[nm*s1l+2]=vl[2];	 /* vx */
	    ul[nm*s1l+3]=vl[3];	 /* vy */
	    ul[nm*s1l+4]=bz[sl]; /* bz */
	    ul[nm*s1l+5]=vl[4];	 /* bx */
	    ul[nm*s1l+6]=vl[5];	 /* by */
	    ul[nm*s1l+7]=vl[6];	 /* pr */
	    /* Right-face @ k-1/2 */
	    ur[nm*s1r+0]=vr[0];	 /* ro */
	    ur[nm*s1r+1]=vr[1];	 /* vz */
	    ur[nm*s1r+2]=vr[2];	 /* vx */
	    ur[nm*s1r+3]=vr[3];	 /* vy */
	    ur[nm*s1r+4]=bz[sr]; /* bz */
	    ur[nm*s1r+5]=vr[4];	 /* bx */
	    ur[nm*s1r+6]=vr[5];	 /* by */
	    ur[nm*s1r+7]=vr[6];	 /* pr */

	    /* Linear interpolation of numerical flux of Bx */
	    mhd_lr_fb(&vz[ss],&vx[ss],&bz[ss],&cx[ss],nxy,lfun_lr,vl,vr);
	    ql[2*s1l+0]=vl[0];	/* bx*vz-bz*vx @ k+1/2 Left */
	    qr[2*s1r+0]=vr[0];	/* bx*vz-bz*vx @ k-1/2 Right */
	    /* Linear interpolation of numerical flux of By */
	    mhd_lr_fb(&vz[ss],&vy[ss],&bz[ss],&cy[ss],nxy,lfun_lr,vl,vr);
	    ql[2*s1l+1]=vl[0];	/* by*vz-bz*vy @ k+1/2 Left */
	    qr[2*s1r+1]=vr[0];	/* by*vz-bz*vy @ k-1/2 Right */
	  }
	  /* Numerical flux */
#pragma simd
	  for (k=3;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1=k;
	    double flux[8]={0},bn=0.5*(ul[nm*s1+4]+ur[nm*s1+4]);
	    double dvtm=min(min(dvx[nx*(ny*(k-1)+j)+i],dvx[nx*(ny*k+j)+i]),
			    min(dvy[nx*(ny*(k-1)+j)+i],dvy[nx*(ny*k+j)+i]));
	    double dvsd[2]={(vz[nx*(ny*k+j)+i]-vz[nx*(ny*(k-1)+j)+i]),dvtm};
	    func_flux(ul[nm*s1+0],ul[nm*s1+1],ul[nm*s1+2],ul[nm*s1+3],ul[nm*s1+5],ul[nm*s1+6],ul[nm*s1+7],
		      ur[nm*s1+0],ur[nm*s1+1],ur[nm*s1+2],ur[nm*s1+3],ur[nm*s1+5],ur[nm*s1+6],ur[nm*s1+7],
		      bn,gamma,dvsd,
		      &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);

	    fz[nm*ss+0]=flux[0]; /* ro */
	    fz[nm*ss+3]=flux[1]; /* mz */
	    fz[nm*ss+1]=flux[2]; /* mx */
	    fz[nm*ss+2]=flux[3]; /* my */
	    fz[nm*ss+6]=0;	 /* bz */
	    fz[nm*ss+4]=flux[5]; /* bx */
	    fz[nm*ss+5]=flux[6]; /* by */
	    fz[nm*ss+7]=flux[7]; /* en */

	    /* Split central and upwind parts in numerical flux of Bx */
	    fc[2*ss+0]=0.5*(ql[2*s1+0]+qr[2*s1+0]); /* Central part */
	    fz[nm*ss+4]-=fc[2*ss+0];		     /* Upwind part */
	    /* Split central and upwind parts in numerical flux of By */
	    fc[2*ss+1]=0.5*(ql[2*s1+1]+qr[2*s1+1]); /* Central part */
	    fz[nm*ss+5]-=fc[2*ss+1];		     /* Upwind part */
	  }
	}
	free(ul);
	free(ur);
	free(ql);
	free(qr);
      }
      /* Numerical flux of Bx at cell corner along X to calculate Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=3;k<nz-2;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*nx);
	ur=(double*)malloc(sizeof(double)*2*nx);
	for (j=0;j<ny;j++){
	  /* Interpolation */
#pragma simd
	  for (i=2;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1l=i+1;
	    s1r=i;
	    sl=nx*(ny*k+j)+(i+1);
	    sr=nx*(ny*k+j)+i;
	    double vl,vr;
	    mhd_lr_single(&fz[nm*ss+4],nm,lfun_lr,&vl,&vr); /* Upwind part */
	    ul[2*s1l+0]=vl;
	    ur[2*s1r+0]=vr;
	    mhd_lr_single(&fc[2*ss+0],2,lfun_lr,&vl,&vr); /* Central part */
	    ul[2*s1l+1]=vl;
	    ur[2*s1r+1]=vr;
	  }
	  /* CT method */
#pragma simd
	  for (i=3;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1=i;
	    ey[ss]+=-0.5*((ul[2*s1+0]+ur[2*s1+0])+(1.0-ct[nxyz*2+ss])*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      /* Numerical flux of By at cell corner along Y to calculate Ex */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=3;k<nz-2;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=(double*)malloc(sizeof(double)*2*ny);
	ur=(double*)malloc(sizeof(double)*2*ny);
	for (i=0;i<nx;i++){
	  /* Interpolation */
#pragma simd
	  for (j=2;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1l=j+1;
	    s1r=j;
	    sl=nx*(ny*k+(j+1))+i;
	    sr=nx*(ny*k+j)+i;
	    double vl,vr;
	    mhd_lr_single(&fz[nm*ss+5],nm*nx,lfun_lr,&vl,&vr); /* Upwind part */
	    ul[2*s1l+0]=vl;
	    ur[2*s1r+0]=vr;
	    mhd_lr_single(&fc[2*ss+1],2*nx,lfun_lr,&vl,&vr); /* Central part */
	    ul[2*s1l+1]=vl;
	    ur[2*s1r+1]=vr;
	  }
	  /* CT method */
#pragma simd
	  for (j=3;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1=j;
	    ex[ss]+=+0.5*((ul[2*s1+0]+ur[2*s1+0])+ct[nxyz*1+ss]*(ul[2*s1+1]+ur[2*s1+1]));
	  }
	}
	free(ul);
	free(ur);
      }
      
      /* Update variable at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=zoff;k<nz-zoff;k++){
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    /* Note: cx,cy,cz are unchanged (dummy) since cflg[4..6]=0 */
	    double *val[]={&ro[ss],&mx[ss],&my[ss],&mz[ss],&cx[ss],&cy[ss],&cz[ss],&en[ss]};
	    double val0[]={ut[0*nxyz+ss],ut[1*nxyz+ss],ut[2*nxyz+ss],ut[3*nxyz+ss],ut[4*nxyz+ss],ut[5*nxyz+ss],ut[6*nxyz+ss],ut[7*nxyz+ss]};
	    mhd_updt3d(val,val0,&fx[nm*ss+0],&fy[nm*ss+0],&fz[nm*ss+0],dtdx,dtdy,dtdz,rk_fac[rk],cflg,1,nx,nxy,func_df);
	  }
	}
      }
      /* Update CT Bx and By */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=zoff;k<nz-zoff;k++){
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff+1;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_updt3d_ctb(&bx[ss],ut[4*nxyz+ss],&ey[ss],&ez[ss],dtdy,dtdz,rk_fac[rk],nx,nxy,func_df);
	  }
	}
	for (j=yoff;j<ny-yoff+1;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_updt3d_ctb(&by[ss],ut[5*nxyz+ss],&ez[ss],&ex[ss],dtdz,dtdx,rk_fac[rk],nxy,1,func_df);
	  }
	}
      }
      /* Update CT Bz */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=zoff;k<nz-zoff+1;k++){
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_updt3d_ctb(&bz[ss],ut[6*nxyz+ss],&ex[ss],&ey[ss],dtdx,dtdy,rk_fac[rk],1,nx,func_df);
	  }
	}
      }

    } /* OpenMP */
    
    /* Boundary condition */
    boundary(p,nm,nx,ny,nz,xoff,yoff,zoff,stxs,dnxs,stys,dnys,stzs,dnzs,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  }
  
  free(ut);
  free(fx);
  free(fy);
  free(fz);
  free(fc);
  free(vx);
  free(vy);
  free(vz);
  free(pr);
  free(cx);
  free(cy);
  free(cz);
  free(ex);
  free(ey);
  free(ez);
  free(ct);
  free(dvx);
  free(dvy);
  free(dvz);
}

void mhd_diff3d(double *p[], double dt, double dx, double dy, double dz,
		double *nu, double *eta,
		int nm, int nx, int ny, int nz, int xoff, int yoff, int zoff, double gamma,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D visco-resistive code */
/* 2nd/4th accuracy in space, CT grid spacing */
/* Need four offsets for boundary */

/* p[0..7] = ro,mx,my,mz,bx,by,bz,en */
/* nu = kinetic viscosity coefficient @ cell center (i,j,k) */
/* eta = resistivity coefficient @ cell center (i,j,k) */
/* nm = Number of variables (8 in multi-D) */
/* gamma: specific heat ratio */
{
  int i,j,k,ss,rk;
  const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny,nxyz=nx*ny*nz;
  const double idx=1.0/dx,idy=1.0/dy,idz=1.0/dz;
  const double dtdx=dt*idx,dtdy=dt*idy,dtdz=dt*idz;
  const double cflg[]={1,1,1,1,0,0,0,1}; /* Flag for cell center variable */
  double (*func_fc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ro=p[0],*mx=p[1],*my=p[2],*mz=p[3],*bx=p[4],*by=p[5],*bz=p[6],*en=p[7];
  double *ut,*fx,*fy,*fz;
  double *vx,*vy,*vz,*ex,*ey,*ez;

  ut=(double*)malloc(sizeof(double)*nm*nxyz);
  fx=(double*)malloc(sizeof(double)*nm*nxyz);
  fy=(double*)malloc(sizeof(double)*nm*nxyz);
  fz=(double*)malloc(sizeof(double)*nm*nxyz);
  vx=(double*)malloc(sizeof(double)*nxyz);
  vy=(double*)malloc(sizeof(double)*nxyz);
  vz=(double*)malloc(sizeof(double)*nxyz);
  ex=(double*)malloc(sizeof(double)*nxyz);
  ey=(double*)malloc(sizeof(double)*nxyz);
  ez=(double*)malloc(sizeof(double)*nxyz);

  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nxyz],p[i],nxyz);
  }

  /* Runge-Kutta stage */
  for (rk=0;rk<R_K;rk++){
    
#ifdef _OPENMP
#pragma omp parallel private(i,j,k,ss)
#endif
    {

      /* Velocity at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxyz;ss++){
	double iro=1.0/ro[ss];
	vx[ss]=mx[ss]*iro;
	vy[ss]=my[ss]*iro;
	vz[ss]=mz[ss]*iro;
      }

      /* Resistive E-field at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=2;k<nz-1;k++){
	for (j=2;j<ny-1;j++){
#pragma simd
	  for (i=2;i<nx-1;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_ct_eres(&bx[ss],&by[ss],&bz[ss],&eta[ss],idx,idy,idz,1,nx,nxy,func_df,&ex[ss],&ey[ss],&ez[ss]);
	  }
	}
      }
      
      /* Numerical fluxes at cell faces */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=2;k<nz-1;k++){
	for (j=2;j<ny-1;j++){
#pragma simd
	  for (i=2;i<nx-1;i++){
	    ss=nx*(ny*k+j)+i;
	    double rnu,flux[2][8]={{0.0},{0.0}};

	    /* i-1/2, X direction */
	    rnu=0.5*(ro[ss-  1]*nu[ss-  1]+ro[ss]*nu[ss]);
	    calc_flux_viscous(&vx[ss],&vy[ss],&vz[ss],rnu,idx,idy,idz,1,nx,nxy,func_fc,func_df,&flux[0][1],&flux[0][2],&flux[0][3],&flux[0][7]);
	    calc_flux_resistive(&by[ss],&bz[ss],&ey[ss],&ez[ss],1,nx,nxy,func_fc,&flux[1][5],&flux[1][6],&flux[1][7]);
	    fx[nm*ss+0]=flux[0][0]+flux[1][0];	/* ro */
	    fx[nm*ss+1]=flux[0][1]+flux[1][1];	/* mx */
	    fx[nm*ss+2]=flux[0][2]+flux[1][2];	/* my */
	    fx[nm*ss+3]=flux[0][3]+flux[1][3];	/* mz */
	    fx[nm*ss+4]=0;	/* bx */
	    fx[nm*ss+5]=flux[0][5]+flux[1][5];	/* by */
	    fx[nm*ss+6]=flux[0][6]+flux[1][6];	/* bz */
	    fx[nm*ss+7]=flux[0][7]+flux[1][7];	/* en */

	    /* j-1/2, Y direction */
	    rnu=0.5*(ro[ss- nx]*nu[ss- nx]+ro[ss]*nu[ss]);
	    calc_flux_viscous(&vy[ss],&vz[ss],&vx[ss],rnu,idy,idz,idx,nx,nxy,1,func_fc,func_df,&flux[0][1],&flux[0][2],&flux[0][3],&flux[0][7]);
	    calc_flux_resistive(&bz[ss],&bx[ss],&ez[ss],&ex[ss],nx,nxy,1,func_fc,&flux[1][5],&flux[1][6],&flux[1][7]);
	    fy[nm*ss+0]=flux[0][0]+flux[1][0];	/* ro */
	    fy[nm*ss+2]=flux[0][1]+flux[1][1];	/* my */
	    fy[nm*ss+3]=flux[0][2]+flux[1][2];	/* mz */
	    fy[nm*ss+1]=flux[0][3]+flux[1][3];	/* mx */
	    fy[nm*ss+5]=0;	/* by */
	    fy[nm*ss+6]=flux[0][5]+flux[1][5];	/* bz */
	    fy[nm*ss+4]=flux[0][6]+flux[1][6];	/* bx */
	    fy[nm*ss+7]=flux[0][7]+flux[1][7];	/* en */

	    /* k-1/2, Z direction */
	    rnu=0.5*(ro[ss-nxy]*nu[ss-nxy]+ro[ss]*nu[ss]);
	    calc_flux_viscous(&vz[ss],&vx[ss],&vy[ss],rnu,idz,idx,idy,nxy,1,nx,func_fc,func_df,&flux[0][1],&flux[0][2],&flux[0][3],&flux[0][7]);
	    calc_flux_resistive(&bx[ss],&by[ss],&ex[ss],&ey[ss],nxy,1,nx,func_fc,&flux[1][5],&flux[1][6],&flux[1][7]);
	    fz[nm*ss+0]=flux[0][0]+flux[1][0];	/* ro */
	    fz[nm*ss+3]=flux[0][1]+flux[1][1];	/* mz */
	    fz[nm*ss+1]=flux[0][2]+flux[1][2];	/* mx */
	    fz[nm*ss+2]=flux[0][3]+flux[1][3];	/* my */
	    fz[nm*ss+6]=0;	/* bz */
	    fz[nm*ss+4]=flux[0][5]+flux[1][5];	/* bx */
	    fz[nm*ss+5]=flux[0][6]+flux[1][6];	/* by */
	    fz[nm*ss+7]=flux[0][7]+flux[1][7];	/* en */
	  }
	}
      }

      /* Update variable at cell center */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=zoff;k<nz-zoff;k++){
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    double dummy;
	    double *val[]={&ro[ss],&mx[ss],&my[ss],&mz[ss],&dummy,&dummy,&dummy,&en[ss]};
	    double val0[]={ut[0*nxyz+ss],ut[1*nxyz+ss],ut[2*nxyz+ss],ut[3*nxyz+ss],ut[4*nxyz+ss],ut[5*nxyz+ss],ut[6*nxyz+ss],ut[7*nxyz+ss]};
	    mhd_updt3d(val,val0,&fx[nm*ss+0],&fy[nm*ss+0],&fz[nm*ss+0],dtdx,dtdy,dtdz,rk_fac[rk],cflg,1,nx,nxy,func_df);
	  }
	}
      }
      /* Update CT Bx and By */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (k=zoff;k<nz-zoff;k++){
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff+1;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_updt3d_ctb(&bx[ss],ut[4*nxyz+ss],&ey[ss],&ez[ss],dtdy,dtdz,rk_fac[rk],nx,nxy,func_df);
	  }
	}
	for (j=yoff;j<ny-yoff+1;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_updt3d_ctb(&by[ss],ut[5*nxyz+ss],&ez[ss],&ex[ss],dtdz,dtdx,rk_fac[rk],nxy,1,func_df);
	  }
	}
      }
      /* Update CT Bz */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=zoff;k<nz-zoff+1;k++){
	for (j=yoff;j<ny-yoff;j++){
#pragma simd
	  for (i=xoff;i<nx-xoff;i++){
	    ss=nx*(ny*k+j)+i;
	    mhd_updt3d_ctb(&bz[ss],ut[6*nxyz+ss],&ex[ss],&ey[ss],dtdx,dtdy,rk_fac[rk],1,nx,func_df);
	  }
	}
      }
      
    } /* OpenMP */
    
    /* Boundary condition */
    boundary(p,nm,nx,ny,nz,xoff,yoff,zoff,stxs,dnxs,stys,dnys,stzs,dnzs,mpi_rank,mpi_numx,mpi_numy,mpi_numz);
  }

  free(ut);
  free(fx);
  free(fy);
  free(fz);
  free(vx);
  free(vy);
  free(vz);
  free(ex);
  free(ey);
  free(ez);
}
