#include "mhd2d_class.hpp"

void MHD2D::ideal(double dt)
{
  // 2D ideal MHD simulation
  int i,j,ss,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  static const int nxy=nx*ny;
  const double dtdx=dt/dx,dtdy=dt/dy;
  void (*func_flux)(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, const double*,
		    double*, double*, double*, double*, double*, double*, double*)=riemann[RMN];
  void (*lfun_lr)(const double *f, double *fl, double *fr)=l_interp[ODR-1];
  void (*func_lr)(const double *f, double *fl, double *fr)=interpol[ODR-1];
  double (*func_bc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ut,*ul,*ur,*fx,*fy,*ql,*qr,*fc;
  double *ez,*ct,*dvx,*dvy;

  ut=new double[nm*nxy];
  ul=new double[nm*nxy];
  ur=new double[nm*nxy];
  fx=new double[nm*nxy];
  fy=new double[nm*nxy];
  ql=new double[nxy];
  qr=new double[nxy];
  fc=new double[nxy];
  ez=new double[nxy];
  ct=new double[nxy];
  dvx=new double[nxy];
  dvy=new double[nxy];

  // Bz @ cell center
  cz=bz;
  
  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nxy],val[i],nxy);
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
	  ss=nx*j+i;
	  cx[ss]=bcell(&bx[ss], 1,func_bc);
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-2;j++){
#pragma simd
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  cy[ss]=bcell(&by[ss],nx,func_bc);
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[2]={0,0},sty[2]={0,0};
	double *pc[]={cx,cy};
	bound(pc,2,stx,&dnxs[4],sty,&dnys[4]);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxy;ss++){
	prmtv(ss);
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
	bound(pv,2,stx,dnx,sty,dny);
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
		      cx[ss],gam,1,func_lr,vl,vr);
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
		    bn,gam,dvsd,
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
		      cy[ss],gam,nx,func_lr,vl,vr);
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
		    bn,gam,dvsd,
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
	  double *val1[]={val[0]+ss,val[1]+ss,val[2]+ss,val[3]+ss,val[4]+ss,val[5]+ss,val[6]+ss,val[7]+ss};
	  double val0[]={ut[0*nxy+ss],ut[1*nxy+ss],ut[2*nxy+ss],ut[3*nxy+ss],ut[4*nxy+ss],ut[5*nxy+ss],ut[6*nxy+ss],ut[7*nxy+ss]};
	  mhd_updt2d(val1[0],val0[0],&fx[nm*ss+0],&fy[nm*ss+0],dtdx,dtdy,rk_fac[rk],nm,nm*nx,func_df); // ro
	  mhd_updt2d(val1[1],val0[1],&fx[nm*ss+1],&fy[nm*ss+1],dtdx,dtdy,rk_fac[rk],nm,nm*nx,func_df); // mx
	  mhd_updt2d(val1[2],val0[2],&fx[nm*ss+2],&fy[nm*ss+2],dtdx,dtdy,rk_fac[rk],nm,nm*nx,func_df); // my
	  mhd_updt2d(val1[3],val0[3],&fx[nm*ss+3],&fy[nm*ss+3],dtdx,dtdy,rk_fac[rk],nm,nm*nx,func_df); // mz
	  mhd_updt2d(val1[6],val0[6],&fx[nm*ss+6],&fy[nm*ss+6],dtdx,dtdy,rk_fac[rk],nm,nm*nx,func_df); // bz
	  mhd_updt2d(val1[7],val0[7],&fx[nm*ss+7],&fy[nm*ss+7],dtdx,dtdy,rk_fac[rk],nm,nm*nx,func_df); // en
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
    bound(val,nm,stxs,dnxs,stys,dnys);
  }
  
  delete[] ut;
  delete[] ul;
  delete[] ur;
  delete[] fx;
  delete[] fy;
  delete[] ql;
  delete[] qr;
  delete[] fc;
  delete[] ez;
  delete[] ct;
  delete[] dvx;
  delete[] dvy;
}

