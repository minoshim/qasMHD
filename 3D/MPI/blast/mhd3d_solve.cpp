#include "mhd3d_class.hpp"

void MHD3D::ideal(double dt)
{
  // 3D ideal MHD simulation
  int i,j,k,ss,s1,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  static const int nxy=nx*ny,nxyz=nx*ny*nz;
  const double dtdx=dt/dx,dtdy=dt/dy,dtdz=dt/dz;
  void (*func_flux)(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
		    double, double, const double*,
		    double*, double*, double*, double*, double*, double*, double*)=riemann[RMN];
  void (*lfun_lr)(const double *f, double *fl, double *fr)=l_interp[ODR-1];
  void (*func_lr)(const double *f, double *fl, double *fr)=interpol[ODR-1];
  double (*func_bc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ut,*fx,*fy,*fz,*fc;
  double *ex,*ey,*ez;
  double *ct,*dvx,*dvy,*dvz;

  ut=new double[nm*nxyz];
  fx=new double[nm*nxyz];
  fy=new double[nm*nxyz];
  fz=new double[nm*nxyz];
  fc=new double[2*nxyz];
  ex=new double[nxyz];
  ey=new double[nxyz];
  ez=new double[nxyz];
  ct=new double[3*nxyz];
  dvx=new double[nxyz];
  dvy=new double[nxyz];
  dvz=new double[nxyz];

  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nxyz],val[i],nxyz);
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
	    ss=nx*(ny*k+j)+i;
	    cx[ss]=bcell(&bx[ss],  1,func_bc);
	  }
	}
	for (j=1;j<ny-2;j++){
#pragma simd
	  for (i=0;i<nx;i++){
	    ss=nx*(ny*k+j)+i;
	    cy[ss]=bcell(&by[ss], nx,func_bc);
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
	    ss=nx*(ny*k+j)+i;
	    cz[ss]=bcell(&bz[ss],nxy,func_bc);
	  }
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	int stx[3]={0,0,0},sty[3]={0,0,0},stz[3]={0,0,0};
	double *pc[]={cx,cy,cz};
	bound(pc,3,stx,&dnxs[4],sty,&dnys[4],stz,&dnzs[4]);
      }

      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxyz;ss++){
	prmtv(ss);
	ex[ss]=ey[ss]=ez[ss]=0.0; // Necessary initialize @ cell corner
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
	bound(pv,3,stx,dnx,sty,dny,stz,dnz);
      }

      /* Numerical flux at cell face along X */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur,*ql,*qr;
	ul=new double[nm*nx];
	ur=new double[nm*nx];
	ql=new double[2*nx];
	qr=new double[2*nx];
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
			cx[ss],gam,1,func_lr,vl,vr);
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
	  }
#pragma simd
	  for (i=2;i<nx-2;i++){
	    ss=nx*(ny*k+j)+i;
	    s1l=i+1;
	    s1r=i;
	    double vl,vr;
	    /* Linear interpolation of numerical flux of By */
	    mhd_lr_fb(&vx[ss],&vy[ss],&bx[ss],&cy[ss],1,lfun_lr,&vl,&vr);
	    ql[2*s1l+0]=vl;	/* by*vx-bx*vy @ i+1/2 Left */
	    qr[2*s1r+0]=vr;	/* by*vx-bx*vy @ i-1/2 Right */
	    /* Linear interpolation of numerical flux of Bz */
	    mhd_lr_fb(&vx[ss],&vz[ss],&bx[ss],&cz[ss],1,lfun_lr,&vl,&vr);
	    ql[2*s1l+1]=vl;	/* bz*vx-bx*vz @ i+1/2 Left */
	    qr[2*s1r+1]=vr;	/* bz*vx-bx*vz @ i-1/2 Right */
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
		      bn,gam,dvsd,
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
	delete[] ul;
	delete[] ur;
	delete[] ql;
	delete[] qr;
      }
      /* Numerical flux of By at cell corner along Y to calculate Ez */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=new double[2*ny];
	ur=new double[2*ny];
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
	delete[] ul;
	delete[] ur;
      }
      /* Numerical flux of Bz at cell corner along Z to calculate Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=new double[2*nz];
	ur=new double[2*nz];
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
	delete[] ul;
	delete[] ur;
      }

      /* Numerical flux at cell face along Y */	  
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur,*ql,*qr;
	ul=new double[nm*ny];
	ur=new double[nm*ny];
	ql=new double[2*ny];
	qr=new double[2*ny];
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
			cy[ss],gam,nx,func_lr,vl,vr);
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
	  }
#pragma simd
	  for (j=2;j<ny-2;j++){
	    ss=nx*(ny*k+j)+i;
	    s1l=j+1;
	    s1r=j;
	    double vl,vr;
	    /* Linear interpolation of numerical flux of Bz */
	    mhd_lr_fb(&vy[ss],&vz[ss],&by[ss],&cz[ss],nx,lfun_lr,&vl,&vr);
	    ql[2*s1l+0]=vl;	/* bz*vy-by*vz @ j+1/2 Left */
	    qr[2*s1r+0]=vr;	/* bz*vy-by*vz @ j-1/2 Right */
	    /* Linear interpolation of numerical flux of Bx */
	    mhd_lr_fb(&vy[ss],&vx[ss],&by[ss],&cx[ss],nx,lfun_lr,&vl,&vr);
	    ql[2*s1l+1]=vl;	/* bx*vy-by*vx @ j+1/2 Left */
	    qr[2*s1r+1]=vr;	/* bx*vy-by*vx @ j-1/2 Right */
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
		      bn,gam,dvsd,
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
	delete[] ul;
	delete[] ur;
	delete[] ql;
	delete[] qr;
      }
      /* Numerical flux of Bz at cell corner along Z to calculate Ex */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=3;j<ny-2;j++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=new double[2*nz];
	ur=new double[2*nz];
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
	delete[] ul;
	delete[] ur;
      }
      /* Numerical flux of Bx at cell corner along X to calculate Ez */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=0;k<nz;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=new double[2*nx];
	ur=new double[2*nx];
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
	delete[] ul;
	delete[] ur;
      }

      /* Numerical flux at cell face along Z */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=0;j<ny;j++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur,*ql,*qr;
	ul=new double[nm*nz];
	ur=new double[nm*nz];
	ql=new double[2*nz];
	qr=new double[2*nz];
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
			cz[ss],gam,nxy,func_lr,vl,vr);
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
	  }
#pragma simd
	  for (k=2;k<nz-2;k++){
	    ss=nx*(ny*k+j)+i;
	    s1l=k+1;
	    s1r=k;
	    double vl,vr;
	    /* Linear interpolation of numerical flux of Bx */
	    mhd_lr_fb(&vz[ss],&vx[ss],&bz[ss],&cx[ss],nxy,lfun_lr,&vl,&vr);
	    ql[2*s1l+0]=vl;	/* bx*vz-bz*vx @ k+1/2 Left */
	    qr[2*s1r+0]=vr;	/* bx*vz-bz*vx @ k-1/2 Right */
	    /* Linear interpolation of numerical flux of By */
	    mhd_lr_fb(&vz[ss],&vy[ss],&bz[ss],&cy[ss],nxy,lfun_lr,&vl,&vr);
	    ql[2*s1l+1]=vl;	/* by*vz-bz*vy @ k+1/2 Left */
	    qr[2*s1r+1]=vr;	/* by*vz-bz*vy @ k-1/2 Right */
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
		      bn,gam,dvsd,
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
	delete[] ul;
	delete[] ur;
	delete[] ql;
	delete[] qr;
      }
      /* Numerical flux of Bx at cell corner along X to calculate Ey */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=3;k<nz-2;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=new double[2*nx];
	ur=new double[2*nx];
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
	delete[] ul;
	delete[] ur;
      }
      /* Numerical flux of By at cell corner along Y to calculate Ex */
#ifdef _OPENMP
#pragma omp for
#endif
      for (k=3;k<nz-2;k++){
	int s1l,s1r;
	int sl,sr;
	double *ul,*ur;
	ul=new double[2*ny];
	ur=new double[2*ny];
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
	delete[] ul;
	delete[] ur;
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
	    double *val1[]={val[0]+ss,val[1]+ss,val[2]+ss,val[3]+ss,val[4]+ss,val[5]+ss,val[6]+ss,val[7]+ss};
	    double val0[]={ut[0*nxyz+ss],ut[1*nxyz+ss],ut[2*nxyz+ss],ut[3*nxyz+ss],ut[4*nxyz+ss],ut[5*nxyz+ss],ut[6*nxyz+ss],ut[7*nxyz+ss]};
	    mhd_updt3d(val1[0],val0[0],&fx[nm*ss+0],&fy[nm*ss+0],&fz[nm*ss+0],dtdx,dtdy,dtdz,rk_fac[rk],nm,nm*nx,nm*nxy,func_df); // ro
	    mhd_updt3d(val1[1],val0[1],&fx[nm*ss+1],&fy[nm*ss+1],&fz[nm*ss+1],dtdx,dtdy,dtdz,rk_fac[rk],nm,nm*nx,nm*nxy,func_df); // mx
	    mhd_updt3d(val1[2],val0[2],&fx[nm*ss+2],&fy[nm*ss+2],&fz[nm*ss+2],dtdx,dtdy,dtdz,rk_fac[rk],nm,nm*nx,nm*nxy,func_df); // my
	    mhd_updt3d(val1[3],val0[3],&fx[nm*ss+3],&fy[nm*ss+3],&fz[nm*ss+3],dtdx,dtdy,dtdz,rk_fac[rk],nm,nm*nx,nm*nxy,func_df); // mz
	    mhd_updt3d(val1[7],val0[7],&fx[nm*ss+7],&fy[nm*ss+7],&fz[nm*ss+7],dtdx,dtdy,dtdz,rk_fac[rk],nm,nm*nx,nm*nxy,func_df); // en
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
    bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
  }
      
  delete[] ut;
  delete[] fx;
  delete[] fy;
  delete[] fz;
  delete[] fc;
  delete[] ex;
  delete[] ey;
  delete[] ez;
  delete[] ct;
  delete[] dvx;
  delete[] dvy;
  delete[] dvz;
}
