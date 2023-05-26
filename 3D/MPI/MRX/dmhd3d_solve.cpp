#include "dmhd3d_class.hpp"

void DMHD3D::dsptv(double dt)
{
  // 3D visco-resistive code for MHD
  int i,j,k,ss,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  static const int nxy=nx*ny,nxyz=nx*ny*nz;
  const double idx=1.0/dx,idy=1.0/dy,idz=1.0/dz;
  const double dtdx=dt*idx,dtdy=dt*idy,dtdz=dt*idz;
  double (*func_fc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ut,*fx,*fy,*fz;
  double *ex,*ey,*ez;

  ut=new double[nm*nxyz];
  fx=new double[nm*nxyz];
  fy=new double[nm*nxyz];
  fz=new double[nm*nxyz];
  ex=new double[nxyz];
  ey=new double[nxyz];
  ez=new double[nxyz];

  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nxyz],val[i],nxyz);
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
	vlcty(ss);
      }

      /* Resistive E-field at CT grid */
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
  delete[] ex;
  delete[] ey;
  delete[] ez;
}
