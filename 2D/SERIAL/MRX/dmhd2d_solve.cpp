#include "dmhd2d_class.hpp"

void DMHD2D::dsptv(double dt)
{
  // 2D visco-resistive code for MHD
  int i,j,ss,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  static const int nxy=nx*ny;
  const double idx=1.0/dx,idy=1.0/dy;
  const double dtdx=dt*idx,dtdy=dt*idy;
  double (*func_fc)(const double *f)=fcen[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];

  double *ut,*fx,*fy;
  double *ex,*ey,*ez;

  ut=new double[nm*nxy];
  fx=new double[nm*nxy];
  fy=new double[nm*nxy];
  ex=new double[nxy];
  ey=new double[nxy];
  ez=new double[nxy];

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

      /* Velocity at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxy;ss++){
	vlcty(ss);
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
  delete[] fx;
  delete[] fy;
  delete[] ex;
  delete[] ey;
  delete[] ez;
}
