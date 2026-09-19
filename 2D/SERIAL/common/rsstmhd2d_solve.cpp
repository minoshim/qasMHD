#include "rsstmhd2d_class.hpp"

#include <cstddef>
#include <vector>

namespace {

void rsst_correct(double du[8], const double v[3], const double b[3],
                  const double db[3], double rho, double pr, double gam, double ixi)
{
  // PVS (Iijima+19, Eqs. 16, 20-21) or PMS (Eqs. A.8-A.9), extended to MHD.
  // du and db are dt times the uncorrected RHS, BEFORE the RK combination.
  // b and db use the same cell-center interpolation of the CT magnetic field.
  if (ixi == 1.0) return;
  const double v2=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  const double vdm=v[0]*du[1]+v[1]*du[2]+v[2]*du[3];
  const double bdb=b[0]*db[0]+b[1]*db[1]+b[2]*db[2];
  const double dp=(gam-1.0)*(du[7]-vdm+0.5*v2*du[0]-bdb);
  const double a2=gam*pr/rho; // Physical sound speed squared, not reduced.
  const double drho=(1.0-ixi*ixi)*dp/a2;
  du[0]-=drho;
#if RSST_FORM == 0
  // PVS: preserve the velocity RHS by correcting momentum as well as density.
  du[1]-=v[0]*drho;
  du[2]-=v[1]*drho;
  du[3]-=v[2]*drho;
  // (E+P-|B|^2/2)/rho = a2/(gam-1)+|v|^2/2 for the ideal-gas EOS.
  du[7]-=(a2/(gam-1.0)+0.5*v2)*drho;
#else
  // PMS: leave momentum RHS unchanged. At fixed momentum, entropy and B,
  // (dE/dP) = (e+P-rho*|v|^2/2)/(rho*a2); e is thermal energy density.
  du[7]-=(a2/(gam-1.0)-0.5*v2)*drho;
#endif
  // All three magnetic increments, including du[6] (Bz), remain unchanged.
}

}

void RSSTMHD2D::ideal(double dt)
{
  // 2D ideal MHD simulation with RSST
  int i,j,ss,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const int nxy=nx*ny;
  const double dtdx=dt/dx,dtdy=dt/dy;
  const auto func_flux=rsstrie[RMN];
  const auto lfun_lr=l_interp[ODR-1];
  const auto func_lr=interpol[ODR-1];
  const auto func_bc=fcen[ODR-1];
  const auto func_df=df1[ODR-1];

  const std::size_t cell_size=static_cast<std::size_t>(nxy);
  const std::size_t work_size=static_cast<std::size_t>(nm)*cell_size;
  const std::size_t required_size=5*work_size+9*cell_size;
  // Reuse the calling thread's storage; worker threads share the pointers below.
  // Resize only before entering OpenMP regions. Same-thread reentrant calls are unsupported.
  static thread_local std::vector<double> work;
  if (work.size() != required_size) work.resize(required_size);
  double *ut=work.data();
  double *ul=ut+work_size;
  double *ur=ul+work_size;
  double *fx=ur+work_size;
  double *fy=fx+work_size;
  double *ql=fy+work_size;
  double *qr=ql+cell_size;
  double *fc=qr+cell_size;
  double *ez=fc+cell_size;
  double *ct=ez+cell_size;
  double *dvx=ct+cell_size;
  double *dvy=dvx+cell_size;
  double *dbx=dvy+cell_size;
  double *dby=dbx+cell_size;

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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=1;i<nx-2;i++){
	  ss=nx*j+i;
	  cx[ss]=bcell(&bx[ss], 1,func_bc);
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=1;j<ny-2;j++){
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
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

      // RSST factor 1/\xi
#ifdef _OPENMP
#pragma omp for
#endif
      for (ss=0;ss<nxy;ss++){
	rsst_(ss);
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
	  // mhd_lrstate(&ro[ss],&vx[ss],&vy[ss],&vz[ss],&cy[ss],&cz[ss],&pr[ss],
	  // 	      cx[ss],gam,1,func_lr,vl,vr); // Nonlinear reconstruction NOT suitable for RSST-LHLLD
	  mhd_lr_single(&ro[ss], 1,func_lr,&vl[0],&vr[0]);
	  mhd_lr_single(&vx[ss], 1,lfun_lr,&vl[1],&vr[1]);
	  mhd_lr_single(&vy[ss], 1,lfun_lr,&vl[2],&vr[2]);
	  mhd_lr_single(&vz[ss], 1,lfun_lr,&vl[3],&vr[3]);
	  mhd_lr_single(&cy[ss], 1,lfun_lr,&vl[4],&vr[4]);
	  mhd_lr_single(&cz[ss], 1,lfun_lr,&vl[5],&vr[5]);
	  mhd_lr_single(&pr[ss], 1,lfun_lr,&vl[6],&vr[6]);
	  
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
#ifdef _OPENMP
#pragma omp simd private(ss,sl,sr)
#endif
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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  double flux[8]={0},bn=0.5*(ul[nm*ss+4]+ur[nm*ss+4]);
	  double dvsd[2]={(vx[nx*j+i]-vx[nx*j+(i-1)]),min(dvy[nx*j+(i-1)],dvy[nx*j+i])};
	  func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+5],ul[nm*ss+6],ul[nm*ss+7],ixi[ss-1],
		    ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+5],ur[nm*ss+6],ur[nm*ss+7],ixi[ss  ],
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
#ifdef _OPENMP
#pragma omp simd private(ss,sl,sr)
#endif
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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
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
	  // mhd_lrstate(&ro[ss],&vy[ss],&vz[ss],&vx[ss],&cz[ss],&cx[ss],&pr[ss],
	  // 	      cy[ss],gam,nx,func_lr,vl,vr); // Nonlinear reconstruction NOT suitable for RSST-LHLLD
	  mhd_lr_single(&ro[ss],nx,func_lr,&vl[0],&vr[0]);
	  mhd_lr_single(&vy[ss],nx,lfun_lr,&vl[1],&vr[1]);
	  mhd_lr_single(&vz[ss],nx,lfun_lr,&vl[2],&vr[2]);
	  mhd_lr_single(&vx[ss],nx,lfun_lr,&vl[3],&vr[3]);
	  mhd_lr_single(&cz[ss],nx,lfun_lr,&vl[4],&vr[4]);
	  mhd_lr_single(&cx[ss],nx,lfun_lr,&vl[5],&vr[5]);
	  mhd_lr_single(&pr[ss],nx,lfun_lr,&vl[6],&vr[6]);
	  
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
#ifdef _OPENMP
#pragma omp simd private(ss,sl,sr)
#endif
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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=0;i<nx;i++){
	  ss=nx*j+i;
	  double flux[8]={0},bn=0.5*(ul[nm*ss+4]+ur[nm*ss+4]);
	  double dvsd[2]={(vy[nx*j+i]-vy[nx*(j-1)+i]),min(dvx[nx*(j-1)+i],dvx[nx*j+i])};
	  func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+5],ul[nm*ss+6],ul[nm*ss+7],ixi[ss-nx],
		    ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+5],ur[nm*ss+6],ur[nm*ss+7],ixi[ss   ],
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
#ifdef _OPENMP
#pragma omp simd private(ss,sl,sr)
#endif
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
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=3;i<nx-2;i++){
	  ss=nx*j+i;
	  ez[ss]+=+0.5*((ul[2*ss+0]+ur[2*ss+0])+ct[ss]*(ul[2*ss+1]+ur[2*ss+1]));
	}
      }

      /* CT increments: dt times the RHS, without an RK update yet. */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=xoff;i<nx-xoff+1;i++){
	  ss=nx*j+i;
	  double stencil[4]={ez[ss-nx],ez[ss],ez[ss+nx],ez[ss+2*nx]};
	  dbx[ss]=-dtdy*func_df(&stencil[1]);
	}
      }
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff+1;j++){
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  double stencil[4]={ez[ss-1],ez[ss],ez[ss+1],ez[ss+2]};
	  dby[ss]=+dtdx*func_df(&stencil[1]);
	}
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	// Current magnetic boundaries are linear and homogeneous, so the same
	// operators apply to increments. Inhomogeneous/time-dependent boundaries
	// would require their own increment boundary conditions.
	double *db[]={dbx,dby};
	bound(db,2,&stxs[4],&dnxs[4],&stys[4],&dnys[4]);
      }

      /* Conservative increments, PVS/PMS correction, then RK combination. */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  double du[nm]={0};
	  for (int m=0;m<nm;m++){
	    if (m == 4 || m == 5) continue; // Bx, By are updated by CT.
	    const double *f=&fx[nm*ss+m],*g=&fy[nm*ss+m];
	    double fs[4]={f[-nm],f[0],f[nm],f[2*nm]};
	    double gs[4]={g[-nm*nx],g[0],g[nm*nx],g[2*nm*nx]};
	    du[m]=-dtdx*func_df(&fs[1])-dtdy*func_df(&gs[1]);
	  }
	  if (ixi[ss] < 1.0){
	    const double v[3]={vx[ss],vy[ss],vz[ss]};
	    const double b[3]={cx[ss],cy[ss],cz[ss]};
	    const double db[3]={bcell(&dbx[ss],1,func_bc),
	                        bcell(&dby[ss],nx,func_bc),du[6]};
	    rsst_correct(du,v,b,db,ro[ss],pr[ss],gam,ixi[ss]);
	  }
	  // Read all stage quantities above before changing any conserved variable.
	  for (int m=0;m<nm;m++){
	    if (m == 4 || m == 5) continue;
	    rk_updt(&val[m][ss],ut[m*nxy+ss],du[m],rk_fac[rk][0],rk_fac[rk][1]);
	  }
	}
      }
      /* Update CT Bx */
#ifdef _OPENMP
#pragma omp for nowait
#endif
      for (j=yoff;j<ny-yoff;j++){
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=xoff;i<nx-xoff+1;i++){
	  ss=nx*j+i;
	  rk_updt(&bx[ss],ut[4*nxy+ss],dbx[ss],rk_fac[rk][0],rk_fac[rk][1]);
	}
      }
      /* Update CT By */
#ifdef _OPENMP
#pragma omp for
#endif
      for (j=yoff;j<ny-yoff+1;j++){
#ifdef _OPENMP
#pragma omp simd private(ss)
#endif
	for (i=xoff;i<nx-xoff;i++){
	  ss=nx*j+i;
	  rk_updt(&by[ss],ut[5*nxy+ss],dby[ss],rk_fac[rk][0],rk_fac[rk][1]);
	}
      }

    } /* OpenMP */

    /* Boundary condition */
    bound(val,nm,stxs,dnxs,stys,dnys);
  }
}
