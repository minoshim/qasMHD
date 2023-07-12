#include "mhd1d_class.hpp"

void MHD1D::ideal(double dt)
{
  // 1D ideal MHD simulation
  int i,ss,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const double dtdx=dt/dx;
  void (*func_flux)(double, double, double, double, double, double, double, 
		    double, double, double, double, double, double, double, 
		    double, double, const double*,
		    double*, double*, double*, double*, double*, double*, double*)=riemann[RMN];
  void (*func_lr)(const double *f, double *fl, double *fr)=interpol[ODR-1];
  double (*func_df)(const double *f)=df1[ODR-1];
  double *ut,*ul,*ur,*fx;
  
  ut=new double[nm*nx];
  ul=new double[nm*nx];
  ur=new double[nm*nx];
  fx=new double[nm*nx];

  // B @ cell center
  cx=bx;
  cy=by;
  cz=bz;
  
  /* Copy current data */
  for (i=0;i<nm;i++){
    cpy_array(&ut[i*nx],val[i],nx);
  }

  /* Runge-kutta stage */
  for (rk=0;rk<R_K;rk++){
#ifdef _OPENMP
#pragma omp parallel private(i,ss)
#endif
    {
      /* Primitive variable at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=0;i<nx;i++){
	ss=i;
	prmtv(ss);
      }
      /* Primitive variable at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=2;i<nx-2;i++){
	ss=i;
	int sl=i+1,sr=i;
	double vl[7],vr[7];
	mhd_lrstate(&ro[ss],&vx[ss],&vy[ss],&vz[ss],&cy[ss],&cz[ss],&pr[ss],
		    cx[ss],gam,1,func_lr,vl,vr);
	/* Left-face @ i+1/2 */
	ul[nm*sl+0]=vl[0];	/* ro */
	ul[nm*sl+1]=vl[1];	/* vx */
	ul[nm*sl+2]=vl[2];	/* vy */
	ul[nm*sl+3]=vl[3];	/* vz */
	ul[nm*sl+4]=bx[sl];	// bx
	ul[nm*sl+5]=vl[4];	/* by */
	ul[nm*sl+6]=vl[5];	/* bz */
	ul[nm*sl+7]=vl[6];	/* pr */
	/* Right-face @ i-1/2 */
	ur[nm*sr+0]=vr[0];	/* ro */
	ur[nm*sr+1]=vr[1];	/* vx */
	ur[nm*sr+2]=vr[2];	/* vy */
	ur[nm*sr+3]=vr[3];	/* vz */
	ur[nm*sr+4]=bx[sr];	// bx
	ur[nm*sr+5]=vr[4];	/* by */
	ur[nm*sr+6]=vr[5];	/* bz */
	ur[nm*sr+7]=vr[6];	/* pr */
      }
      /* Numerical flux at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=3;i<nx-2;i++){
	ss=i;	
	double flux[8]={0},bn=0.5*(ul[nm*ss+4]+ur[nm*ss+4]);
	double dvsd[2]={(vx[i]-vx[i-1]),0.0};
	func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+5],ul[nm*ss+6],ul[nm*ss+7],
		  ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+5],ur[nm*ss+6],ur[nm*ss+7],
		  bn,gam,dvsd,
		  &flux[0],&flux[1],&flux[2],&flux[3],&flux[5],&flux[6],&flux[7]);
	fx[nm*ss+0]=flux[0];	/* ro */
	fx[nm*ss+1]=flux[1];	/* mx */
	fx[nm*ss+2]=flux[2];	/* my */
	fx[nm*ss+3]=flux[3];	/* mz */
	fx[nm*ss+4]=0;		// bx
	fx[nm*ss+5]=flux[5];	/* by */
	fx[nm*ss+6]=flux[6];	/* bz */
	fx[nm*ss+7]=flux[7];	/* en */
      }
      /* Update */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=xoff;i<nx-xoff;i++){
	ss=i;
	double *val1[]={val[0]+ss,val[1]+ss,val[2]+ss,val[3]+ss,val[4]+ss,val[5]+ss,val[6]+ss,val[7]+ss};
	double val0[]={ut[0*nx+ss],ut[1*nx+ss],ut[2*nx+ss],ut[3*nx+ss],ut[4*nx+ss],ut[5*nx+ss],ut[6*nx+ss],ut[7*nx+ss]};
	mhd_updt1d(val1[0],val0[0],&fx[nm*ss+0],dtdx,rk_fac[rk],nm,func_df); // ro
	mhd_updt1d(val1[1],val0[1],&fx[nm*ss+1],dtdx,rk_fac[rk],nm,func_df); // mx
	mhd_updt1d(val1[2],val0[2],&fx[nm*ss+2],dtdx,rk_fac[rk],nm,func_df); // my
	mhd_updt1d(val1[3],val0[3],&fx[nm*ss+3],dtdx,rk_fac[rk],nm,func_df); // mz
	mhd_updt1d(val1[5],val0[5],&fx[nm*ss+5],dtdx,rk_fac[rk],nm,func_df); // by
	mhd_updt1d(val1[6],val0[6],&fx[nm*ss+6],dtdx,rk_fac[rk],nm,func_df); // bz
	mhd_updt1d(val1[7],val0[7],&fx[nm*ss+7],dtdx,rk_fac[rk],nm,func_df); // en
      }
      
    } // OpenMP
    // Boundary condition
    bound(val,nm,dnxs);
  }

  delete[] ut;
  delete[] ul;
  delete[] ur;
  delete[] fx;
}
