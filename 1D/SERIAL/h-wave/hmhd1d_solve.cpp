#include "hmhd1d_class.hpp"

void HMHD1D::ideal(double dt)
{
  // 1D ideal MHD simulation
  // "Artificial" electron inertia effect incuded
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

      // E-field correction by "artificial" electron inertia
#ifdef _OPENMP
#pragma omp single
#endif
      if (de != 0){
	double dummy[nx];
	for (i=0;i<nx;i++) dummy[i]=1.0;
	for (int m=5;m<7;m++){
	  double *etmp=new double[2*nx];
	  for (i=3;i<nx-2;i++){
	    ss=i;
	    etmp[0*nx+ss]=etmp[1*nx+ss]=fx[nm*ss+m]; // -ez,+ey
	  }
	  eorg2enew(&etmp[0*nx],&etmp[1*nx],dummy,-dnxs[m]);
	  for (i=3;i<nx-2;i++){
	    ss=i;
	    double de=etmp[1*nx+ss]-etmp[0*nx+ss];   // -ez,+ey
	    double bb=0.5*(val[m][ss-1]+val[m][ss]); // +by,+bz
	    fx[nm*ss+m]=etmp[1*nx+ss];
	    fx[nm*ss+7]+=de*bb;	// This is NOT parallelized
	  }
	  delete[] etmp;
	}
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

void HMHD1D::hall_(double dt)
{
  // 1D Hall term simulation
  // "Artificial" electron inertia effect incuded
  int i,ss,rk;
  static const double rk_fac[3][2]={{0.0,1.0},{0.5+(R_K-2)*0.25,0.5-(R_K-2)*0.25},{1./3.,2./3.}};
  const double dtdx=dt/dx;
  void (*func_flux)(double, double, double, double, double, double, double, 
		    double, double, double, double, double, double, double, 
		    double, double,
		    double*, double*, double*)=&calc_flux_hall_lf;
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
      /* Hall velocity at cell center */
#ifdef _OPENMP
#pragma omp for
#endif
#pragma simd
      for (i=1;i<nx-1;i++){
	ss=i;
	hallv(ss);
      }
#ifdef _OPENMP
#pragma omp single
#endif
      {
	double *pv[3]={hx,hy,hz};
	bound(pv,3,&dnxs[1]);
      }
      /* Primitive variable at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=2;i<nx-2;i++){
	ss=i;
	int sl=i+1,sr=i;
	double vl[7],vr[7];
	mhd_lr_single(&ro[ss],1,func_lr,&vl[0],&vr[0]);
	mhd_lr_single(&hx[ss],1,func_lr,&vl[1],&vr[1]);
	mhd_lr_single(&hy[ss],1,func_lr,&vl[2],&vr[2]);
	mhd_lr_single(&hz[ss],1,func_lr,&vl[3],&vr[3]);
	mhd_lr_single(&cy[ss],1,func_lr,&vl[4],&vr[4]);
	mhd_lr_single(&cz[ss],1,func_lr,&vl[5],&vr[5]);
	mhd_lr_single(&en[ss],1,func_lr,&vl[6],&vr[6]);
	/* Left-face @ i+1/2 */
	ul[nm*sl+0]=vl[0];	/* ro */
	ul[nm*sl+1]=vl[1];	/* hx */
	ul[nm*sl+2]=vl[2];	/* hy */
	ul[nm*sl+3]=vl[3];	/* hz */
	ul[nm*sl+4]=bx[sl];	// bx
	ul[nm*sl+5]=vl[4];	/* by */
	ul[nm*sl+6]=vl[5];	/* bz */
	ul[nm*sl+7]=vl[6];	/* en */
	/* Right-face @ i-1/2 */
	ur[nm*sr+0]=vr[0];	/* ro */
	ur[nm*sr+1]=vr[1];	/* hx */
	ur[nm*sr+2]=vr[2];	/* hy */
	ur[nm*sr+3]=vr[3];	/* hz */
	ur[nm*sr+4]=bx[sr];	// bx
	ur[nm*sr+5]=vr[4];	/* by */
	ur[nm*sr+6]=vr[5];	/* bz */
	ur[nm*sr+7]=vr[6];	/* en */
      }
      /* Numerical flux at cell face */
#ifdef _OPENMP
#pragma omp for
#endif
      for (i=3;i<nx-2;i++){
	ss=i;	
	double flux[8]={0},bn=0.5*(ul[nm*ss+4]+ur[nm*ss+4]);
	double rn=min(ul[nm*ss+0],ur[nm*ss+0]);
	double smax=0.0;
	smax+=max(fabs(ul[nm*ss+1]),fabs(ur[nm*ss+1])); // Hall velocity
	smax+=vphix*fabs(bn/rn); // Whistler velocity
	func_flux(ul[nm*ss+0],ul[nm*ss+1],ul[nm*ss+2],ul[nm*ss+3],ul[nm*ss+5],ul[nm*ss+6],ul[nm*ss+7],
		  ur[nm*ss+0],ur[nm*ss+1],ur[nm*ss+2],ur[nm*ss+3],ur[nm*ss+5],ur[nm*ss+6],ur[nm*ss+7],
		  bn,smax,
		  &flux[5],&flux[6],&flux[7]);
	fx[nm*ss+0]=flux[0];	/* ro */
	fx[nm*ss+1]=flux[1];	/* mx */
	fx[nm*ss+2]=flux[2];	/* my */
	fx[nm*ss+3]=flux[3];	/* mz */
	fx[nm*ss+4]=0;		// bx
	fx[nm*ss+5]=flux[5];	/* by */
	fx[nm*ss+6]=flux[6];	/* bz */
	fx[nm*ss+7]=flux[7];	/* en */
      }

      // E-field correction by "artificial" electron inertia
#ifdef _OPENMP
#pragma omp single
#endif
      if (de != 0){
	double dummy[nx];
	for (i=0;i<nx;i++) dummy[i]=1.0;
	for (int m=5;m<7;m++){
	  double *etmp=new double[2*nx];
	  for (i=3;i<nx-2;i++){
	    ss=i;
	    etmp[0*nx+ss]=etmp[1*nx+ss]=fx[nm*ss+m]; // -ez,+ey
	  }
	  eorg2enew(&etmp[0*nx],&etmp[1*nx],dummy,-dnxs[m]);
	  for (i=3;i<nx-2;i++){
	    ss=i;
	    double de=etmp[1*nx+ss]-etmp[0*nx+ss];   // -ez,+ey
	    double bb=0.5*(val[m][ss-1]+val[m][ss]); // +by,+bz
	    fx[nm*ss+m]=etmp[1*nx+ss];
	    fx[nm*ss+7]+=de*bb;	// This is NOT parallelized
	  }
	  delete[] etmp;
	}
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
	// mhd_updt1d(val1[0],val0[0],&fx[nm*ss+0],dtdx,rk_fac[rk],nm,func_df); // ro
	// mhd_updt1d(val1[1],val0[1],&fx[nm*ss+1],dtdx,rk_fac[rk],nm,func_df); // mx
	// mhd_updt1d(val1[2],val0[2],&fx[nm*ss+2],dtdx,rk_fac[rk],nm,func_df); // my
	// mhd_updt1d(val1[3],val0[3],&fx[nm*ss+3],dtdx,rk_fac[rk],nm,func_df); // mz
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

// [IMPORTANT] Following functions relate electric field Enew and its base value Eorg
// Their relation is (ro-de^2 \nabla^2) Enew = ro*Eorg
// ro: number density (normalized by ro0)
// de: electron inertia length in density of ro0
// Eorg: electric field WITHOUT the electron inertia effect
// Enew: electric field WITH the electron inertia effect
// Reference: Amano, 2015, JCP

void HMHD1D::enew2eorg(double *enew, double *eorg, const double *ro, int dnx)
/* Convert Enew => Eorg */
/* ro is number density and de is electron inertia length (both are normalized) */
/* Boundary condition included */
{
  const double fac=(de*de)/(dx*dx);
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={eorg,enew};
  int dnxs[]={dnx};

  /* Boundary condition for Enew */
  bound(&p[1],1,dnxs);

  for (int i=xoff;i<nx-xoff;i++){
    int ss=i;
    eorg[ss]=enew[ss]-fac*func_d2f(&enew[ss])/ro[ss];
  }
  
  /* Boundary condition for Eorg (output) */
  bound(&p[0],1,dnxs);
}

void HMHD1D::eorg2enew(double *eorg, double *enew, const double *ro, int dnx)
/* Wrapper of emhd_eorg2enew functions */
{
  int cnt;
  cnt=eorg2enew_cg(eorg,enew,ro,dnx);
  if (!cnt) puts("emhd_eorg2enew: Not converged!!");
}

int HMHD1D::eorg2enew_cg(double *eorg, double *enew, const double *ro, int dnx)
/* Convert Eorg => Enew with Conjugate Gradient method */
/* ro is number density and de is electron inertia length (both are normalized) */
/* Boundary condition included */
{
  int i,ss;
  int cnt=0;
  const int cmax=nx;
  const double eps=1e-8;
  const double fac=(de*de)/(dx*dx);
  double numer,denom,coef,anorm,anormf=0.0;
  double *rk,*pk,*ap;
  rk=new double[nx];
  pk=new double[nx];
  ap=new double[nx];
  double (*func_d2f)(const double *f)=df2[ODR-1];
  double *p[]={eorg,enew,pk};
  int dnxs[]={dnx};

  /* Initialize */
  bound(&p[1],1,dnxs); /* Boundary condition for Enew */
  numer=0.0;
  for (i=xoff;i<nx-xoff;i++){
    ss=i;
    rk[ss]=ro[ss]*eorg[ss]-(ro[ss]*enew[ss]-fac*func_d2f(&enew[ss]));
    pk[ss]=rk[ss];
    numer+=rk[ss]*rk[ss];
    anormf+=fabs(rk[ss]);
  }

  /* Iteration */
  do{

    bound(&p[2],1,dnxs); /* Boundary condition for Pk */
    denom=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      ap[ss]=(ro[ss]*pk[ss]-fac*func_d2f(&pk[ss]));
      denom+=pk[ss]*ap[ss];
    }
    coef=(denom == 0)?0:numer/denom;

    anorm=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      enew[ss]+=coef*pk[ss];
      rk[ss]-=coef*ap[ss];
      anorm+=fabs(rk[ss]);
    }
    
    denom=numer;
    numer=0.0;
    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      numer+=rk[ss]*rk[ss];
    }
    coef=(denom == 0)?0:numer/denom;

    for (i=xoff;i<nx-xoff;i++){
      ss=i;
      pk[ss]=rk[ss]+coef*pk[ss];
    }
    
  } while( (anorm > eps*anormf) && (++cnt < cmax) );

  /* Boundary condition for Enew (output) */
  bound(&p[1],1,dnxs);

  delete[] rk;
  delete[] pk;
  delete[] ap;

  return cmax-cnt;		/* If zero, iteration does NOT converge */
}
