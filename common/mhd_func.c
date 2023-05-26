#include "mhd_func.h"
#include <math.h>
#include "common_func.h"
#include "mhd_eigen.h"

double mhd_ieos(double ro, double vx, double vy, double vz,
		double bx, double by, double bz, double en,
		double gamma)
/* Equation of State for Ideal MHD */
{
  double v2,b2,pr;
  v2=vx*vx+vy*vy+vz*vz;
  b2=bx*bx+by*by+bz*bz;
  pr=(gamma-1)*(en-0.5*(ro*v2+b2));
  return pr;
}

double mhd_energy(double ro, double vx, double vy, double vz,
		  double bx, double by, double bz, double pr,
		  double gamma)
/* Calculate MHD total energy */
{
  double v2,b2,en;
  v2=vx*vx+vy*vy+vz*vz;
  b2=bx*bx+by*by+bz*bz;
  en=pr/(gamma-1)+0.5*(ro*v2+b2);
  return en;
}

void mhd_prmtv(double ro, double mx, double my, double mz,
	       double bx, double by, double bz, double en,
	       double *vx, double *vy, double *vz, double *pr,
	       double gamma)
/* Calculate MHD primitive variables */
{
  double iro;
  iro=1.0/ro;
  (*vx)=mx*iro;
  (*vy)=my*iro;
  (*vz)=mz*iro;
  (*pr)=mhd_ieos(ro,*vx,*vy,*vz,bx,by,bz,en,gamma);
}

void mhd_cnsvt(double ro, double vx, double vy, double vz,
	       double bx, double by, double bz, double pr,
	       double *mx, double *my, double *mz, double *en,
	       double gamma)
/* Calculate MHD conservative variables */
{
  (*mx)=ro*vx;
  (*my)=ro*vy;
  (*mz)=ro*vz;
  (*en)=mhd_energy(ro,vx,vy,vz,bx,by,bz,pr,gamma);
}

double mhd_cuct_weight(const double *ro, const double *vx, const double *vy,
		       const double *bx, const double *by,
		       int xoffset, int yoffset)
/* Weight for Central Upwind CT (Minoshima+19, Eq. (42)) */
/* Upwinding with respect to Alfven mode in x-y plane for Ez */
/* ro,vx,vy @ cell center, bx and by @ cell face */
/* xoffset and yoffset are (1,nx) for x-y plane, (nx,nx*ny) for y-z plane, (nx*ny,1) for z-x plane */
{
  double eps=1e-12;
  double roc,vxc,vyc,bxc,byc;
  int ss[4]={-xoffset-yoffset,-yoffset,-xoffset,0};
  roc=0.25*(ro[ss[0]]+ro[ss[1]]+ro[ss[2]]+ro[ss[3]]);
  roc=1.0/sqrt(roc);
  vxc=0.25*(vx[ss[0]]+vx[ss[1]]+vx[ss[2]]+vx[ss[3]]);
  vyc=0.25*(vy[ss[0]]+vy[ss[1]]+vy[ss[2]]+vy[ss[3]]);
  bxc=0.5*(bx[ss[1]]+bx[ss[3]]);
  byc=0.5*(by[ss[2]]+by[ss[3]]);
  vxc=fabs(vxc)+fabs(bxc*roc);
  vyc=fabs(vyc)+fabs(byc*roc);
  return (vxc+0.5*eps)/(vxc+vyc+eps);
}

void mhd_ct_eres(const double *bx, const double *by, const double *bz,
		 const double *eta,
		 double idx, double idy, double idz,
		 int xoffset, int yoffset, int zoffset,
		 double (*func_df)(const double *f),
		 double *ex, double *ey, double *ez)
/* Calculate resistive E-field (eta*j) @ CT grid */
/* B are defined @ Bx(i-1/2,j,k), By(i,j-1/2,k), Bz(i,j,k-1/2) */
/* E are defined @ Ex(i,j-1/2,k-1/2), Ey(i-1/2,j,k-1/2), Ez(i-1/2,j-1/2,k) */
/* eta is resistivity defined @ cell center (i,j,k) */
/* idx=1.0/dx */
/* xoffset=1, yoffset=nx (or 0 in 1D), zoffset=nx*ny (or 0 in 2D) */
{
  int si[4]={-2*xoffset,-xoffset,0,+xoffset}; /* i-2,i-1,i,i+1 */
  int sj[4]={-2*yoffset,-yoffset,0,+yoffset}; /* j-2,j-1,j,j+1 */
  int sk[4]={-2*zoffset,-zoffset,0,+zoffset}; /* k-2,k-1,k,k+1 */
  double bxj[4]={bx[sj[0]],bx[sj[1]],bx[sj[2]],bx[sj[3]]};
  double bxk[4]={bx[sk[0]],bx[sk[1]],bx[sk[2]],bx[sk[3]]};
  double byk[4]={by[sk[0]],by[sk[1]],by[sk[2]],by[sk[3]]};
  double byi[4]={by[si[0]],by[si[1]],by[si[2]],by[si[3]]};
  double bzi[4]={bz[si[0]],bz[si[1]],bz[si[2]],bz[si[3]]};
  double bzj[4]={bz[sj[0]],bz[sj[1]],bz[sj[2]],bz[sj[3]]};
  double jx=func_df(&bzj[1])*idy-func_df(&byk[1])*idz; /* i,j-1/2,k-1/2 */
  double jy=func_df(&bxk[1])*idz-func_df(&bzi[1])*idx; /* i-1/2,j,k-1/2 */
  double jz=func_df(&byi[1])*idx-func_df(&bxj[1])*idy; /* i-1/2,j-1/2,k */
  *ex=0.25*(eta[-yoffset-zoffset]+eta[-zoffset]+eta[-yoffset]+eta[0])*jx;
  *ey=0.25*(eta[-zoffset-xoffset]+eta[-xoffset]+eta[-zoffset]+eta[0])*jy;
  *ez=0.25*(eta[-xoffset-yoffset]+eta[-yoffset]+eta[-xoffset]+eta[0])*jz;
}

void mhd_lrstate(const double *ro, const double *vx, const double *vy, const double *vz,
		 const double *by, const double *bz, const double *pr,
		 double bx, double gamma, int offset,
		 void (*func_lr)(const double*, double*, double*),
		 double *vl, double *vr)
/* Calculate left and right state at the interface  */
/* offset should be 1 (in x), nx (in y), nx*ny (in z) */
{
  const int ns=5;		/* Number of stencil */
  const double eps=1e-12;
  int sm2,sm1,ss0,sp1,sp2;
  sm2=-2*offset;
  sm1=-1*offset;
  ss0=0;
  sp1=+1*offset;
  sp2=+2*offset;
  /* Primitive variables in the stencil */
  double ros[ns]={ro[sm2],ro[sm1],ro[ss0],ro[sp1],ro[sp2]};
  double vxs[ns]={vx[sm2],vx[sm1],vx[ss0],vx[sp1],vx[sp2]};
  double vys[ns]={vy[sm2],vy[sm1],vy[ss0],vy[sp1],vy[sp2]};
  double vzs[ns]={vz[sm2],vz[sm1],vz[ss0],vz[sp1],vz[sp2]};
  double bys[ns]={by[sm2],by[sm1],by[ss0],by[sp1],by[sp2]};
  double bzs[ns]={bz[sm2],bz[sm1],bz[ss0],bz[sp1],bz[sp2]};
  double prs[ns]={pr[sm2],pr[sm1],pr[ss0],pr[sp1],pr[sp2]};

  /* Characteristic decomposition for supersonic stencil */
  /* int flg[3]={(ros[ns/2-1]*(vxs[ns/2-1]*vxs[ns/2-1]+vys[ns/2-1]*vys[ns/2-1]+vzs[ns/2-1]*vzs[ns/2-1]) > gamma*prs[ns/2-1]), */
  /* 	      (ros[ns/2+0]*(vxs[ns/2+0]*vxs[ns/2+0]+vys[ns/2+0]*vys[ns/2+0]+vzs[ns/2+0]*vzs[ns/2+0]) > gamma*prs[ns/2+0]), */
  /* 	      (ros[ns/2+1]*(vxs[ns/2+1]*vxs[ns/2+1]+vys[ns/2+1]*vys[ns/2+1]+vzs[ns/2+1]*vzs[ns/2+1]) > gamma*prs[ns/2+1])}; */
  /* if (flg[0]+flg[1]+flg[2] != 0){ */
  /*   /\* mhd_c_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr); *\/ */
  /*   mhd_m_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr); */
  /* } else{ */
  /*   mhd_a_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr); */
  /* } */
  int flg=(ros[ns/2+0]*(vxs[ns/2+0]*vxs[ns/2+0]+vys[ns/2+0]*vys[ns/2+0]+vzs[ns/2+0]*vzs[ns/2+0]) > gamma*prs[ns/2+0]);
  if (flg){
    mhd_m_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr);
  } else{
    mhd_a_reconst(ros,vxs,vys,vzs,bys,bzs,prs,bx,gamma,ns,vl,vr,func_lr);
  }

  /* Positivity preservation */
  if (vl[0] <= eps || vl[6] <= eps){
    vl[0]=ros[ns/2];
    vl[1]=vxs[ns/2];
    vl[6]=prs[ns/2];
  }
  if (vr[0] <= eps || vr[6] <= eps){
    vr[0]=ros[ns/2];
    vr[1]=vxs[ns/2];
    vr[6]=prs[ns/2];
  }
}

void mhd_lr_fb(const double *vx, const double *vy, const double *bx, const double *by,
	       int offset,
	       void (*func_lr)(const double*, double*, double*),
	       double *vl, double *vr)
/* Calculate left and right state of numerical flux of B, by*vx-bx*vy, required for CUCT  */
/* Bx @ cell face, By @ cell center */
/* offset should be 1 (in x), nx (in y), nx*ny (in z) */
{
  const int ns=5;		/* Number of stencil */
  int sm2,sm1,ss0,sp1,sp2;
  sm2=-2*offset;
  sm1=-1*offset;
  ss0=0;
  sp1=+1*offset;
  sp2=+2*offset;
  double vxs[ns]={vx[sm2],vx[sm1],vx[ss0],vx[sp1],vx[sp2]};
  double vys[ns]={vy[sm2],vy[sm1],vy[ss0],vy[sp1],vy[sp2]};
  double bxs[ns]={bx[sm2],bx[sm1],bx[ss0],bx[sp1],bx[sp2]};
  double bys[ns]={by[sm2],by[sm1],by[ss0],by[sp1],by[sp2]};
  double val[ns]={bys[0]*vxs[0],bys[1]*vxs[1],bys[2]*vxs[2],bys[3]*vxs[3],bys[4]*vxs[4]};
  double fl[2],fr[2];
  func_lr(&val[ns/2],&fl[0],&fr[0]);
  func_lr(&vys[ns/2],&fl[1],&fr[1]);
  (*vl)=fl[0]-bxs[ns/2+1]*fl[1];
  (*vr)=fr[0]-bxs[ns/2+0]*fr[1];
}

void mhd_lr_single(const double *f, int offset,
		   void (*func_lr)(const double*, double*, double*),
		   double *vl, double *vr)
/* Calculate left and right state of single variable  */
/* offset should be 1 (in x), nx (in y), nx*ny (in z) */
{
  const int ns=5;		/* Number of stencil */
  int sm2,sm1,ss0,sp1,sp2;
  sm2=-2*offset;
  sm1=-1*offset;
  ss0=0;
  sp1=+1*offset;
  sp2=+2*offset;
  double val[ns]={f[sm2],f[sm1],f[ss0],f[sp1],f[sp2]};
  func_lr(&val[ns/2],vl,vr);
}

void mhd_updt1d(double *val, double val0, const double *fx,
		double dtdx, const double rk_fac[2], int xoffset,
		double (*func_df)(const double*))
/* Update 1D MHD variables */
{
  int sm1,ss0,sp1,sp2;
  sm1=-1*xoffset;
  ss0=0;
  sp1=+1*xoffset;
  sp2=+2*xoffset;
  double valx[4]={fx[sm1],fx[ss0],fx[sp1],fx[sp2]};
  double du=func_df(&valx[1]);
  rk_updt(val,val0,-dtdx*du,rk_fac[0],rk_fac[1]);
}

void mhd_updt2d(double *val, double val0,
		const double *fx, const double *fy,
		double dtdx, double dtdy, const double rk_fac[2],
		int xoffset, int yoffset,
		double (*func_df)(const double*))
/* Update 2D MHD variables @ cell center */
{
  int sim,si0,sip,sip2;
  int sjm,sj0,sjp,sjp2;
  sim=-xoffset;
  si0=0;
  sip=+xoffset;
  sip2=2*xoffset;
  sjm=-yoffset;
  sj0=0;
  sjp=+yoffset;
  sjp2=2*yoffset;
  double du[2];
  /* dF/dx */
  double valx[4]={fx[sim],fx[si0],fx[sip],fx[sip2]};
  du[0]=func_df(&valx[1]);
  /* dG/dy */
  double valy[4]={fy[sjm],fy[sj0],fy[sjp],fy[sjp2]};
  du[1]=func_df(&valy[1]);

  rk_updt(val,val0,(-dtdx*du[0]-dtdy*du[1]),rk_fac[0],rk_fac[1]);
}

void mhd_updt3d(double *val, double val0,
		const double *fx, const double *fy, const double *fz,
		double dtdx, double dtdy, double dtdz, const double rk_fac[2],
		int xoffset, int yoffset, int zoffset,
		double (*func_df)(const double*))
/* Update 3D MHD variables @ cell center */
{
  int sim,si0,sip,sip2;
  int sjm,sj0,sjp,sjp2;
  int skm,sk0,skp,skp2;
  sim=-xoffset;
  si0=0;
  sip=+xoffset;
  sip2=2*xoffset;
  sjm=-yoffset;
  sj0=0;
  sjp=+yoffset;
  sjp2=2*yoffset;
  skm=-zoffset;
  sk0=0;
  skp=+zoffset;
  skp2=2*zoffset;
  double du[3];
  /* dF/dx */
  double valx[4]={fx[sim],fx[si0],fx[sip],fx[sip2]};
  du[0]=func_df(&valx[1]);
  /* dG/dy */
  double valy[4]={fy[sjm],fy[sj0],fy[sjp],fy[sjp2]};
  du[1]=func_df(&valy[1]);
  /* dH/dz */
  double valz[4]={fz[skm],fz[sk0],fz[skp],fz[skp2]};
  du[2]=func_df(&valz[1]);
  
  rk_updt(val,val0,(-dtdx*du[0]-dtdy*du[1]-dtdz*du[2]),rk_fac[0],rk_fac[1]);
}

void mhd_updt2d_ctb(double *val, double val0, const double *ez,
		    double dtdy, const double rk_fac[2], int offset,
		    double (*func_df)(const double*))
/* Update 2D MHD B @ cell face (bx or by) by CT method */
/* dtdy should include the direction. For exsample, it is negative when val=by */
{
  int sm1,ss0,sp1,sp2;
  sm1=-1*offset;
  ss0=0;
  sp1=+1*offset;
  sp2=+2*offset;
  double ezs[4]={ez[sm1],ez[ss0],ez[sp1],ez[sp2]};
  double du=func_df(&ezs[1]);
  rk_updt(val,val0,-dtdy*du,rk_fac[0],rk_fac[1]);
}
		    
void mhd_updt3d_ctb(double *bx, double bx0, const double *ey, const double *ez,
		    double dtdy, double dtdz, const double rk_fac[2],
		    int yoffset, int zoffset,
		    double (*func_df)(const double*))
/* Update 3D MHD Bx @ cell face by CT method */
/* yoffset and zoffset are nx and nx*ny */
{
  int sjm1,sj00,sjp1,sjp2;
  int skm1,sk00,skp1,skp2;
  sjm1=-1*yoffset;
  sj00=0;
  sjp1=+1*yoffset;
  sjp2=+2*yoffset;
  skm1=-1*zoffset;
  sk00=0;
  skp1=+1*zoffset;
  skp2=+2*zoffset;
  double eys[4]={ey[skm1],ey[sk00],ey[skp1],ey[skp2]};
  double ezs[4]={ez[sjm1],ez[sj00],ez[sjp1],ez[sjp2]};
  double du[2];
  du[0]=+func_df(&ezs[1]);
  du[1]=-func_df(&eys[1]);
  rk_updt(bx,bx0,-dtdy*du[0]-dtdz*du[1],rk_fac[0],rk_fac[1]);
}

