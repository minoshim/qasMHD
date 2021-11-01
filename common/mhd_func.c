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
  double ros[]={ro[sm2],ro[sm1],ro[ss0],ro[sp1],ro[sp2]};
  double vxs[]={vx[sm2],vx[sm1],vx[ss0],vx[sp1],vx[sp2]};
  double vys[]={vy[sm2],vy[sm1],vy[ss0],vy[sp1],vy[sp2]};
  double vzs[]={vz[sm2],vz[sm1],vz[ss0],vz[sp1],vz[sp2]};
  double bys[]={by[sm2],by[sm1],by[ss0],by[sp1],by[sp2]};
  double bzs[]={bz[sm2],bz[sm1],bz[ss0],bz[sp1],bz[sp2]};
  double prs[]={pr[sm2],pr[sm1],pr[ss0],pr[sp1],pr[sp2]};

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
  double vxs[]={vx[sm2],vx[sm1],vx[ss0],vx[sp1],vx[sp2]};
  double vys[]={vy[sm2],vy[sm1],vy[ss0],vy[sp1],vy[sp2]};
  double bxs[]={bx[sm2],bx[sm1],bx[ss0],bx[sp1],bx[sp2]};
  double bys[]={by[sm2],by[sm1],by[ss0],by[sp1],by[sp2]};
  double val[]={bys[0]*vxs[0],bys[1]*vxs[1],bys[2]*vxs[2],bys[3]*vxs[3],bys[4]*vxs[4]};
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
  double val[]={f[sm2],f[sm1],f[ss0],f[sp1],f[sp2]};
  func_lr(&val[ns/2],vl,vr);
}

void mhd_updt1d(double *val[], const double *val0,
		const double *fx, double dtdx, const double rk_fac[2],
		double (*func_df)(const double*))
/* Update 1D MHD variables */
{
  const int nm=7;		/* Number of variables */
  int sm1,ss0,sp1,sp2;
  sm1=-1;
  ss0=0;
  sp1=+1;
  sp2=+2;
  double du[nm];
  double valx[nm][4]={{fx[nm*sm1+0],fx[nm*ss0+0],fx[nm*sp1+0],fx[nm*sp2+0]},
		      {fx[nm*sm1+1],fx[nm*ss0+1],fx[nm*sp1+1],fx[nm*sp2+1]},
		      {fx[nm*sm1+2],fx[nm*ss0+2],fx[nm*sp1+2],fx[nm*sp2+2]},
		      {fx[nm*sm1+3],fx[nm*ss0+3],fx[nm*sp1+3],fx[nm*sp2+3]},
		      {fx[nm*sm1+4],fx[nm*ss0+4],fx[nm*sp1+4],fx[nm*sp2+4]},
		      {fx[nm*sm1+5],fx[nm*ss0+5],fx[nm*sp1+5],fx[nm*sp2+5]},
		      {fx[nm*sm1+6],fx[nm*ss0+6],fx[nm*sp1+6],fx[nm*sp2+6]}};
  du[0]=func_df(&valx[0][1]);
  du[1]=func_df(&valx[1][1]);
  du[2]=func_df(&valx[2][1]);
  du[3]=func_df(&valx[3][1]);
  du[4]=func_df(&valx[4][1]);
  du[5]=func_df(&valx[5][1]);
  du[6]=func_df(&valx[6][1]);

  rk_updt(val[0],val0[0],-dtdx*du[0],rk_fac[0],rk_fac[1]);
  rk_updt(val[1],val0[1],-dtdx*du[1],rk_fac[0],rk_fac[1]);
  rk_updt(val[2],val0[2],-dtdx*du[2],rk_fac[0],rk_fac[1]);
  rk_updt(val[3],val0[3],-dtdx*du[3],rk_fac[0],rk_fac[1]);
  rk_updt(val[4],val0[4],-dtdx*du[4],rk_fac[0],rk_fac[1]);
  rk_updt(val[5],val0[5],-dtdx*du[5],rk_fac[0],rk_fac[1]);
  rk_updt(val[6],val0[6],-dtdx*du[6],rk_fac[0],rk_fac[1]);
}
		
void mhd_updt2d(double *val[], const double *val0,
		const double *fx, const double *fy,
		double dtdx, double dtdy, const double rk_fac[2], const double flag[],
		int xoffset, int yoffset,
		double (*func_df)(const double*))
/* Update 2D MHD variables @ cell center */
/* Flag = 1 or 0. Update ith-variable if flag[i]=1 */
/* xoffset and yoffset are 1 and nx */
{
  const int nm=8;		/* Number of variables */
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
  double du[2][nm];
  /* dF/dx */
  double valx[nm][4]={{fx[nm*sim+0],fx[nm*si0+0],fx[nm*sip+0],fx[nm*sip2+0]},
		      {fx[nm*sim+1],fx[nm*si0+1],fx[nm*sip+1],fx[nm*sip2+1]},
		      {fx[nm*sim+2],fx[nm*si0+2],fx[nm*sip+2],fx[nm*sip2+2]},
		      {fx[nm*sim+3],fx[nm*si0+3],fx[nm*sip+3],fx[nm*sip2+3]},
		      {fx[nm*sim+4],fx[nm*si0+4],fx[nm*sip+4],fx[nm*sip2+4]},
		      {fx[nm*sim+5],fx[nm*si0+5],fx[nm*sip+5],fx[nm*sip2+5]},
		      {fx[nm*sim+6],fx[nm*si0+6],fx[nm*sip+6],fx[nm*sip2+6]},
		      {fx[nm*sim+7],fx[nm*si0+7],fx[nm*sip+7],fx[nm*sip2+7]}};
  du[0][0]=func_df(&valx[0][1]);
  du[0][1]=func_df(&valx[1][1]);
  du[0][2]=func_df(&valx[2][1]);
  du[0][3]=func_df(&valx[3][1]);
  du[0][4]=func_df(&valx[4][1]);
  du[0][5]=func_df(&valx[5][1]);
  du[0][6]=func_df(&valx[6][1]);
  du[0][7]=func_df(&valx[7][1]);
  /* dG/dy */
  double valy[nm][4]={{fy[nm*sjm+0],fy[nm*sj0+0],fy[nm*sjp+0],fy[nm*sjp2+0]},
		      {fy[nm*sjm+1],fy[nm*sj0+1],fy[nm*sjp+1],fy[nm*sjp2+1]},
		      {fy[nm*sjm+2],fy[nm*sj0+2],fy[nm*sjp+2],fy[nm*sjp2+2]},
		      {fy[nm*sjm+3],fy[nm*sj0+3],fy[nm*sjp+3],fy[nm*sjp2+3]},
		      {fy[nm*sjm+4],fy[nm*sj0+4],fy[nm*sjp+4],fy[nm*sjp2+4]},
		      {fy[nm*sjm+5],fy[nm*sj0+5],fy[nm*sjp+5],fy[nm*sjp2+5]},
		      {fy[nm*sjm+6],fy[nm*sj0+6],fy[nm*sjp+6],fy[nm*sjp2+6]},
		      {fy[nm*sjm+7],fy[nm*sj0+7],fy[nm*sjp+7],fy[nm*sjp2+7]}};
  du[1][0]=func_df(&valy[0][1]);
  du[1][1]=func_df(&valy[1][1]);
  du[1][2]=func_df(&valy[2][1]);
  du[1][3]=func_df(&valy[3][1]);
  du[1][4]=func_df(&valy[4][1]);
  du[1][5]=func_df(&valy[5][1]);
  du[1][6]=func_df(&valy[6][1]);
  du[1][7]=func_df(&valy[7][1]);

  rk_updt(val[0],val0[0],flag[0]*(-dtdx*du[0][0]-dtdy*du[1][0]),rk_fac[0],rk_fac[1]);
  rk_updt(val[1],val0[1],flag[1]*(-dtdx*du[0][1]-dtdy*du[1][1]),rk_fac[0],rk_fac[1]);
  rk_updt(val[2],val0[2],flag[2]*(-dtdx*du[0][2]-dtdy*du[1][2]),rk_fac[0],rk_fac[1]);
  rk_updt(val[3],val0[3],flag[3]*(-dtdx*du[0][3]-dtdy*du[1][3]),rk_fac[0],rk_fac[1]);
  rk_updt(val[4],val0[4],flag[4]*(-dtdx*du[0][4]-dtdy*du[1][4]),rk_fac[0],rk_fac[1]);
  rk_updt(val[5],val0[5],flag[5]*(-dtdx*du[0][5]-dtdy*du[1][5]),rk_fac[0],rk_fac[1]);
  rk_updt(val[6],val0[6],flag[6]*(-dtdx*du[0][6]-dtdy*du[1][6]),rk_fac[0],rk_fac[1]);
  rk_updt(val[7],val0[7],flag[7]*(-dtdx*du[0][7]-dtdy*du[1][7]),rk_fac[0],rk_fac[1]);
}

void mhd_updt3d(double *val[], const double *val0,
		const double *fx, const double *fy, const double *fz,
		double dtdx, double dtdy, double dtdz, const double rk_fac[2], const double flag[],
		int xoffset, int yoffset, int zoffset,
		double (*func_df)(const double*))
/* Update 3D MHD variables @ cell center */
/* Flag = 1 or 0. Update ith-variable if flag[i]=1 */
/* xoffset,yoffset,zoffset are 1,nx,nx*ny */
{
  const int nm=8;		/* Number of variables */
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
  double du[3][nm];
  /* dF/dx */
  double valx[nm][4]={{fx[nm*sim+0],fx[nm*si0+0],fx[nm*sip+0],fx[nm*sip2+0]},
		      {fx[nm*sim+1],fx[nm*si0+1],fx[nm*sip+1],fx[nm*sip2+1]},
		      {fx[nm*sim+2],fx[nm*si0+2],fx[nm*sip+2],fx[nm*sip2+2]},
		      {fx[nm*sim+3],fx[nm*si0+3],fx[nm*sip+3],fx[nm*sip2+3]},
		      {fx[nm*sim+4],fx[nm*si0+4],fx[nm*sip+4],fx[nm*sip2+4]},
		      {fx[nm*sim+5],fx[nm*si0+5],fx[nm*sip+5],fx[nm*sip2+5]},
		      {fx[nm*sim+6],fx[nm*si0+6],fx[nm*sip+6],fx[nm*sip2+6]},
		      {fx[nm*sim+7],fx[nm*si0+7],fx[nm*sip+7],fx[nm*sip2+7]}};
  du[0][0]=func_df(&valx[0][1]);
  du[0][1]=func_df(&valx[1][1]);
  du[0][2]=func_df(&valx[2][1]);
  du[0][3]=func_df(&valx[3][1]);
  du[0][4]=func_df(&valx[4][1]);
  du[0][5]=func_df(&valx[5][1]);
  du[0][6]=func_df(&valx[6][1]);
  du[0][7]=func_df(&valx[7][1]);
  /* dG/dy */
  double valy[nm][4]={{fy[nm*sjm+0],fy[nm*sj0+0],fy[nm*sjp+0],fy[nm*sjp2+0]},
		      {fy[nm*sjm+1],fy[nm*sj0+1],fy[nm*sjp+1],fy[nm*sjp2+1]},
		      {fy[nm*sjm+2],fy[nm*sj0+2],fy[nm*sjp+2],fy[nm*sjp2+2]},
		      {fy[nm*sjm+3],fy[nm*sj0+3],fy[nm*sjp+3],fy[nm*sjp2+3]},
		      {fy[nm*sjm+4],fy[nm*sj0+4],fy[nm*sjp+4],fy[nm*sjp2+4]},
		      {fy[nm*sjm+5],fy[nm*sj0+5],fy[nm*sjp+5],fy[nm*sjp2+5]},
		      {fy[nm*sjm+6],fy[nm*sj0+6],fy[nm*sjp+6],fy[nm*sjp2+6]},
		      {fy[nm*sjm+7],fy[nm*sj0+7],fy[nm*sjp+7],fy[nm*sjp2+7]}};
  du[1][0]=func_df(&valy[0][1]);
  du[1][1]=func_df(&valy[1][1]);
  du[1][2]=func_df(&valy[2][1]);
  du[1][3]=func_df(&valy[3][1]);
  du[1][4]=func_df(&valy[4][1]);
  du[1][5]=func_df(&valy[5][1]);
  du[1][6]=func_df(&valy[6][1]);
  du[1][7]=func_df(&valy[7][1]);
  /* dH/dz */
  double valz[nm][4]={{fz[nm*skm+0],fz[nm*sk0+0],fz[nm*skp+0],fz[nm*skp2+0]},
		      {fz[nm*skm+1],fz[nm*sk0+1],fz[nm*skp+1],fz[nm*skp2+1]},
		      {fz[nm*skm+2],fz[nm*sk0+2],fz[nm*skp+2],fz[nm*skp2+2]},
		      {fz[nm*skm+3],fz[nm*sk0+3],fz[nm*skp+3],fz[nm*skp2+3]},
		      {fz[nm*skm+4],fz[nm*sk0+4],fz[nm*skp+4],fz[nm*skp2+4]},
		      {fz[nm*skm+5],fz[nm*sk0+5],fz[nm*skp+5],fz[nm*skp2+5]},
		      {fz[nm*skm+6],fz[nm*sk0+6],fz[nm*skp+6],fz[nm*skp2+6]},
		      {fz[nm*skm+7],fz[nm*sk0+7],fz[nm*skp+7],fz[nm*skp2+7]}};
  du[2][0]=func_df(&valz[0][1]);
  du[2][1]=func_df(&valz[1][1]);
  du[2][2]=func_df(&valz[2][1]);
  du[2][3]=func_df(&valz[3][1]);
  du[2][4]=func_df(&valz[4][1]);
  du[2][5]=func_df(&valz[5][1]);
  du[2][6]=func_df(&valz[6][1]);
  du[2][7]=func_df(&valz[7][1]);
  
  rk_updt(val[0],val0[0],flag[0]*(-dtdx*du[0][0]-dtdy*du[1][0]-dtdz*du[2][0]),rk_fac[0],rk_fac[1]);
  rk_updt(val[1],val0[1],flag[1]*(-dtdx*du[0][1]-dtdy*du[1][1]-dtdz*du[2][1]),rk_fac[0],rk_fac[1]);
  rk_updt(val[2],val0[2],flag[2]*(-dtdx*du[0][2]-dtdy*du[1][2]-dtdz*du[2][2]),rk_fac[0],rk_fac[1]);
  rk_updt(val[3],val0[3],flag[3]*(-dtdx*du[0][3]-dtdy*du[1][3]-dtdz*du[2][3]),rk_fac[0],rk_fac[1]);
  rk_updt(val[4],val0[4],flag[4]*(-dtdx*du[0][4]-dtdy*du[1][4]-dtdz*du[2][4]),rk_fac[0],rk_fac[1]);
  rk_updt(val[5],val0[5],flag[5]*(-dtdx*du[0][5]-dtdy*du[1][5]-dtdz*du[2][5]),rk_fac[0],rk_fac[1]);
  rk_updt(val[6],val0[6],flag[6]*(-dtdx*du[0][6]-dtdy*du[1][6]-dtdz*du[2][6]),rk_fac[0],rk_fac[1]);
  rk_updt(val[7],val0[7],flag[7]*(-dtdx*du[0][7]-dtdy*du[1][7]-dtdz*du[2][7]),rk_fac[0],rk_fac[1]);
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
