#include "mhd3d_class.hpp"

inline double cal_r(double x, double y, double z, double cost, double sint);
inline double cal_a3(double r, double rad, double a0);

void MHD3D::init_()
{
  // 3D loop advection based on Lee+13
  int i,j,k;
  // Initial condition parameters
  const double ro0=1.0;
  const double pr0=1e0;
  const double v_0=1.0;
  const double a0=1e-3;
  const double rad=0.3;
  const double tilt=atan(0.5);	// Loop tilt around y-axis
  const double theta=0.01; // Advection angle (radian)
  const double vx0=v_0*cos(theta);
  const double vy0=v_0*sin(theta);
  const double vz0=v_0*2.0;
  const double cost=cos(tilt);
  const double sint=sin(tilt);
  const double idx=1.0/dx,idy=1.0/dy,idz=1.0/dz;
  
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;

	ro[ss]=ro0;
	vx[ss]=vx0;
	vy[ss]=vy0;
	vz[ss]=vz0;
	pr[ss]=pr0;

	double rr[4],a3[4];
	// CT Bx
	rr[0]=cal_r(x[i]-0.5*dx,y[j]-1.5*dy,z[k],cost,sint);
	rr[1]=cal_r(x[i]-0.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[2]=cal_r(x[i]-0.5*dx,y[j]+0.5*dy,z[k],cost,sint);
	rr[3]=cal_r(x[i]-0.5*dx,y[j]+1.5*dy,z[k],cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	bx[ss]=+cost*idy*df1[ODR-1](&a3[1]);
	// CT By
	rr[0]=cal_r(x[i]-1.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[1]=cal_r(x[i]-0.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[2]=cal_r(x[i]+0.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	rr[3]=cal_r(x[i]+1.5*dx,y[j]-0.5*dy,z[k],cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	by[ss]=-cost*idx*df1[ODR-1](&a3[1]);
	rr[0]=cal_r(x[i],y[j]-0.5*dy,z[k]-1.5*dz,cost,sint);
	rr[1]=cal_r(x[i],y[j]-0.5*dy,z[k]-0.5*dz,cost,sint);
	rr[2]=cal_r(x[i],y[j]-0.5*dy,z[k]+0.5*dz,cost,sint);
	rr[3]=cal_r(x[i],y[j]-0.5*dy,z[k]+1.5*dz,cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	by[ss]-=sint*idz*df1[ODR-1](&a3[1]);
	// CT Bz
	rr[0]=cal_r(x[i],y[j]-1.5*dy,z[k]-0.5*dz,cost,sint);
	rr[1]=cal_r(x[i],y[j]-0.5*dy,z[k]-0.5*dz,cost,sint);
	rr[2]=cal_r(x[i],y[j]+0.5*dy,z[k]-0.5*dz,cost,sint);
	rr[3]=cal_r(x[i],y[j]+1.5*dy,z[k]-0.5*dz,cost,sint);
	a3[0]=cal_a3(rr[0],rad,a0);
	a3[1]=cal_a3(rr[1],rad,a0);
	a3[2]=cal_a3(rr[2],rad,a0);
	a3[3]=cal_a3(rr[3],rad,a0);
	bz[ss]=+sint*idy*df1[ODR-1](&a3[1]);
      }
    }
  }
  // Conservative variables
  for (k=1;k<nz-1;k++){
    for (j=1;j<ny-1;j++){
      for (i=1;i<nx-1;i++){
	int ss=nx*(ny*k+j)+i;
	bb2cc(ss,fcen[ODR-1],1*stxs[4],nx*stys[5],nx*ny*stzs[6]);
	cnsvt(ss);
      }
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
}

inline double cal_r(double x, double y, double z, double cost, double sint)
{
  double x1=+x*cost+z*sint;
  double x2=y;
  double x3=-x*sint+z*cost;
  // Satisfy periodic boundary condition (see gardiner+08)
  while(x1 > +0.5*cost){
    x1-=cost;
  }
  while(x1 < -0.5*cost){
    x1+=cost;
  }
  while(x2 > +0.5){
    x2-=1.0;
  }
  while(x2 < -0.5){
    x2+=1.0;
  }
  return(sqrt(x1*x1+x2*x2));
}

inline double cal_a3(double r, double rad, double a0)
{
  double ans=0;
  if (r <= rad) ans=a0*(rad-r);
  return(ans);
}
