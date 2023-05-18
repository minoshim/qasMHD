#include "mhd3d_class.hpp"

#define NUM (0)
// Flag for initial condition
// 0,1,2: X-Y plane, Y-Z plane, Z-X plane
// 3,4,5: Transpose of 0,1,2
// 6: 3D

void MHD3D::init_()
{
  int i,j,k;
  int isum=0,jsum=0,ksum=0,m;
  for (m=0;m<mpi_ranx;m++){
    isum+=(XMESH+m)/mpi_numx;
  }
  for (m=0;m<mpi_rany;m++){
    jsum+=(YMESH+m)/mpi_numy;
  }
  for (m=0;m<mpi_ranz;m++){
    ksum+=(ZMESH+m)/mpi_numz;
  }
  for (i=0;i<nx;i++) x[i]=(i-xoff+isum)*dx+xmin;
  for (j=0;j<ny;j++) y[j]=(j-yoff+jsum)*dy+ymin;
  for (k=0;k<nz;k++) z[k]=(k-zoff+ksum)*dz+zmin;

  // Orszag-Tang vortex
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;

	ro[ss]=gam*gam;
	pr[ss]=gam;

#if (NUM == 0)
	// X-Y plane
	vx[ss]=-sin(y[j]);
	vy[ss]=+sin(x[i]);
	vz[ss]=0.0;
	bx[ss]=-sin(y[j]);
	by[ss]=+sin(2*x[i]);
	bz[ss]=0.0;
#elif (NUM ==1)
	// Y-Z plane
	vy[ss]=-sin(z[k]);
	vz[ss]=+sin(y[j]);
	vx[ss]=0.0;
	by[ss]=-sin(z[k]);
	bz[ss]=+sin(2*y[j]);
	bx[ss]=0.0;
#elif (NUM ==2)
	// Z-X plane
	vz[ss]=-sin(x[i]);
	vx[ss]=+sin(z[k]);
	vy[ss]=0.0;
	bz[ss]=-sin(x[i]);
	bx[ss]=+sin(2*z[k]);
	by[ss]=0.0;
#elif (NUM ==3)
	// X-Y plane (transpose)
	vx[ss]=+sin(y[j]);
	vy[ss]=-sin(x[i]);
	vz[ss]=0.0;
	bx[ss]=+sin(2*y[j]);
	by[ss]=-sin(x[i]);
	bz[ss]=0.0;
#elif (NUM ==4)
	// Y-Z plane (transpose)
	vy[ss]=+sin(z[k]);
	vz[ss]=-sin(y[j]);
	vx[ss]=0.0;
	by[ss]=+sin(2*z[k]);
	bz[ss]=-sin(y[j]);
	bx[ss]=0.0;
#elif (NUM ==5)
	// Z-X plane (transpose)
	vz[ss]=+sin(x[i]);
	vx[ss]=-sin(z[k]);
	vy[ss]=0.0;
	bz[ss]=+sin(2*x[i]);
	bx[ss]=-sin(z[k]);
	by[ss]=0.0;
#else
	vx[ss]=0.5*(+sin(2*z[k])-sin(y[j]));
	vy[ss]=0.5*(+sin(3*x[i])-sin(z[k]));
	vz[ss]=0.5*(+sin(4*y[j])-sin(x[i]));
	bx[ss]=0.5*(+sin(3*z[k])-sin(y[j]));
	by[ss]=0.5*(+sin(4*x[i])-sin(z[k]));
	bz[ss]=0.5*(+sin(2*y[j])-sin(x[i]));
#endif

	cx[ss]=bx[ss];
	cy[ss]=by[ss];
	cz[ss]=bz[ss];
	
	cnsvt(ss);
      }
    }
  }

  // Boundary condition
  bound(val,nm,stxs,dnxs,stys,dnys,stzs,dnzs);
}
