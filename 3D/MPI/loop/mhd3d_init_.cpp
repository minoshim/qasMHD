#include "mhd3d_class.hpp"

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
  for (i=0;i<nx;i++) x[i]=(i-xoff+isum+0.5)*dx+xmin;
  for (j=0;j<ny;j++) y[j]=(j-yoff+jsum+0.5)*dy+ymin;
  for (k=0;k<nz;k++) z[k]=(k-zoff+ksum+0.5)*dz+zmin;

  // Blast wave

  // Initial condition parameters
  const double r_0=0.25;       // Radius of imposed high-P cylinder
  const double ro0=1e0;         // Ambient density
  const double ro1=1e0;         // Density in cylinder
  const double pr0=1e0;         // Ambient pressure
  const double pr1=1e2;         // Pressure in cylinder
  const double b_0=10.0;        // Ambient |B|
  const double bthe=90.0;	// B angle relative to z axis
  const double bphi=30.0;       // B angle relative to x axis
  
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double rr=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
	
	ro[ss]=(rr <= r_0)?ro1:ro0;
	vx[ss]=0.0;
	vy[ss]=0.0;
	vz[ss]=0.0;
	bx[ss]=b_0*sin(bthe*dtor)*cos(bphi*dtor);
	by[ss]=b_0*sin(bthe*dtor)*sin(bphi*dtor);
	bz[ss]=b_0*cos(bthe*dtor);
	pr[ss]=(rr <= r_0)?pr1:pr0;
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
