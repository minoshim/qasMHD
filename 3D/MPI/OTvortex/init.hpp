#define IFLD (0)
// Flag for initial condition
// 0,1,2: X-Y plane, Y-Z plane, Z-X plane
// 3,4,5: Transpose of 0,1,2
// 6: 3D

void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X,Y, and Z coordinates
{
  int i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  for (i=0;i<nx;i++){
    // x[i]=(i-xoff+0.5+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
    x[i]=(i-xoff+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
  }
  for (j=0;j<ny;j++){
    // y[j]=(j-yoff+0.5+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
    y[j]=(j-yoff+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
  }
  for (k=0;k<nz;k++){
    // z[k]=(k-zoff+0.5+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz+zmin;
    z[k]=(k-zoff+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz+zmin;
  }
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // Orszag-Tang vortex
  int i,j,k;
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double vx,vy,vz,pr;

	ro[ss]=gam*gam;
	pr=gam;

#if (IFLD == 0)
	// X-Y plane
	vx=-sin(y[j]);
	vy=+sin(x[i]);
	vz=0.0;
	bx[ss]=-sin(y[j]);
	by[ss]=+sin(2*x[i]);
	bz[ss]=0.0;
#elif (IFLD ==1)
	// Y-Z plane
	vy=-sin(z[k]);
	vz=+sin(y[j]);
	vx=0.0;
	by[ss]=-sin(z[k]);
	bz[ss]=+sin(2*y[j]);
	bx[ss]=0.0;
#elif (IFLD ==2)
	// Z-X plane
	vz=-sin(x[i]);
	vx=+sin(z[k]);
	vy=0.0;
	bz[ss]=-sin(x[i]);
	bx[ss]=+sin(2*z[k]);
	by[ss]=0.0;
#elif (IFLD ==3)
	// X-Y plane (transpose)
	vx=+sin(y[j]);
	vy=-sin(x[i]);
	vz=0.0;
	bx[ss]=+sin(2*y[j]);
	by[ss]=-sin(x[i]);
	bz[ss]=0.0;
#elif (IFLD ==4)
	// Y-Z plane (transpose)
	vy=+sin(z[k]);
	vz=-sin(y[j]);
	vx=0.0;
	by[ss]=+sin(2*z[k]);
	bz[ss]=-sin(y[j]);
	bx[ss]=0.0;
#elif (IFLD ==5)
	// Z-X plane (transpose)
	vz=+sin(x[i]);
	vx=-sin(z[k]);
	vy=0.0;
	bz[ss]=+sin(2*x[i]);
	bx[ss]=-sin(z[k]);
	by[ss]=0.0;
#else
	vx=0.5*(+sin(2*z[k])-sin(y[j]));
	vy=0.5*(+sin(3*x[i])-sin(z[k]));
	vz=0.5*(+sin(4*y[j])-sin(x[i]));
	bx[ss]=0.5*(+sin(3*z[k])-sin(y[j]));
	by[ss]=0.5*(+sin(4*x[i])-sin(z[k]));
	bz[ss]=0.5*(+sin(2*y[j])-sin(x[i]));
#endif

	mhd_cnsvt(ro[ss],vx,vy,vz,bx[ss],by[ss],bz[ss],pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
      }
    }
  }
}
