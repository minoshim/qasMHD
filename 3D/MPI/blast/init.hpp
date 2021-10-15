void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X,Y, and Z coordinates
{
  int i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  for (i=0;i<nx;i++){
    x[i]=(i-xoff+0.5+((mpi_rank%m_xy)%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
  }
  for (j=0;j<ny;j++){
    y[j]=(j-yoff+0.5+((mpi_rank%m_xy)/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
  }
  for (k=0;k<nz;k++){
    z[k]=(k-zoff+0.5+(mpi_rank/m_xy)*(ZMESH/MNP_Z))*dz+zmin;
  }
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // Blast wave
  int i,j,k;
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	int ss=nx*(ny*k+j)+i;
	double rr,vx,vy,vz,pr;
	rr=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
      
	ro[ss]=(rr <= r_0)?ro1:ro0;
	vx=0.0;
	vy=0.0;
	vz=0.0;
	pr=(rr <= r_0)?pr1:pr0;
      
	bx[ss]=b_0*sin(bthe*dtor)*cos(bphi*dtor);
	by[ss]=b_0*sin(bthe*dtor)*sin(bphi*dtor);
	bz[ss]=b_0*cos(bthe*dtor);

	mhd_cnsvt(ro[ss],vx,vy,vz,bx[ss],by[ss],bz[ss],pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
      }
    }
  }
}
