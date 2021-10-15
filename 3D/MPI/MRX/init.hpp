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
  // MRX
  int i,j,k,ss;
  int m_xy=mpi_numx*mpi_numy;
  double para[2]={0,lambda};
  double dvy[nz][nx];
  double stim;			// Time at node 0, used for seed of rand_noise
  unsigned seed;

  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  seed=(unsigned)(stim*( 1 + ((mpi_rank%m_xy)%mpi_numx) + mpi_numx*(mpi_rank/m_xy) ));

  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      dvy[k][i]=0.0;
#if (RANDOM)
      double dvpara[2]={0,dv};
      dvy[k][i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#endif
    }
  }

  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	double vx,vy,vz,cx,cy,cz,pr;
	double xm,ym,zm;

	ss=nx*(ny*k+j)+i;
	xm=x[i]-0.5*dx;
	ym=y[j]-0.5*dy;
	zm=z[k]-0.5*dz;
	
	ro[ss]=(ro0-ro1)*harris_density(y[j],para)+ro1;

	vx=0.0;
	vy=0.0;
	vz=0.0;
	// Perturbation to vy
	vy+=dvy[k][i]*exp(-(y[j]*y[j])/(4*lambda*lambda));

	bx[ss]=b0*harris_field(y[j],para);
	by[ss]=0.0;
	bz[ss]=bg;
	cx=bx[ss];
	cy=by[ss];
	cz=bz[ss];
	pr=(1.0+beta)*(b0*b0+bg*bg)-(cx*cx+cy*cy+cz*cz);
	pr*=0.5;
	
	// Mag field perturbation by Zenitani
	bx[ss]-=b1*(y[j]/lambda)*exp(-(xm*xm+y[j]*y[j])/(4*lambda*lambda));
	by[ss]+=b1*(x[i]/lambda)*exp(-(x[i]*x[i]+ym*ym)/(4*lambda*lambda));
	cx-=b1*(y[j]/lambda)*exp(-(x[i]*x[i]+y[j]*y[j])/(4*lambda*lambda));
	cy+=b1*(x[i]/lambda)*exp(-(x[i]*x[i]+y[j]*y[j])/(4*lambda*lambda));

	mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);      
      }
    }
  }
}
