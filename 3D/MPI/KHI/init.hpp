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
  // KHI
  int i,j,k,ss;
  int m_xy=mpi_numx*mpi_numy;
  double dvy[nz][nx];
  double stim;			// Time at node 0, used for seed of rand_noise
  unsigned seed;
  double dvparax[2]={0,dv};	// Perturbation parameter for vx
  double dvparay[2]={0,dv};	// Perturbation paeameter for vy

  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //Random seed depending on x and z, used for vy perturbation
  seed=(unsigned)(stim*( 1 + ((mpi_rank%m_xy)%mpi_numx) + mpi_numx*(mpi_rank/m_xy) ));
  
  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      dvy[k][i]=0.0;
#if (RANDOM)
      dvy[k][i]+=rand_noise(dvparay,seed); // Multiple mode perturbation
#else
      dvy[k][i]+=dv*sin(2*pi*x[i]/wlen); // Single mode perturbation
#endif    
    }
  }

  //Reset random seed depending on z, used for vx perturbation
  seed=(unsigned)(stim*(1+(mpi_rank/m_xy)));
  srandom(seed);
  
  for (k=0;k<nz;k++){
    double vamp=(global::vamp)+rand_noise(dvparax,seed); // Perturbation to vx is independent of x and y    
    for (j=0;j<ny;j++){
      double angle=0.5*((angle_u+angle_l)+(angle_u-angle_l)*tanh((y[j]-s0)/lambda));
      for (i=0;i<nx;i++){
	double vx,vy,vz,cx,cy,cz,pr;
	double xm,ym,zm;

	ss=nx*(ny*k+j)+i;
	xm=x[i]-0.5*dx;
	ym=y[j]-0.5*dy;
	zm=z[k]-0.5*dz;

	ro[ss]=0.5*((ro_u+ro_l)+(ro_u-ro_l)*tanh((y[j]-s0)/lambda));

	vx=+vamp*tanh((y[j]-s0)/lambda);
	vy=dvy[k][i]*exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda));
	vz=0.0;
	
	if (angle == 90){
	  bx[ss]=cx=0.0;
	  by[ss]=cy=0.0;
	  bz[ss]=cz=b0;
	} else{
	  bx[ss]=cx=b0*cos(angle*dtor);
	  by[ss]=cy=0.0;
	  bz[ss]=cz=b0*sin(angle*dtor);
	}
	pr=pr0;
	
	mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);      
      }
    }
  }
}
