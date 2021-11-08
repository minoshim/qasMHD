void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);
double cal_pressure(double y0, double y1, double pr0, int n);

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
  // RTI
  int i,j,k,ss;
  int m_xy=mpi_numx*mpi_numy;
  double dvz[ny][nx];
  double stim;			// Time at node 0, used for seed of rand_noise
  unsigned seed;
  double prmin=1e-4;

  if (mpi_rank == 0) stim=MPI_Wtime();
  MPI_Bcast(&stim,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  seed=(unsigned)(stim*(1+(mpi_rank%m_xy)));

  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      dvz[j][i]=0.0;
#if (RANDOM)
      double dvpara[2]={0,dv};
      dvz[j][i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#else
      // dvz[j][i]+=dv*cos(2*pi*x[i]/lx); // Single mode perturbation
      // dvz[j][i]+=dv*cos(2*pi*y[j]/ly); // Single mode perturbation
      dvz[j][i]+=dv*cos(2*pi*x[i]/wlen)*cos(2*pi*y[j]/wlen); // Single mode perturbation
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
	
	ro[ss]=sqrwave2(ro_u,ro_l,(z[k]-s0)/lambda,(z[k]+s0)/lambda);
	vx=sqrwave2(vamp,-vamp,(z[k]-s0)/lambda,(z[k]+s0)/lambda);
	vy=0.0;	
	vz=dvz[j][i]*(exp(-((z[k]-s0)*(z[k]-s0))/(4*lambda*lambda))-exp(-((z[k]+s0)*(z[k]+s0))/(4*lambda*lambda)));

	bx[ss]=cx=b0/sqrt(2.0);
	by[ss]=cy=b0/sqrt(2.0);
	bz[ss]=cz=0.0;
	
	// Gravitational potential, phi_g = g0*lg*log(cosh(z/lg))
	// d(phi_g)/dz = g0*tanh(z/lg)
	phi_g[ss]=g_potential(z[k],g0,lg,0);
      
	// Pressure is obtained by integrating the equilibrium
	pr=cal_pressure(0,z[k],pr0,1+(int)(fabs(z[k])/dz+0.5));
	pr=max(pr,prmin);

	mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
	en[ss]+=ro[ss]*phi_g[ss];	// Add gravitational potential
      }
    }
  }
}

double cal_pressure(double y0, double y1, double pr0, int n)
// Calculate pressure profile from hydrostatic equilibrium dP/dz = -rho dPhi_g/dz
{
  int j;
  double dyc=(y1-y0)/n;
  double ans=pr0;
  for (j=0;j<n;j++){
    double yc,roc,dphi_gc;
    yc=y0+(j+0.5)*dyc;
    roc=sqrwave2(ro_u,ro_l,(yc-s0)/lambda,(yc+s0)/lambda);
    dphi_gc=g_potential(yc,g0,lg,1);
    ans+=-roc*dphi_gc*dyc;
  }
  return ans;
}
