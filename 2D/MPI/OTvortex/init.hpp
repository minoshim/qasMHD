void init_grid(int mpi_rank);
void init_plasma(int mpi_rank);

void init_grid(int mpi_rank)
// Define X and Y coordinates
{
  int i,j;
  for (i=0;i<nx;i++){
    // x[i]=(i-xoff+0.5+(mpi_rank%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
    x[i]=(i-xoff+(mpi_rank%mpi_numx)*(XMESH/MNP_X))*dx+xmin;
  }
  for (j=0;j<ny;j++){
    // y[j]=(j-yoff+0.5+(mpi_rank/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
    y[j]=(j-yoff+(mpi_rank/mpi_numx)*(YMESH/MNP_Y))*dy+ymin;
  }
}

void init_plasma(int mpi_rank)
// Set initial condition
{
  // Orszag-Tang vortex
  int i,j;
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double vx,vy,vz,pr;
      ro[ss]=gam*gam;
      vx=-sin(y[j]);
      vy=+sin(x[i]);
      vz=0.0;
      bx[ss]=-sin(y[j]);
      by[ss]=+sin(2*x[i]);
      bz[ss]=0.0;
      pr=gam;

      mhd_cnsvt(ro[ss],vx,vy,vz,bx[ss],by[ss],bz[ss],pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
    }
  }
}
