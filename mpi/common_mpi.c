#include "common_mpi.h"
#include <stdlib.h>
#include "mpi.h"

void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff,
		int dnx, int dny,
		int mpi_rank, int mpi_numx, int mpi_numy)
/* MPI SendRecv for 2D variables */
/* Set dn=0 for Periodic boundary */
/* For Neumann or Dirichlet condition, call mpi_xbc2d and mpi_ybc2d later */
{
  int i,j,n;
  int mpi_tag=0;
  int rankl,rankh;
  int ntot,ntot2;
  MPI_Status r_stat;
  double *fold,*fcpy;

  /* XBC */
  if (dnx == 0){
    rankl=((mpi_rank % mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
    rankh=((mpi_rank % mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
  } else{
    rankl=((mpi_rank % mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
    rankh=((mpi_rank % mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
  }
  if (mpi_numx != 1){
    ntot=nn*ny*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (j=0;j<ny;j++){	/* Transpose */
	for (n=0;n<nn;n++){
	  fold[nn*(ny*i+j)+n]=f[n][nx*j+(xoff+i)];
	  fold[nn*(ny*(2*xoff-1-i)+j)+n]=f[n][nx*j+(nx-xoff-1-i)];
	  fcpy[nn*(ny*i+j)+n]=f[n][nx*j+i];
	  fcpy[nn*(ny*(2*xoff-1-i)+j)+n]=f[n][nx*j+(nx-1-i)];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	for (n=0;n<nn;n++){
	  f[n][nx*j+i]=fcpy[nn*(ny*i+j)+n];
	  f[n][nx*j+(nx-1-i)]=fcpy[nn*(ny*(2*xoff-1-i)+j)+n];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
    if (dnx == 0){
      /* Periodic. avoid communication to myself */
      for (n=0;n<nn;n++){
	for (j=0;j<ny;j++){
	  for (i=0;i<xoff;i++){
	    f[n][nx*j+(nx-1-i)]=f[n][nx*j+(2*xoff-1-i)];
	    f[n][nx*j+i]=f[n][nx*j+(nx-2*xoff+i)];
	  }
	}
      }
    }
  }

  /* YBC */
  if (dny == 0){
    rankl=((mpi_rank / mpi_numx) == 0)?(mpi_rank+mpi_numx*(mpi_numy-1)):(mpi_rank-mpi_numx);
    rankh=((mpi_rank / mpi_numx) == (mpi_numy-1))?(mpi_rank-mpi_numx*(mpi_numy-1)):(mpi_rank+mpi_numx);
  } else{
    rankl=((mpi_rank / mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-mpi_numx);
    rankh=((mpi_rank / mpi_numx) == (mpi_numy-1))?(MPI_PROC_NULL):(mpi_rank+mpi_numx);
  }
  if (mpi_numy != 1){
    ntot=nn*nx*yoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (n=0;n<nn;n++){
	  fold[nn*(nx*j+i)+n]=f[n][nx*(yoff+j)+i];
	  fold[nn*(nx*(2*yoff-1-j)+i)+n]=f[n][nx*(ny-yoff-1-j)+i];
	  fcpy[nn*(nx*j+i)+n]=f[n][nx*j+i];
	  fcpy[nn*(nx*(2*yoff-1-j)+i)+n]=f[n][nx*(ny-1-j)+i];
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (n=0;n<nn;n++){
	  f[n][nx*j+i]=fcpy[nn*(nx*j+i)+n];
	  f[n][nx*(ny-1-j)+i]=fcpy[nn*(nx*(2*yoff-1-j)+i)+n];
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
    if (dny == 0){
      /* Periodic. avoid communication to myself */
      for (n=0;n<nn;n++){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++){
	    f[n][nx*(ny-1-j)+i]=f[n][nx*(2*yoff-1-j)+i];
	    f[n][nx*j+i]=f[n][nx*(ny-2*yoff+j)+i];
	  }
	}
      }
    }
  }
}

void mpi_xbc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D X Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    int i,j;
    /* Left */
    if ((mpi_rank % mpi_numx) == 0){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++) f[nx*j+i]=dn*f[nx*j+(2*xoff-1+st)-i];
      }
    }
    /* Right */
    if ((mpi_rank % mpi_numx) == (mpi_numx-1)){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff-st;i++) f[nx*j+(nx-1-i)]=dn*f[nx*j+(nx-2*xoff+st)+i];
      }
    }
  }
}

void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D Y Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    int i,j;
    /* Left */
    if (mpi_rank/mpi_numx == 0){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++) f[nx*j+i]=dn*f[nx*((2*yoff-1+st)-j)+i];
      }
    }
    /* Right */
    if (mpi_rank/mpi_numx == (mpi_numy-1)){
      for (j=0;j<yoff-st;j++){
	for (i=0;i<nx;i++) f[nx*(ny-1-j)+i]=dn*f[nx*((ny-2*yoff+st)+j)+i];
      }
    }
  }
}

void mpi_sdrv3d(double *f[], int nn, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int dnx, int dny, int dnz,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* MPI SendRecv for 3D variables */
/* Set dn=0 for Periodic boundary */
/* For Neumann or Dirichlet condition, call mpi_x(y,z)bc3d later */
{
  int i,j,k,n;
  int m_xy=mpi_numx*mpi_numy;
  int mpi_tag=0;
  int rankl,rankh;
  int ntot,ntot2;
  MPI_Status r_stat;
  double *fold,*fcpy;

  /* XBC */
  if (dnx == 0){
    rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
    rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
  } else{
    rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
    rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
  }
  if (mpi_numx != 1){
    ntot=nn*ny*nz*xoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (i=0;i<xoff;i++){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){	/* Transpose */
	  for (n=0;n<nn;n++){
	    fold[nn*(ny*(nz*i+k)+j)+n]=f[n][nx*(ny*k+j)+(xoff+i)];
	    fold[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+n]=f[n][nx*(ny*k+j)+(nx-xoff-1-i)];
	    fcpy[nn*(ny*(nz*i+k)+j)+n]=f[n][nx*(ny*k+j)+i];
	    fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+n]=f[n][nx*(ny*k+j)+(nx-1-i)];
	  }
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[nn*(ny*(nz*i+k)+j)+n];
	    f[n][nx*(ny*k+j)+(nx-1-i)]=fcpy[nn*(ny*(nz*(2*xoff-1-i)+k)+j)+n];
	  }
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
    if (dnx == 0){
      /* Periodic. avoid communication to myself */
      for (n=0;n<nn;n++){
	for (k=0;k<nz;k++){
	  for (j=0;j<ny;j++){
	    for (i=0;i<xoff;i++){
	      f[n][nx*(ny*k+j)+(nx-1-i)]=f[n][nx*(ny*k+j)+(2*xoff-1-i)];
	      f[n][nx*(ny*k+j)+i]=f[n][nx*(ny*k+j)+(nx-2*xoff+i)];
	    }
	  }
	}
      }
    }
  }
  
  /* YBC */
  if (dny == 0){
    rankl=(((mpi_rank%m_xy)/mpi_numx) == 0)?(mpi_rank+mpi_numx*(mpi_numy-1)):(mpi_rank-mpi_numx);
    rankh=(((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1))?(mpi_rank-mpi_numx*(mpi_numy-1)):(mpi_rank+mpi_numx);
  } else{
    rankl=(((mpi_rank%m_xy)/mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-mpi_numx);
    rankh=(((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1))?(MPI_PROC_NULL):(mpi_rank+mpi_numx);
  }
  if (mpi_numy != 1){
    ntot=nn*nz*nx*yoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (k=0;k<nz;k++){	/* Transpose */
	  for (n=0;n<nn;n++){
	    fold[nn*(nz*(nx*j+i)+k)+n]=f[n][nx*(ny*k+(yoff+j))+i];
	    fold[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+n]=f[n][nx*(ny*k+(ny-yoff-1-j))+i];
	    fcpy[nn*(nz*(nx*j+i)+k)+n]=f[n][nx*(ny*k+j)+i];
	    fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+n]=f[n][nx*(ny*k+(ny-1-j))+i];
	  }
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[nn*(nz*(nx*j+i)+k)+n];
	    f[n][nx*(ny*k+(ny-1-j))+i]=fcpy[nn*(nz*(nx*(2*yoff-1-j)+i)+k)+n];
	  }
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
    if (dny == 0){
      /* Periodic. avoid communication to myself */
      for (n=0;n<nn;n++){
	for (k=0;k<nz;k++){
	  for (j=0;j<yoff;j++){
	    for (i=0;i<nx;i++){
	      f[n][nx*(ny*k+(ny-1-j))+i]=f[n][nx*(ny*k+(2*yoff-1-j))+i];
	      f[n][nx*(ny*k+j)+i]=f[n][nx*(ny*k+(ny-2*yoff+j))+i];
	    }
	  }
	}
      }
    }
  }

  /* ZBC */
  if (dnz == 0){
    rankl=((mpi_rank/m_xy) == 0)?(mpi_rank+m_xy*(mpi_numz-1)):(mpi_rank-m_xy);
    rankh=((mpi_rank/m_xy) == (mpi_numz-1))?(mpi_rank-m_xy*(mpi_numz-1)):(mpi_rank+m_xy);
  } else{
    rankl=((mpi_rank/m_xy) == 0)?(MPI_PROC_NULL):(mpi_rank-m_xy);
    rankh=((mpi_rank/m_xy) == (mpi_numz-1))?(MPI_PROC_NULL):(mpi_rank+m_xy);
  }
  if (mpi_numz != 1){
    ntot=nn*nx*ny*zoff;
    ntot2=ntot*2;
    fold=(double*)malloc(sizeof(double)*ntot2);
    fcpy=(double*)malloc(sizeof(double)*ntot2);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  for (n=0;n<nn;n++){
	    fold[nn*(nx*(ny*k+j)+i)+n]=f[n][nx*(ny*(zoff+k)+j)+i];
	    fold[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+n]=f[n][nx*(ny*(nz-zoff-1-k)+j)+i];
	    fcpy[nn*(nx*(ny*k+j)+i)+n]=f[n][nx*(ny*k+j)+i];
	    fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+n]=f[n][nx*(ny*(nz-1-k)+j)+i];
	  }
	}
      }
    }
    MPI_Sendrecv(&fold[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 &fcpy[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    MPI_Sendrecv(&fold[ntot],ntot,MPI_DOUBLE,rankh,mpi_tag,
		 &fcpy[0],ntot,MPI_DOUBLE,rankl,mpi_tag,
		 MPI_COMM_WORLD,&r_stat);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[nn*(nx*(ny*k+j)+i)+n];
	    f[n][nx*(ny*(nz-1-k)+j)+i]=fcpy[nn*(nx*(ny*(2*zoff-1-k)+j)+i)+n];
	  }
	}
      }
    }
    free(fold);
    free(fcpy);
  } else{
    if (dnz == 0){
      /* Periodic. avoid communication to myself */
      for (n=0;n<nn;n++){
	for (k=0;k<zoff;k++){
	  for (j=0;j<ny;j++){
	    for (i=0;i<nx;i++){
	      f[n][nx*(ny*(nz-1-k)+j)+i]=f[n][nx*(ny*(2*zoff-1-k)+j)+i];
	      f[n][nx*(ny*k+j)+i]=f[n][nx*(ny*(nz-2*zoff+k)+j)+i];
	    }
	  }
	}
      }
    }
  }

}

void mpi_xbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D X Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    int i,j,k;
    int m_xy=mpi_numx*mpi_numy;
    /* Left */
    if (((mpi_rank%m_xy)%mpi_numx) == 0){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<xoff;i++) f[nx*(ny*k+j)+i]=dn*f[nx*(ny*k+j)+(2*xoff-1+st)-i];
	}
      }
    }
    /* Right */
    if (((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1)){
      for (k=0;k<nz;k++){
        for (j=0;j<ny;j++){
	  for (i=0;i<xoff-st;i++) f[nx*(ny*k+j)+(nx-1-i)]=dn*f[nx*(ny*k+j)+(nx-2*xoff+st)+i];
        }
      }
    }
  }
}

void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D Y Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    int i,j,k;
    int m_xy=mpi_numx*mpi_numy;
    /* Left */
    if (((mpi_rank%m_xy)/mpi_numx) == 0){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+j)+i]=dn*f[nx*(ny*k+(2*yoff-1+st)-j)+i];
	}
      }
    }
    /* Right */
    if (((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1)){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff-st;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+(ny-1-j))+i]=dn*f[nx*(ny*k+(ny-2*yoff+st)+j)+i];
	}
      }
    }
  }
}

void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D Z Dirichlet or Neumann BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1) or Neumann (+1). if dn==0, nothing to do */
{
  if (dn != 0){
    int i,j,k;
    int m_xy=mpi_numx*mpi_numy;
    /* Left */
    if ((mpi_rank/m_xy) == 0){
      for (k=0;k<zoff;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+j)+i]=dn*f[nx*(ny*((2*zoff-1+st)-k)+j)+i];
	}
      }
    }
    /* Right */
    if ((mpi_rank/m_xy) == (mpi_numz-1)){
      for (k=0;k<zoff-st;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*(nz-1-k)+j)+i]=dn*f[nx*(ny*((nz-2*zoff+st)+k)+j)+i];
	}
      }
    }
  }
}
