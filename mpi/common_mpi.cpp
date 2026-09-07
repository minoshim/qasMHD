#include "common_mpi.h"
#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <vector>
#include "mpi.h"

namespace {

[[noreturn]] void mpi_abort_all(int error_code)
{
  MPI_Abort(MPI_COMM_WORLD,error_code);
  std::abort();
}

std::size_t checked_mpi_count(std::initializer_list<int> factors)
{
  const std::size_t max_count=static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t count=1;

  for (int factor : factors){
    if ((factor <= 0) || (count > max_count/static_cast<std::size_t>(factor))){
      mpi_abort_all(MPI_ERR_COUNT);
    }
    count*=static_cast<std::size_t>(factor);
  }
  return(count);
}

std::size_t buffer_index_2d(int n, int nn, int outer, int inner, int inner_size)
{
  return(static_cast<std::size_t>(nn)*
	 (static_cast<std::size_t>(inner_size)*static_cast<std::size_t>(outer)+
	  static_cast<std::size_t>(inner))+static_cast<std::size_t>(n));
}

std::size_t buffer_index_3d(int n, int nn, int outer, int middle, int inner,
			    int middle_size, int inner_size)
{
  return(static_cast<std::size_t>(nn)*
	 (static_cast<std::size_t>(inner_size)*
	  (static_cast<std::size_t>(middle_size)*static_cast<std::size_t>(outer)+
	   static_cast<std::size_t>(middle))+static_cast<std::size_t>(inner))+
	 static_cast<std::size_t>(n));
}

void resize_buffers(std::vector<double>& fold, std::vector<double>& fcpy, std::size_t count)
{
  if ((count > fold.max_size()/2) || (count > fcpy.max_size()/2)){
    mpi_abort_all(MPI_ERR_NO_MEM);
  }
  try{
    fold.resize(2*count);
    fcpy.resize(2*count);
  } catch (...){
    mpi_abort_all(MPI_ERR_NO_MEM);
  }
}

void check_mpi_error(int error_code)
{
  if (error_code != MPI_SUCCESS) mpi_abort_all(error_code);
}

}

int mpi_sync_status(int local_status, int *global_status)
{
  int status_min,status_max;
  int mpi_error;

  if (global_status == NULL) return(MPI_ERR_ARG);

  mpi_error=MPI_Allreduce(&local_status,&status_min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if (mpi_error != MPI_SUCCESS) return(mpi_error);
  mpi_error=MPI_Allreduce(&local_status,&status_max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if (mpi_error != MPI_SUCCESS) return(mpi_error);

  *global_status=status_min;
  if (status_min != status_max) return(MPI_ERR_OTHER);
  return(MPI_SUCCESS);
}

void mpi_sdrv2d(double *f[], int nn, int nx, int ny, int xoff, int yoff,
		int dnx, int dny,
		int mpi_rank, int mpi_numx, int mpi_numy)
/* MPI SendRecv for 2D variables */
/* Set dn=0 for Periodic boundary */
/* For other condition, call mpi_xbc2d and mpi_ybc2d later */
{
  int i,j,n;
  int mpi_tag=0;
  int rankl,rankh;
  std::size_t ntot;
  int mpi_count;
  MPI_Status r_stat;
  std::vector<double> fold,fcpy;

  /* XBC */
  if (dnx == 0){
    rankl=((mpi_rank % mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
    rankh=((mpi_rank % mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
  } else{
    rankl=((mpi_rank % mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
    rankh=((mpi_rank % mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
  }
  if (mpi_numx != 1){
    ntot=checked_mpi_count({nn,ny,xoff});
    mpi_count=static_cast<int>(ntot);
    resize_buffers(fold,fcpy,ntot);
    for (i=0;i<xoff;i++){
      for (j=0;j<ny;j++){	/* Transpose */
	for (n=0;n<nn;n++){
	  fold[buffer_index_2d(n,nn,i,j,ny)]=f[n][nx*j+(xoff+i)];
	  fold[buffer_index_2d(n,nn,2*xoff-1-i,j,ny)]=f[n][nx*j+(nx-xoff-1-i)];
	  fcpy[buffer_index_2d(n,nn,i,j,ny)]=f[n][nx*j+i];
	  fcpy[buffer_index_2d(n,nn,2*xoff-1-i,j,ny)]=f[n][nx*j+(nx-1-i)];
	}
      }
    }
    check_mpi_error(MPI_Sendrecv(fold.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 fcpy.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    check_mpi_error(MPI_Sendrecv(fold.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 fcpy.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	for (n=0;n<nn;n++){
	  f[n][nx*j+i]=fcpy[buffer_index_2d(n,nn,i,j,ny)];
	  f[n][nx*j+(nx-1-i)]=fcpy[buffer_index_2d(n,nn,2*xoff-1-i,j,ny)];
	}
      }
    }
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
    ntot=checked_mpi_count({nn,nx,yoff});
    mpi_count=static_cast<int>(ntot);
    resize_buffers(fold,fcpy,ntot);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (n=0;n<nn;n++){
	  fold[buffer_index_2d(n,nn,j,i,nx)]=f[n][nx*(yoff+j)+i];
	  fold[buffer_index_2d(n,nn,2*yoff-1-j,i,nx)]=f[n][nx*(ny-yoff-1-j)+i];
	  fcpy[buffer_index_2d(n,nn,j,i,nx)]=f[n][nx*j+i];
	  fcpy[buffer_index_2d(n,nn,2*yoff-1-j,i,nx)]=f[n][nx*(ny-1-j)+i];
	}
      }
    }
    check_mpi_error(MPI_Sendrecv(fold.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 fcpy.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    check_mpi_error(MPI_Sendrecv(fold.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 fcpy.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (n=0;n<nn;n++){
	  f[n][nx*j+i]=fcpy[buffer_index_2d(n,nn,j,i,nx)];
	  f[n][nx*(ny-1-j)+i]=fcpy[buffer_index_2d(n,nn,2*yoff-1-j,i,nx)];
	}
      }
    }
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
/* 2D X BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
{
  int i,j;
  if (std::abs(dn) == 1){
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
  } else if (std::abs(dn) == 2){
    /* Left */
    if ((mpi_rank % mpi_numx) == 0){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++) f[nx*j+i]=0.25*(2+dn)*f[nx*j+xoff];
      }
    }
    /* Right */
    if ((mpi_rank % mpi_numx) == (mpi_numx-1)){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff-st;i++) f[nx*j+(nx-1-i)]=0.25*(2+dn)*f[nx*j+(nx-1-xoff+st)];
      }
    }
  }
}

void mpi_ybc2d(double *f, int nx, int ny, int xoff, int yoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy)
/* 2D Y BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
{
  int i,j;
  if (std::abs(dn) == 1){
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
  } else if (std::abs(dn) == 2){
    /* Left */
    if (mpi_rank/mpi_numx == 0){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++) f[nx*j+i]=0.25*(2+dn)*f[nx*yoff+i];
      }
    }
    /* Right */
    if (mpi_rank/mpi_numx == (mpi_numy-1)){
      for (j=0;j<yoff-st;j++){
	for (i=0;i<nx;i++) f[nx*(ny-1-j)+i]=0.25*(2+dn)*f[nx*(ny-1-yoff+st)+i];
      }
    }
  }
}

void mpi_sdrv3d(double *f[], int nn, int nx, int ny, int nz, int xoff, int yoff, int zoff,
		int dnx, int dny, int dnz,
		int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* MPI SendRecv for 3D variables */
/* Set dn=0 for Periodic boundary */
/* For other condition, call mpi_x(y,z)bc3d later */
{
  int i,j,k,n;
  int m_xy=mpi_numx*mpi_numy;
  int mpi_tag=0;
  int rankl,rankh;
  std::size_t ntot;
  int mpi_count;
  MPI_Status r_stat;
  std::vector<double> fold,fcpy;

  /* XBC */
  if (dnx == 0){
    rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(mpi_rank+(mpi_numx-1)):(mpi_rank-1);
    rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(mpi_rank-(mpi_numx-1)):(mpi_rank+1);
  } else{
    rankl=(((mpi_rank%m_xy)%mpi_numx) == 0)?(MPI_PROC_NULL):(mpi_rank-1);
    rankh=(((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1))?(MPI_PROC_NULL):(mpi_rank+1);
  }
  if (mpi_numx != 1){
    ntot=checked_mpi_count({nn,ny,nz,xoff});
    mpi_count=static_cast<int>(ntot);
    resize_buffers(fold,fcpy,ntot);
    for (i=0;i<xoff;i++){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){	/* Transpose */
	  for (n=0;n<nn;n++){
	    fold[buffer_index_3d(n,nn,i,k,j,nz,ny)]=f[n][nx*(ny*k+j)+(xoff+i)];
	    fold[buffer_index_3d(n,nn,2*xoff-1-i,k,j,nz,ny)]=f[n][nx*(ny*k+j)+(nx-xoff-1-i)];
	    fcpy[buffer_index_3d(n,nn,i,k,j,nz,ny)]=f[n][nx*(ny*k+j)+i];
	    fcpy[buffer_index_3d(n,nn,2*xoff-1-i,k,j,nz,ny)]=f[n][nx*(ny*k+j)+(nx-1-i)];
	  }
	}
      }
    }
    check_mpi_error(MPI_Sendrecv(fold.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 fcpy.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    check_mpi_error(MPI_Sendrecv(fold.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 fcpy.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[buffer_index_3d(n,nn,i,k,j,nz,ny)];
	    f[n][nx*(ny*k+j)+(nx-1-i)]=fcpy[buffer_index_3d(n,nn,2*xoff-1-i,k,j,nz,ny)];
	  }
	}
      }
    }
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
    ntot=checked_mpi_count({nn,nz,nx,yoff});
    mpi_count=static_cast<int>(ntot);
    resize_buffers(fold,fcpy,ntot);
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	for (k=0;k<nz;k++){	/* Transpose */
	  for (n=0;n<nn;n++){
	    fold[buffer_index_3d(n,nn,j,i,k,nx,nz)]=f[n][nx*(ny*k+(yoff+j))+i];
	    fold[buffer_index_3d(n,nn,2*yoff-1-j,i,k,nx,nz)]=f[n][nx*(ny*k+(ny-yoff-1-j))+i];
	    fcpy[buffer_index_3d(n,nn,j,i,k,nx,nz)]=f[n][nx*(ny*k+j)+i];
	    fcpy[buffer_index_3d(n,nn,2*yoff-1-j,i,k,nx,nz)]=f[n][nx*(ny*k+(ny-1-j))+i];
	  }
	}
      }
    }
    check_mpi_error(MPI_Sendrecv(fold.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 fcpy.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    check_mpi_error(MPI_Sendrecv(fold.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 fcpy.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[buffer_index_3d(n,nn,j,i,k,nx,nz)];
	    f[n][nx*(ny*k+(ny-1-j))+i]=fcpy[buffer_index_3d(n,nn,2*yoff-1-j,i,k,nx,nz)];
	  }
	}
      }
    }
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
    ntot=checked_mpi_count({nn,nx,ny,zoff});
    mpi_count=static_cast<int>(ntot);
    resize_buffers(fold,fcpy,ntot);
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  for (n=0;n<nn;n++){
	    fold[buffer_index_3d(n,nn,k,j,i,ny,nx)]=f[n][nx*(ny*(zoff+k)+j)+i];
	    fold[buffer_index_3d(n,nn,2*zoff-1-k,j,i,ny,nx)]=f[n][nx*(ny*(nz-zoff-1-k)+j)+i];
	    fcpy[buffer_index_3d(n,nn,k,j,i,ny,nx)]=f[n][nx*(ny*k+j)+i];
	    fcpy[buffer_index_3d(n,nn,2*zoff-1-k,j,i,ny,nx)]=f[n][nx*(ny*(nz-1-k)+j)+i];
	  }
	}
      }
    }
    check_mpi_error(MPI_Sendrecv(fold.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 fcpy.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    check_mpi_error(MPI_Sendrecv(fold.data()+ntot,mpi_count,MPI_DOUBLE,rankh,mpi_tag,
				 fcpy.data(),mpi_count,MPI_DOUBLE,rankl,mpi_tag,
				 MPI_COMM_WORLD,&r_stat));
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  for (n=0;n<nn;n++){
	    f[n][nx*(ny*k+j)+i]=fcpy[buffer_index_3d(n,nn,k,j,i,ny,nx)];
	    f[n][nx*(ny*(nz-1-k)+j)+i]=fcpy[buffer_index_3d(n,nn,2*zoff-1-k,j,i,ny,nx)];
	  }
	}
      }
    }
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
/* 3D X BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
{
  int i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  if (std::abs(dn) == 1){
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
  } else if (std::abs(dn) == 2){
    /* Left */
    if (((mpi_rank%m_xy)%mpi_numx) == 0){
      for (k=0;k<nz;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<xoff;i++) f[nx*(ny*k+j)+i]=0.25*(2+dn)*f[nx*(ny*k+j)+xoff];
	}
      }
    }
    /* Right */
    if (((mpi_rank%m_xy)%mpi_numx) == (mpi_numx-1)){
      for (k=0;k<nz;k++){
        for (j=0;j<ny;j++){
	  for (i=0;i<xoff-st;i++) f[nx*(ny*k+j)+(nx-1-i)]=0.25*(2+dn)*f[nx*(ny*k+j)+(nx-1-xoff+st)];
        }
      }
    }
  }
}

void mpi_ybc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D Y BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
{
  int i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  if (std::abs(dn) == 1){
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
  } else if (std::abs(dn) == 2){
    /* Left */
    if (((mpi_rank%m_xy)/mpi_numx) == 0){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+j)+i]=0.25*(2+dn)*f[nx*(ny*k+yoff)+i];
	}
      }
    }
    /* Right */
    if (((mpi_rank%m_xy)/mpi_numx) == (mpi_numy-1)){
      for (k=0;k<nz;k++){
	for (j=0;j<yoff-st;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+(ny-1-j))+i]=0.25*(2+dn)*f[nx*(ny*k+(ny-1-yoff+st))+i];
	}
      }
    }
  }
}

void mpi_zbc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff, int st, int dn,
	       int mpi_rank, int mpi_numx, int mpi_numy, int mpi_numz)
/* 3D Z BC under MPI */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center) */
/* dn: Factor of Dirichlet (-1), Neumann (+1), Zero-fix (-2), Open (+2). if dn==0, nothing to do */
{
  int i,j,k;
  int m_xy=mpi_numx*mpi_numy;
  if (std::abs(dn) == 1){
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
  } else if (std::abs(dn) == 2){
    /* Left */
    if ((mpi_rank/m_xy) == 0){
      for (k=0;k<zoff;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*k+j)+i]=0.25*(2+dn)*f[nx*(ny*zoff+j)+i];
	}
      }
    }
    /* Right */
    if ((mpi_rank/m_xy) == (mpi_numz-1)){
      for (k=0;k<zoff-st;k++){
	for (j=0;j<ny;j++){
	  for (i=0;i<nx;i++) f[nx*(ny*(nz-1-k)+j)+i]=0.25*(2+dn)*f[nx*(ny*(nz-1-zoff+st)+j)+i];
	}
      }
    }
  }
}
