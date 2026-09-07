#include "common_func.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <string>

double rand_noise(const double *params, unsigned seed)
{
  static int r_flag=0;
  if (r_flag == 0){
    srandom(seed);
    r_flag=1;
  }
  return(params[0]+params[1]*((double)random()/RAND_MAX-0.5)*2.0);
}

double rand_noise_mt(const double *params, std::mt19937 &engine)
{
  const double range=static_cast<double>(std::mt19937::max())+1.0;
  const double unit=(static_cast<double>(engine())+0.5)/range;
  return(params[0]+params[1]*(2.0*unit-1.0));
}

void cpy_array(double *a, const double *b, int n)
{
  int i;
  for (i=0;i<n;i++) a[i]=b[i];
}

void conv_d2f(float *valo, const double *vali, int n)
{
  int i;
  for (i=0;i<n;i++) valo[i]=(float)vali[i];
}

void conv_f2d(double *valo, const float *vali, int n)
{
  int i;
  for (i=0;i<n;i++) valo[i]=(double)vali[i];
}

void bc1d(double *f, int nx, int xoff, int dnx)
/* 1D boundary condition */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann, -2 for zero-fix, +2 for open condition */
{
  int i;
  if (dnx == 0){
    /* Periodic */
    for (i=0;i<xoff;i++){
      f[i]=f[nx-2*xoff+i];      
      f[nx-1-i]=f[2*xoff-1-i];
    }
  } else if (abs(dnx) == 1){
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (i=0;i<xoff;i++){
      f[i]=dnx*f[2*xoff-1-i];
      f[nx-1-i]=dnx*f[nx-2*xoff+i];
    }
  } else if (abs(dnx) == 2){
    /* Zero-fix (dn = -2) or Open (dn = +2) */
    for (i=0;i<xoff;i++){
      f[i]=0.25*(2+dnx)*f[xoff];
      f[nx-1-i]=0.25*(2+dnx)*f[nx-1-xoff];
    }
  }
}

void bc2d(double *f, int nx, int ny, int xoff, int yoff,
	  int stx, int dnx, int sty, int dny)
/* 2D Boundary condition */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center), else 0 */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann, -2 for zero-fix, +2 for open condition */
{
  int i,j;
  if (dnx == 0){
    /* Periodic */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f[nx*j+(nx-1-i)]=f[nx*j+(2*xoff-1-i)];
	f[nx*j+i]=f[nx*j+(nx-2*xoff+i)];
      }
    }
  } else if (abs(dnx) == 1){
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	/* Left */
	f[nx*j+i]=dnx*f[nx*j+(2*xoff-1+stx)-i];
      }
      for (i=0;i<xoff-stx;i++){
	/* Right */
	f[nx*j+(nx-1-i)]=dnx*f[nx*j+(nx-2*xoff+stx)+i];
      }
    }
  } else if (abs(dnx) == 2){
    /* Zero-fix (dn = -2) or Open (dn = +2) */
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	/* Left */
	f[nx*j+i]=0.25*(2+dnx)*f[nx*j+xoff];
      }
      for (i=0;i<xoff-stx;i++){
	/* Right */
	f[nx*j+(nx-1-i)]=0.25*(2+dnx)*f[nx*j+(nx-1-xoff+stx)];
      }
    }
  }
  
  if (dny == 0){
    /* Periodic */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*(ny-1-j)+i]=f[nx*(2*yoff-1-j)+i];
	f[nx*j+i]=f[nx*(ny-2*yoff+j)+i];
      }
    }
  } else if (abs(dny) == 1){
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	/* Left */
	f[nx*j+i]=dny*f[nx*(2*yoff-1+sty-j)+i];	
      }
    }
    for (j=0;j<yoff-sty;j++){
      for (i=0;i<nx;i++){
	/* Right */
	f[nx*(ny-1-j)+i]=dny*f[nx*(ny-2*yoff+sty+j)+i];
      }
    }
  } else if (abs(dny) == 2){
    /* Zero-fix (dn = -2) or Open (dn = +2) */
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	/* Left */
	f[nx*j+i]=0.25*(2+dny)*f[nx*yoff+i];
      }
    }
    for (j=0;j<yoff-sty;j++){
      for (i=0;i<nx;i++){
	/* Right */
	f[nx*(ny-1-j)+i]=0.25*(2+dny)*f[nx*(ny-1-yoff+sty)+i];
      }
    }
  }
}

void bc3d(double *f, int nx, int ny, int nz, int xoff, int yoff, int zoff,
	  int stx, int dnx, int sty, int dny, int stz, int dnz)
/* 3D Boundary condition */
/* st: Flag for staggered grid. Set 1 when f is @ cell face (not center), else 0 */
/* dn: 0 for periodic, -1 for Dirichlet, +1 for Neumann, -2 for zero-fix, +2 for open condition */
{
  int i,j,k;
  if (dnx == 0){
    /* Periodic */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  f[nx*(ny*k+j)+(nx-1-i)]=f[nx*(ny*k+j)+( 2*xoff-1-i)];
	  f[nx*(ny*k+j)+(     i)]=f[nx*(ny*k+j)+(nx-2*xoff+i)];
	}
      }
    }
  } else if (abs(dnx) == 1){
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  /* Left */
	  f[nx*(ny*k+j)+(     i)]=dnx*f[nx*(ny*k+j)+( 2*xoff-1+stx-i)];
	}
	for (i=0;i<xoff-stx;i++){
	  /* Right */
	  f[nx*(ny*k+j)+(nx-1-i)]=dnx*f[nx*(ny*k+j)+(nx-2*xoff+stx+i)];
	}
      }
    }
  } else if (abs(dnx) == 2){
    /* Zero-fix (dn = -2) or Open (dn = +2) */
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<xoff;i++){
	  /* Left */
	  f[nx*(ny*k+j)+(     i)]=0.25*(2+dnx)*f[nx*(ny*k+j)+(         xoff)];
	}
	for (i=0;i<xoff-stx;i++){
	  /* Right */
	  f[nx*(ny*k+j)+(nx-1-i)]=0.25*(2+dnx)*f[nx*(ny*k+j)+(nx-1-xoff+stx)];
	}
      }
    }
  }
  
  if (dny == 0){
    /* Periodic */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*k+(ny-1-j))+i]=f[nx*(ny*k+( 2*yoff-1-j))+i];
	  f[nx*(ny*k+(     j))+i]=f[nx*(ny*k+(ny-2*yoff+j))+i];
	}
      }
    }
  } else if (abs(dny) == 1){
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  /* Left */
	  f[nx*(ny*k+(     j))+i]=dny*f[nx*(ny*k+( 2*yoff-1+sty-j))+i];	
	}
      }
      for (j=0;j<yoff-sty;j++){
	for (i=0;i<nx;i++){
	  /* Right */
	  f[nx*(ny*k+(ny-1-j))+i]=dny*f[nx*(ny*k+(ny-2*yoff+sty+j))+i];
	}
      }
    }
  } else if (abs(dny) == 2){
    /* Zero-fix (dn = -2) or Open (dn = +2) */
    for (k=0;k<nz;k++){
      for (j=0;j<yoff;j++){
	for (i=0;i<nx;i++){
	  /* Left */
	  f[nx*(ny*k+(     j))+i]=0.25*(2+dny)*f[nx*(ny*k+(         yoff))+i];
	}
      }
      for (j=0;j<yoff-sty;j++){
	for (i=0;i<nx;i++){
	  /* Right */
	  f[nx*(ny*k+(ny-1-j))+i]=0.25*(2+dny)*f[nx*(ny*k+(ny-1-yoff+sty))+i];
	}
      }
    }
  }

  if (dnz == 0){
    /* Periodic */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  f[nx*(ny*(nz-1-k)+j)+i]=f[nx*(ny*( 2*zoff-1-k)+j)+i];
	  f[nx*(ny*(     k)+j)+i]=f[nx*(ny*(nz-2*zoff+k)+j)+i];
	}
      }
    }
  } else if (abs(dnz) == 1){
    /* Dirichlet (dn = -1) or Neumann (dn = +1) */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  /* Left */
	  f[nx*(ny*(     k)+j)+i]=dnz*f[nx*(ny*( 2*zoff-1+stz-k)+j)+i];
	}
      }
    }
    for (k=0;k<zoff-stz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  /* Right */
	  f[nx*(ny*(nz-1-k)+j)+i]=dnz*f[nx*(ny*(nz-2*zoff+stz+k)+j)+i];	  
	}
      }
    }
  } else if (abs(dnz) == 2){
    /* Zero-fix (dn = -2) or Open (dn = +2) */
    for (k=0;k<zoff;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  /* Left */
	  f[nx*(ny*(     k)+j)+i]=0.25*(2+dnz)*f[nx*(ny*(         zoff)+j)+i];
	}
      }
    }
    for (k=0;k<zoff-stz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  /* Right */
	  f[nx*(ny*(nz-1-k)+j)+i]=0.25*(2+dnz)*f[nx*(ny*(nz-1-zoff+stz)+j)+i];
	}
      }
    }
  }
}

static std::string bkup_path(const char *fildir, const char *filename)
{
  std::string path(fildir);
  if (!path.empty() && path[path.size()-1] != '/') path += '/';
  path += filename;
  return(path);
}

static int bkup_check_array(double *p[], int nm)
{
  int i;
  if (p == NULL) return(0);
  for (i=0;i<nm;i++){
    if (p[i] == NULL) return(0);
  }
  return(1);
}

static int bkup_close(FILE *file, const std::string &path)
{
  if (fclose(file) != 0){
    fprintf(stderr,"Cannot close backup file %s: %s\n",path.c_str(),strerror(errno));
    return(BKUP_ERROR);
  }
  return(BKUP_OK);
}

int bkup_load(double *p[], int nm, int nd, int *n, int *cnt,
              double *tim, double *dt, double *trec,
              int mpi_rank, const char *fildir)
// Load bkup files
{
  int i;
  FILE *infil;
  int vali[2];
  double vald[3];
  char filename[32];
  char extra;
  int itmp;
  size_t ttmp;
  std::string path;

  if (nm <= 0 || nd <= 0 || mpi_rank < 0 || fildir == NULL || fildir[0] == '\0' ||
      n == NULL || cnt == NULL || tim == NULL || dt == NULL || trec == NULL ||
      !bkup_check_array(p,nm)){
    fprintf(stderr,"Invalid argument passed to bkup_load\n");
    return(BKUP_ERROR);
  }

  path=bkup_path(fildir,"bkup_stamp.dat");
  infil=fopen(path.c_str(),"r");
  if (infil == NULL){
    if (errno == ENOENT) return(BKUP_NOT_FOUND);
    fprintf(stderr,"Cannot open backup file %s: %s\n",path.c_str(),strerror(errno));
    return(BKUP_ERROR);
  }
  itmp=fscanf(infil,"%d %d %c",&vali[0],&vali[1],&extra);
  if (itmp != 2){
    fprintf(stderr,"Invalid backup stamp in %s\n",path.c_str());
    bkup_close(infil,path);
    return(BKUP_ERROR);
  }
  if (bkup_close(infil,path) != BKUP_OK) return(BKUP_ERROR);

  path=bkup_path(fildir,"bkup_time.dat");
  infil=fopen(path.c_str(),"rb");
  if (infil == NULL){
    fprintf(stderr,"Cannot open backup file %s: %s\n",path.c_str(),strerror(errno));
    return(BKUP_ERROR);
  }
  ttmp=fread(vald,sizeof(double),3,infil);
  if (ttmp != 3 || fgetc(infil) != EOF || ferror(infil)){
    fprintf(stderr,"Invalid backup time data in %s (expected 3 values, read %zu)\n",
            path.c_str(),ttmp);
    bkup_close(infil,path);
    return(BKUP_ERROR);
  }
  if (bkup_close(infil,path) != BKUP_OK) return(BKUP_ERROR);

  itmp=snprintf(filename,sizeof(filename),"bkup_data_%05d.dat",mpi_rank);
  if (itmp < 0 || (size_t)itmp >= sizeof(filename)){
    fprintf(stderr,"Cannot construct backup filename for rank %d\n",mpi_rank);
    return(BKUP_ERROR);
  }
  path=bkup_path(fildir,filename);
  infil=fopen(path.c_str(),"rb");
  if (infil == NULL){
    fprintf(stderr,"Cannot open backup file %s: %s\n",path.c_str(),strerror(errno));
    return(BKUP_ERROR);
  }
  for (i=0;i<nm;i++){
    ttmp=fread(p[i],sizeof(double),nd,infil);
    if (ttmp != (size_t)nd){
      fprintf(stderr,"Invalid backup data in %s for variable %d (expected %d values, read %zu)\n",
              path.c_str(),i,nd,ttmp);
      bkup_close(infil,path);
      return(BKUP_ERROR);
    }
  }
  if (fgetc(infil) != EOF || ferror(infil)){
    fprintf(stderr,"Invalid backup data size in %s\n",path.c_str());
    bkup_close(infil,path);
    return(BKUP_ERROR);
  }
  if (bkup_close(infil,path) != BKUP_OK) return(BKUP_ERROR);

  *n=vali[0];
  *cnt=vali[1];
  *tim=vald[0];
  *dt=vald[1];
  *trec=vald[2];

  if (mpi_rank == 0) printf("Load backup files at %d steps (T = %f)\n",vali[0],vald[0]);
  return(BKUP_OK);
}

int bkup_save(double *p[], int nm, int nd, int n, int cnt,
              double tim, double dt, double trec,
              int mpi_rank, const char *fildir)
// Save bkup files
{
  int i;
  FILE *outfil;
  double vald[]={tim,dt,trec};
  char filename[32];
  int itmp;
  size_t ttmp;
  std::string path;

  if (nm <= 0 || nd <= 0 || mpi_rank < 0 || fildir == NULL || fildir[0] == '\0' ||
      !bkup_check_array(p,nm)){
    fprintf(stderr,"Invalid argument passed to bkup_save\n");
    return(BKUP_ERROR);
  }

  if (mpi_rank == 0){
    path=bkup_path(fildir,"bkup_stamp.dat");
    outfil=fopen(path.c_str(),"w");
    if (outfil == NULL){
      fprintf(stderr,"Cannot open backup file %s: %s\n",path.c_str(),strerror(errno));
      return(BKUP_ERROR);
    }
    if (fprintf(outfil,"%d %d\n",n,cnt) < 0){
      fprintf(stderr,"Cannot write backup file %s: %s\n",path.c_str(),strerror(errno));
      bkup_close(outfil,path);
      return(BKUP_ERROR);
    }
    if (bkup_close(outfil,path) != BKUP_OK) return(BKUP_ERROR);

    path=bkup_path(fildir,"bkup_time.dat");
    outfil=fopen(path.c_str(),"wb");
    if (outfil == NULL){
      fprintf(stderr,"Cannot open backup file %s: %s\n",path.c_str(),strerror(errno));
      return(BKUP_ERROR);
    }
    ttmp=fwrite(vald,sizeof(double),3,outfil);
    if (ttmp != 3){
      fprintf(stderr,"Cannot write backup file %s (expected 3 values, wrote %zu)\n",
              path.c_str(),ttmp);
      bkup_close(outfil,path);
      return(BKUP_ERROR);
    }
    if (bkup_close(outfil,path) != BKUP_OK) return(BKUP_ERROR);
  }

  itmp=snprintf(filename,sizeof(filename),"bkup_data_%05d.dat",mpi_rank);
  if (itmp < 0 || (size_t)itmp >= sizeof(filename)){
    fprintf(stderr,"Cannot construct backup filename for rank %d\n",mpi_rank);
    return(BKUP_ERROR);
  }
  path=bkup_path(fildir,filename);
  outfil=fopen(path.c_str(),"wb");
  if (outfil == NULL){
    fprintf(stderr,"Cannot open backup file %s: %s\n",path.c_str(),strerror(errno));
    return(BKUP_ERROR);
  }
  for (i=0;i<nm;i++){
    ttmp=fwrite(p[i],sizeof(double),nd,outfil);
    if (ttmp != (size_t)nd){
      fprintf(stderr,"Cannot write backup file %s for variable %d (expected %d values, wrote %zu)\n",
              path.c_str(),i,nd,ttmp);
      bkup_close(outfil,path);
      return(BKUP_ERROR);
    }
  }
  if (bkup_close(outfil,path) != BKUP_OK) return(BKUP_ERROR);
  return(BKUP_OK);
}
