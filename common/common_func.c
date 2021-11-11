#include "common_func.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double max(double a, double b)
{
  return( (a >= b)?a:b );
}

double min(double a, double b)
{
  return( (a >= b)?b:a );
}

double minmod(double a, double b)
{
  return(+min(max(a,b),0)
	 +max(min(a,b),0));
}

double minmod3(double a, double b, double c)
{
  return(+min(max(max(a,b),c),0)
	 +max(min(min(a,b),c),0));
}

double rand_noise(const double *params, unsigned seed)
{
  static int r_flag=0;
  if (r_flag == 0){
    srandom(seed);
    r_flag=1;
  }
  return(params[0]+params[1]*((double)random()/RAND_MAX-0.5)*2.0);
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

void rk_updt(double *f1, double f0, double df, double rkfac0, double rkfac1)
/* TVD Runge-Kutta update */
{
  (*f1)=rkfac0*f0+rkfac1*((*f1)+df);
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

void bkup_load(double *p[], int nm, int nd, int *n, int *cnt, double *tim, double *dt, double *trec, int mpi_rank, const char *fildir)
// Load bkup files
{
  int i;
  FILE *infil;
  int vali[2];
  double vald[3];
  char filname[100];
  int itmp;
  size_t ttmp;
  
  sprintf(filname,"%s/bkup_stamp.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    itmp=fscanf(infil,"%d %d\n",&vali[0],&vali[1]);
    fclose(infil);
    *n=vali[0];
    *cnt=vali[1];

    sprintf(filname,"%s/bkup_time.dat",fildir);
    infil=fopen(filname,"rb");
    ttmp=fread(&vald,sizeof(double),3,infil);
    fclose(infil);
    *tim=vald[0];
    *dt=vald[1];
    *trec=vald[2];

    sprintf(filname,"%s/bkup_data_%05d.dat",fildir,mpi_rank);
    infil=fopen(filname,"rb");
    for (i=0;i<nm;i++){
      ttmp=fread(p[i],sizeof(double),nd,infil);
    }
    fclose(infil);
    
    if (mpi_rank == 0) printf("Load backup files at %d steps (T = %f)\n",vali[0],vald[0]);
  }

}

void bkup_save(double *p[], int nm, int nd, int n, int cnt, double tim, double dt, double trec, int mpi_rank, const char *fildir)
// Save bkup files
{
  int i;
  FILE *outfil;
  double vald[]={tim,dt,trec};
  char filname[100];

  if (mpi_rank == 0){
    sprintf(filname,"%s/bkup_stamp.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%d %d\n",n,cnt);
    fclose(outfil);

    sprintf(filname,"%s/bkup_time.dat",fildir);
    outfil=fopen(filname,"wb");
    fwrite(&vald,sizeof(double),3,outfil);
    fclose(outfil);
  }

  sprintf(filname,"%s/bkup_data_%05d.dat",fildir,mpi_rank);
  outfil=fopen(filname,"wb");
  for (i=0;i<nm;i++){
    fwrite(p[i],sizeof(double),nd,outfil);
  }
  fclose(outfil);
}
