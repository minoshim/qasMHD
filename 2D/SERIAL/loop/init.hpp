void init_grid();
void init_plasma();

void init_grid()
// Define X and Y coordinates
{
  int i,j;
  char filname[100];
  FILE *outfil;

  sprintf(filname,"%s/x.dat",fildir);
  outfil=fopen(filname,"w");
  for (i=0;i<nx;i++){
    x[i]=(i-xoff+0.5)*dx+xmin;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (j=0;j<ny;j++){
    y[j]=(j-yoff+0.5)*dy+ymin;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}

void init_plasma()
// Set initial condition
{
  // Magnetic loop advection in high beta plasma
  int i,j;
  double *azp,*azc;
  azp=new double[nd];
  azc=new double[nd];
  // Vector potential
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double r,x0,y0;
      x0=x[i]-0.5*dx;
      y0=y[j]-0.5*dy;
      r=sqrt(x0*x0+y0*y0);
      azp[ss]=(r <= rad)?a0*(rad-r):0;
      x0=x[i];
      y0=y[j];
      r=sqrt(x0*x0+y0*y0);
      azc[ss]=(r <= rad)?a0*(rad-r):0;
    }
  }
  // Conservative variables
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double vx,vy,vz,pr,cx,cy;
      ro[ss]=ro0;
      vx=vx0;
      vy=vy0;
      vz=vz0;
      pr=pr0;

      bx[ss]=by[ss]=bz[ss]=0.0;
      cx=cy=0.0;
      if (j >= 1 && j <= ny-2){
	bx[ss]=+(azp[nx*(j+1)+i]-azp[nx*j+i])*idy;
	cx=+(azc[nx*(j+1)+i]-azc[nx*(j-1)+i])*idy*0.5;
      }
      if (i >= 1 && i <= nx-2){
	by[ss]=-(azp[nx*j+(i+1)]-azp[nx*j+i])*idx;
	cy=-(azc[nx*j+(i+1)]-azc[nx*j+(i-1)])*idx*0.5;
      }
      
      mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,bz[ss],pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
    }
  }

  delete[] azp;
  delete[] azc;
}
