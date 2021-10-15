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
    // x[i]=(i-xoff+0.5)*dx+xmin;
    x[i]=(i-xoff)*dx+xmin;
    fprintf(outfil,"%.12f\n",x[i]);
  }
  fclose(outfil);

  sprintf(filname,"%s/y.dat",fildir);
  outfil=fopen(filname,"w");
  for (j=0;j<ny;j++){
    // y[j]=(j-yoff+0.5)*dy+ymin;
    y[j]=(j-yoff)*dy+ymin;
    fprintf(outfil,"%.12f\n",y[j]);
  }
  fclose(outfil);

}

void init_plasma()
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
