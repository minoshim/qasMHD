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
  // Blast wave
  int i,j;
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      double rr,vx,vy,vz,pr;
      rr=sqrt(x[i]*x[i]+y[j]*y[j]);
      
      ro[ss]=(rr <= r_0)?ro1:ro0;
      vx=0.0;
      vy=0.0;
      vz=0.0;
      pr=(rr <= r_0)?pr1:pr0;
      
      bx[ss]=b_0*cos(ban*dtor);
      by[ss]=b_0*sin(ban*dtor);
      bz[ss]=0.0;

      mhd_cnsvt(ro[ss],vx,vy,vz,bx[ss],by[ss],bz[ss],pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
    }
  }
}
