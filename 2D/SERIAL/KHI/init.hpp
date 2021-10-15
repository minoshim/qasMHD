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
  // KHI
  int i,j,ss;
  double dvy[nx];
  unsigned seed;

  seed=(unsigned)time(NULL);

  for (i=0;i<nx;i++){
    dvy[i]=0.0;
#if (RANDOM)
    double dvpara[2]={0,dv};
    dvy[i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#else
    dvy[i]+=dv*sin(2*pi*x[i]/wlen); // Single mode perturbation
#endif    
  }

  for (j=0;j<ny;j++){
    double angle=0.5*((angle_u+angle_l)+(angle_u-angle_l)*tanh((y[j]-s0)/lambda));
    for (i=0;i<nx;i++){
      double vx,vy,vz,cx,cy,cz,pr;
      double xm,ym;

      ss=nx*j+i;
      xm=x[i]-0.5*dx;
      ym=y[j]-0.5*dy;

      ro[ss]=0.5*((ro_u+ro_l)+(ro_u-ro_l)*tanh((y[j]-s0)/lambda));

      vx=+vamp*tanh((y[j]-s0)/lambda);
      vy=dvy[i]*exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda));
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
      pr=0.5*beta*(cx*cx+cy*cy+cz*cz);

      mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);      
    }
  }
}
