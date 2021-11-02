void init_grid();
void init_plasma();
double cal_pressure(double y0, double y1, double pr0, int n);

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
  // RTI
  int i,j,ss;
  double dvy[nx];
  unsigned seed;
  double prmin=1e-4;

  seed=(unsigned)time(NULL);

  for (i=0;i<nx;i++){
    dvy[i]=0.0;
#if (RANDOM)
    double dvpara[2]={0,dv};
    dvy[i]+=rand_noise(dvpara,seed); // Multiple mode perturbation
#else
    dvy[i]+=dv*cos(2*pi*x[i]/wlen); // Single mode perturbation
#endif    
  }

  for (j=0;j<ny;j++){
    double angle=sqrwave2(angle_u,angle_l,(y[j]-s0)/lambda,(y[j]+s0)/lambda);
    for (i=0;i<nx;i++){
      double vx,vy,vz,cx,cy,cz,pr;
      double xm,ym;

      ss=nx*j+i;
      xm=x[i]-0.5*dx;
      ym=y[j]-0.5*dy;

      ro[ss]=sqrwave2(ro_u,ro_l,(y[j]-s0)/lambda,(y[j]+s0)/lambda);
      vx=sqrwave2(vamp,-vamp,(y[j]-s0)/lambda,(y[j]+s0)/lambda);
      vy=dvy[i]*(exp(-((y[j]-s0)*(y[j]-s0))/(4*lambda*lambda))-exp(-((y[j]+s0)*(y[j]+s0))/(4*lambda*lambda)));
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

      // Gravitational potential, phi_g = g0*lg*log(cosh(y/lg))
      // d(phi_g)/dy = g0*tanh(y/lg)
      phi_g[ss]=g_potential(y[j],g0,lg,0);
      
      // Pressure is obtained by integrating the equilibrium
      pr=cal_pressure(0,y[j],pr0,1+(int)(fabs(y[j])/dy+0.5));
      pr=max(pr,prmin);
      
      mhd_cnsvt(ro[ss],vx,vy,vz,cx,cy,cz,pr,&mx[ss],&my[ss],&mz[ss],&en[ss],gam);
      en[ss]+=ro[ss]*phi_g[ss];	// Add gravitational potential
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
