void get_dt(double *dt, double cfl, double dr, int comment);
double get_dcoefs(double *nu, double *eta);

void get_dt(double *dt, double cfl, double dr, int comment)
{
  // Modify time step
  int i,j,k,ss;
  double vmax=0.0;

  for (k=zoff;k<nz-zoff;k++){
    double vx,vy,vz,cx,cy,cz,pr,cf,vtmp;
    for (j=yoff;j<ny-yoff;j++){
      for (i=xoff;i<nx-xoff;i++){
	ss=nx*(ny*k+j)+i;
	cx=0.5*(bx[ss]+bx[nx*(ny*k+j)+(i+1)]); // 3D CT
	cy=0.5*(by[ss]+by[nx*(ny*k+(j+1))+i]); // 3D CT
	cz=0.5*(bz[ss]+bz[nx*(ny*(k+1)+j)+i]); // 3D CT
	mhd_prmtv(ro[ss],mx[ss],my[ss],mz[ss],cx,cy,cz,en[ss],
		  &vx,&vy,&vz,&pr,gam);
	cf=sqrt((gam*pr+(cx*cx+cy*cy+cz*cz))/ro[ss]);
	vtmp=sqrt(vx*vx+vy*vy+vz*vz)+cf; // 3D
	if (vtmp > vmax) vmax=vtmp;
      }
    }
  }

  // MPI Allreduce
  double vmax_a;
  MPI_Allreduce(&vmax,&vmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  vmax=vmax_a;

  if (finite(vmax)) (*dt)=cfl*dr/vmax;
  if (comment){
    printf("vmax = %f, dt = %.9f, dx = %f\n",vmax,(*dt),dr);
    printf("CFL number = %f\n",vmax*(*dt)/dr);
    printf("CFL condition for diffusion is %f\n ",max(eta0,2*nu0)*(*dt)/(dr*dr));
    // Factor 2 is multiplied in viscous coef. for robust estimation
  }
}

double get_dcoefs(double *nu, double *eta)
{
  // Get diffusion coefficients
  int i,j,k,ss;
  double dcmax=0.0,dtmp;
  for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      for (i=0;i<nx;i++){
	ss=nx*(ny*k+j)+i;
	nu[ss]=nu0;
	eta[ss]=eta0;

	// Localized diffusion
	// double r=sqrt(x[i]*x[i]+y[j]*y[j]);
	// double r=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
	// nu[ss]=nu0*exp(-r/1.0);
	// eta[ss]=eta0*exp(-r/1.0);
      
	dtmp=max(2*nu[ss],eta[ss]); // Factor 2 is multiplied in viscous coef. for robust estimation
	if (dtmp > dcmax) dcmax=dtmp;
      }
    }
  }

  // MPI Allreduce
  double dcmax_a;
  MPI_Allreduce(&dcmax,&dcmax_a,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  dcmax=dcmax_a;
  
  return dcmax;
}
