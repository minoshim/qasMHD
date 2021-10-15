void dataio(double *p[], int nm, int nd, int n, int cnt, double tim, int msg, int mpi_rank)
// Output simulation data
{
  int i;
  int m_xy=mpi_numx*mpi_numy;
  char filname[100];
  FILE *outfil;
  float *val;
  val=new float[nd];

  if (n == 0){
    if (mpi_rank == 0){
      sprintf(filname,"%s/params.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%f\n",gam);
      fclose(outfil);

      sprintf(filname,"%s/offsets.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%d %d %d\n",xoff,yoff,zoff);
      fclose(outfil);

      sprintf(filname,"%s/t.dat",fildir);
      outfil=fopen(filname,"w");
      fprintf(outfil,"%f\n",tim);
      fclose(outfil);
    }
    if ((mpi_rank/m_xy == 0) && (((mpi_rank%m_xy)/mpi_numx) == 0)){
      sprintf(filname,"%s/x_%05d.dat",fildir,((mpi_rank%m_xy)%mpi_numx));
      outfil=fopen(filname,"w");
      for (i=0;i<nx;i++) fprintf(outfil,"%.12f\n",x[i]);
      fclose(outfil);
    }
    if ((mpi_rank/m_xy == 0) && (((mpi_rank%m_xy)%mpi_numx) == 0)){
      sprintf(filname,"%s/y_%05d.dat",fildir,((mpi_rank%m_xy)/mpi_numx));
      outfil=fopen(filname,"w");
      for (i=0;i<ny;i++) fprintf(outfil,"%.12f\n",y[i]);
      fclose(outfil);
    }
    if ((mpi_rank%m_xy) == 0){
      sprintf(filname,"%s/z_%05d.dat",fildir,(mpi_rank/m_xy));
      outfil=fopen(filname,"w");
      for (i=0;i<nz;i++) fprintf(outfil,"%.12f\n",z[i]);
      fclose(outfil);
    }
  } else{
    if (mpi_rank == 0){
      sprintf(filname,"%s/t.dat",fildir);
      outfil=fopen(filname,"a");
      fprintf(outfil,"%f\n",tim);
      fclose(outfil);
    }
  }
  sprintf(filname,"%s/outdat_%05d_%05d.dat",fildir,cnt,mpi_rank);
  outfil=fopen(filname,"wb");
  for (i=0;i<nm;i++){
    conv_d2f(val,p[i],nd);
    fwrite(val,sizeof(float),nd,outfil);
  }
  fclose(outfil);
  if (msg){
    printf("Output data at t=%.4f (%d iterations).\n",tim,n);
  }

  delete[] val;
}
