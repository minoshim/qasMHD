void dataio(double *p[], int nm, int nd, int n, int cnt, double tim, int msg)
// Output simulation data
{
  int i;
  char filname[100];
  FILE *outfil;

  if (n == 0){
    sprintf(filname,"%s/params.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%f\n",gam);
    fclose(outfil);

    sprintf(filname,"%s/offsets.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%d %d\n",xoff,yoff);
    fclose(outfil);
    
    sprintf(filname,"%s/t.dat",fildir);
    outfil=fopen(filname,"w");
    fprintf(outfil,"%f\n",tim);
    fclose(outfil);
  } else{
    sprintf(filname,"%s/t.dat",fildir);
    outfil=fopen(filname,"a");
    fprintf(outfil,"%f\n",tim);
    fclose(outfil);
  }
  sprintf(filname,"%s/outdat_%05d.dat",fildir,cnt);  
  outfil=fopen(filname,"wb");
  for (i=0;i<nm;i++){
    fwrite(p[i],sizeof(double),nd,outfil);
  }
  fclose(outfil);
  if (msg){
    printf("Output data at t=%.4f (%d iterations).\n",tim,n);
  }

}
