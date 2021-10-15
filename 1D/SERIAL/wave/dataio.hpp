void dataio(int n, int cnt, double tim)
// Output simulation data
{
  int i;
  char filname[100];
  FILE *outfil;

  if (n == 0){
    sprintf(filname,"%s/params.dat",fildir);  
    outfil=fopen(filname,"w");
    fprintf(outfil,"%.12f %.12f\n",bx,gam);
    fclose(outfil);

    sprintf(filname,"%s/xoff.dat",fildir);  
    outfil=fopen(filname,"w");
    fprintf(outfil,"%d\n",xoff);
    fclose(outfil);

    sprintf(filname,"%s/t.dat",fildir);  
    outfil=fopen(filname,"w");
    fprintf(outfil,"%.12f\n",tim);
    fclose(outfil);
  } else{
    sprintf(filname,"%s/t.dat",fildir);  
    outfil=fopen(filname,"a");
    fprintf(outfil,"%.12f\n",tim);
    fclose(outfil);
  }
  sprintf(filname,"%s/outdat_%05d.dat",fildir,cnt);  
  outfil=fopen(filname,"wb");
  double *p[]={ro,mx,my,mz,by,bz,en};
  for (i=0;i<sizeof(p)/sizeof(p[0]);i++){
    fwrite(p[i],sizeof(double),nx,outfil);
  }
  fclose(outfil);

}
