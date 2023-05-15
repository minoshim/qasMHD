#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
/* Merge 2D MPI-decomposed data */
/* > ./merge.out FILE_DIRECTORY OUTPUT_DIRECTORY */
{
  if (argc !=3){
    fputs("Set file directory and output directory\n",stderr);
    fputs("Command Syntax: ./a.out FILDIR OUTDIR\n",stderr);
    return 0;
  }
  char *chtmp;
  char fildir[100];
  char outdir[100];
  chtmp=argv[1];
  strcpy(fildir,chtmp);
  strcat(fildir,"/");
  chtmp=argv[2];
  strcpy(outdir,chtmp);
  strcat(outdir,"/");
  printf("File directory is %s\n",fildir);
  printf("Output directory is %s\n",outdir);
  puts("");

  const int nstt=0;		/* Start time index */
  const int smax=8;		/* Number of MHD variables */
  int i,j,m,n,mx,my,s;
  int nt=0;
  int nxall=0,nyall=0;
  int xoff=0,yoff=0;
  int mpix=1,mpiy=1;
  int isum=0,jsum=0;
  int itmp,flag;

  FILE *infil,*outfil;
  char filname[300];
  char command[300];
  double dummy;

  /* Read offset data */
  sprintf(filname,"%soffsets.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    int tmp[2]={0};
    itmp=fscanf(infil,"%d %d\n",&tmp[0],&tmp[1]);
    xoff=tmp[0];
    yoff=tmp[1];
    fclose(infil);

    puts("Offset data is found.");
    printf("Offsets in X and Y directions are %d and %d.\n",xoff,yoff);
    puts("");
  } else{
    puts("No offset data is found.");
    return 0;
  }

  /* Read MPI number data */
  sprintf(filname,"%smpinum.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    int tmp[2]={0};
    itmp=fscanf(infil,"%d %d\n",&tmp[0],&tmp[1]);
    mpix=tmp[0];
    mpiy=tmp[1];
    fclose(infil);

    puts("MPI number is found.");
    printf("MPI number of processes in X and Y directions are %d and %d.\n",mpix,mpiy);
    puts("");
  } else{
    puts("No MPI number is found.");
    return 0;
  }

  /* Check X-array data */
  int nx[mpix];
  flag=0;
  for (m=0;m<mpix;m++){
    nx[m]=0;
    sprintf(filname,"%sx_%05d.dat",fildir,m);
    if ((infil=fopen(filname,"r"))!=NULL){
      while(fscanf(infil,"%lf",&dummy) !=EOF) nx[m]++;
      nxall+=(nx[m]-2*xoff);
      fclose(infil);
      flag++;
    }
  }
  if (flag == mpix){
    puts("X-array data is found.");
    printf("Total X grid (w.o. offset) is %d\n",nxall);
    puts("");
  } else{
    puts("No X-array data is found.");
    return 0;
  }

  /* Read, merge, and write X-array data (removing offset) */
  double x[nxall];
  isum=0;
  for (m=0;m<mpix;m++){
    double xtmp[nx[m]];
    sprintf(filname,"%sx_%05d.dat",fildir,m);
    infil=fopen(filname,"r");
    i=0;
    while(fscanf(infil,"%lf",&dummy) != EOF) xtmp[i++]=dummy;
    fclose(infil);
    for (i=0;i<nx[m]-2*xoff;i++){
      x[i+isum]=xtmp[i+xoff];
    }
    isum+=i;
  }
  sprintf(filname,"%smerge_x.dat",outdir);
  outfil=fopen(filname,"w");
  for (i=0;i<nxall;i++) fprintf(outfil,"%f\n",x[i]);
  fclose(outfil);

  /* Check Y-array data */
  int ny[mpiy];
  flag=0;
  for (m=0;m<mpiy;m++){
    ny[m]=0;
    sprintf(filname,"%sy_%05d.dat",fildir,m);
    if ((infil=fopen(filname,"r"))!=NULL){
      while(fscanf(infil,"%lf",&dummy) !=EOF) ny[m]++;
      nyall+=(ny[m]-2*yoff);
      fclose(infil);
      flag++;
    }
  }
  if (flag == mpiy){
    puts("Y-array data is found.");
    printf("Total Y grid (w.o. offset) is %d\n",nyall);
    puts("");
  } else{
    puts("No Y-array data is found.");
    return 0;
  }

  /* Read, merge, and write Y-array data (removing offset) */
  double y[nyall];
  jsum=0;
  for (m=0;m<mpiy;m++){
    double ytmp[ny[m]];
    sprintf(filname,"%sy_%05d.dat",fildir,m);
    infil=fopen(filname,"r");
    j=0;
    while(fscanf(infil,"%lf",&dummy) != EOF) ytmp[j++]=dummy;
    fclose(infil);
    for (j=0;j<ny[m]-2*yoff;j++){
      y[j+jsum]=ytmp[j+yoff];
    }
    jsum+=j;
  }
  sprintf(filname,"%smerge_y.dat",outdir);
  outfil=fopen(filname,"w");
  for (j=0;j<nyall;j++) fprintf(outfil,"%f\n",y[j]);
  fclose(outfil);

  /* Check Time-array data */
  sprintf(filname,"%st.dat",fildir);
  if ((infil=fopen(filname,"r"))!=NULL){
    while(fscanf(infil,"%lf",&dummy) != EOF) nt++;
    fclose(infil);

    puts("Time-array data is found");
    printf("Number of time step is %d.\n",nt);
    puts("");
  } else{
    puts("No Time-array data is found.");
    return 0;
  }

  /* Copy fildir/t.dat and fildir/params.dat to outdir */
  if (strcmp(fildir,outdir) != 0){
    sprintf(command,"cp %st.dat %s",fildir,outdir);
    itmp=system(command);
    sprintf(command,"cp %sparams.dat %s",fildir,outdir);
    itmp=system(command);
  }

  /* Data read, merge, write */
  int mpi=mpix*mpiy;
  float *merge_val[smax];
  for (s=0;s<smax;s++){
    merge_val[s]=(float*)malloc(sizeof(float)*nxall*nyall);
  }

  printf("Data read, merge, write ");
  for (n=nstt;n<nt;n++){
    
    putchar('.');

    /* Initialize */
    for (s=0;s<smax;s++){
      for (j=0;j<nyall;j++){
	for (i=0;i<nxall;i++) merge_val[s][nxall*j+i]=0;
      }
    }
    /* Read and merge */
    flag=0;    
    isum=0;
    jsum=0;
    for (m=0;m<mpi;m++){
      mx=m%mpix;
      my=m/mpix;

      sprintf(filname,"%soutdat_%05d_%05d.dat",fildir,n,m);
      if ((infil=fopen(filname,"rb"))!=NULL){
	flag++;
	/* Read */
	float *val=(float*)malloc(sizeof(float)*nx[mx]*ny[my]);
	for (s=0;s<smax;s++){
	  itmp=fread(val,sizeof(float),nx[mx]*ny[my],infil);
	  for (j=0;j<ny[my]-2*yoff;j++){
	    for (i=0;i<nx[mx]-2*xoff;i++){
	      merge_val[s][nxall*(j+jsum)+(i+isum)]=val[nx[mx]*(j+yoff)+(i+xoff)];
	    }
	  }
	}
	free(val);
	fclose(infil);

	/* Start position of merge_val for next step */
	int mxnext=(m+1)%mpix;
	int mynext=(m+1)/mpix;
	if ((mxnext-mx) == 1){
	  isum+=i;
	} else if ((mxnext-mx) < 0){
	  isum=0;
	}
	if ((mynext-my) == 1){
	  jsum+=j;
	} else if ((mynext-my) < 0){
	  jsum=0;
	}

      }
    }
    /* Write */
    if (flag == mpi){
      sprintf(filname,"%smerge_outdat_%05d.dat",outdir,n);
      outfil=fopen(filname,"wb");
      for (s=0;s<smax;s++) fwrite(merge_val[s],sizeof(float),nxall*nyall,outfil);
      fclose(outfil);
    }
  }
  puts(" finished.");

  for (s=0;s<smax;s++){
    free(merge_val[s]);
  }

  return 0;
}
