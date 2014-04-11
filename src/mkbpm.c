/* mkbpm
 *
 * Creates bad pixel mask from biascor and flatcor 
 *
 * 02/27/2013: J. Marriner <marriner@fnal.gov>
 *   Give each bad pixel test a numeric code
 *   Separate bias pixels into "hot" and "warm"
 *   Warm pixels are masked, hot pixels trigger masking the entire column
 *   Change defaults to loose tolerance on flats, 20 pixels at edge
 *
 *
 * $Id: mkbpm.c 20184 2014-03-21 22:13:36Z kadrlica $ 
 *
 */

#include "imsupport.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include "argutils.h"

//generic bad pixel code
#define BADPIX 1
//test specific bad pixel codes
#define BADPIX_FLAT_MIN 1
#define BADPIX_FLAT_MAX 2
#define BADPIX_FLAT_MASK 4
#define BADPIX_BIAS_HOT 8
#define BADPIX_BIAS_WARM 16
#define BADPIX_BIAS_MASK 32
#define BADPIX_BIAS_COL 64
#define BADPIX_EDGE 128
#define BADPIX_CORR 256

int badcol(char *fname,int ccd,desimage out);


int check_for_fits_list(char *filename,int *nim)

{

  FILE *inp = NULL;

  char imagename[1000];

  char event[2000];

  int numim = 0;

  void reportevt();

  inp=fopen(filename,"r");

  if (inp==NULL) {

    sprintf(event,"File not found: %s",filename);

    reportevt(1,STATUS,5,event);

    exit(1);

  }

  while (fscanf(inp,"%s",imagename)!=EOF) {

    if (strncmp(&(imagename[strlen(imagename)-5]),".fits",5) &&

	strncmp(&(imagename[strlen(imagename)-8]),".fits.gz",8)) {

      sprintf(event,"File must contain list of FITS or compressed FITS images");

      reportevt(1,STATUS,5,event);

      exit(1);

    }

    numim++;

  }

  if (fclose(inp)) {

    sprintf(event,"Input image list didn't close: %s",filename);

    reportevt(1,STATUS,5,event);

    exit(1);

  }

  if(nim)

    *nim = numim;

  return(numim);

}



int populate_list_from_file(char **list,char *filename)

{

  FILE *inp = NULL;

  char imagename[1000];

  char event[2000];

  int numim = 0;

  void reportevt();

  inp=fopen(filename,"r");

  if (inp==NULL) {

    sprintf(event,"File not found: %s",filename);

    reportevt(1,STATUS,5,event);

    return(1);

  }

  while (fscanf(inp,"%s",imagename)!=EOF) {

    if(!list[numim]){

      sprintf(event,"List did not have enough space (%d), to hold list from %s.",numim,filename);

      reportevt(1,STATUS,5,event);

      return(1);

    }    

    sprintf(list[numim],"%s",imagename);

    numim++;

  }

  if (fclose(inp)) {

    sprintf(event,"Input image list didn't close: %s",filename);

    reportevt(1,STATUS,5,event);

    return(1);

  }

  return(0);

}

int mkbpm(int argc,char *argv[])

{

  

  static int status=0;

  void	printerror();

  int	i,j,hdutype,hdunum,flag_verbose=1,flag_bias=YES,edgesize=20,

    flag_output=YES,flag_edgemask=NO,flag_list=NO,imnum,imoutnum,im,

    mkpath(),flag_image_compare=NO,ncompare=0,ccdnum=0,xpos,ypos,badpix=0,flag_mask=NO;

  char	comment[200],imagename[500],longcomment[10000],

    event[10000];

  float	min=0.5,max=1.5,biashot=500.0,biaswarm=20.0,maxdev,rms,offset;

  desimage flat,bias,output,template;

  char **flatnames = NULL;

  char **biasnames = NULL;

  desimage *flatlist = NULL;

  desimage *biaslist = NULL;

  int nflat = 1;

  int nbias = 1;

  int initialize = 1;

  void	rd_desimage(),reportevt(),image_compare(),destroy_desimage(),headercheck();

  FILE	*inp,*out;

  time_t	tm;

  int doGeneric=0;

  int nbadcol;
  int kount[16];
  int ix, iy, iyp, iyhot, io;
  int nhotbias;
  double badfract, corrfract;
  int *mintal;
  int *maxtal;
  int *hottal;
  int *warmtal;

  int corrpix = 0;
  int num_flat_min=1;
  int num_flat_max=3;
  int num_bias_hot=2;
  int num_bias_warm=3;

  init_desimage(&flat);

  init_desimage(&bias);

  init_desimage(&output);

  init_desimage(&template);



  if (argc<3) {

    printf("mkbpm <flatcor file | list> <biascor file | list> <bpm output> <options>\n");

    printf("  -flatmax <#,1.5>\n");
    printf(" -numfmax <#,1>");

    printf("  -flatmin <#,0.5>\n");
    printf(" -numfmin <#,1>");

    printf("  -biashot <#,500>\n");
    printf(" -numbhot <#,1>");

    printf("  -biaswarm <#,20>\n");
    printf(" -numbwarm <#,1>");

    printf("  -image_compare <template>\n");

    printf("  -verbose <0-3>\n");

    printf("  -ignorebiasimage \n");

    printf("  -mask_edges\n");	

    printf("  -edgesize <#,20> \n");

    printf("  -respectmask\n");
 
    printf(" -doGeneric (set all BPM codes to BADPIX=1)");

    return(1);;

  }

  

  /* ****************************************************************** */

  /* ****************** Extract verbose flag ************************** */

  /* ****************************************************************** */

  

  for (i=4;i<argc;i++) 

    if (!strcmp(argv[i],"-verbose")) {

      sscanf(argv[++i],"%d",&flag_verbose);

      if (flag_verbose<0 || flag_verbose>3) {

	sprintf(event,"Verbose level out of range %d . Reset to 2",

		flag_verbose);

	flag_verbose=2;

	reportevt(2,STATUS,3,event);

      }

    }

  

  /* ****************************************************************** */

  /* ************** Test input and output images ********************** */

  /* ****************************************************************** */

  

  /* copy input flatcor image name if FITS file*/

  if (!strncmp(&(argv[1][strlen(argv[1])-5]),".fits",5) 

      || !strncmp(&(argv[1][strlen(argv[1])-8]),".fits.gz",8))  {

    flatnames    = (char **)calloc(1,sizeof(char *));

    flatnames[0] = (char *)calloc(1000,sizeof(char));

    sprintf(flatnames[0],"%s",argv[1]);

    imnum=1;

  }

  else if(check_for_fits_list(argv[1],&nflat)){

    flatnames = (char **)calloc(nflat,sizeof(char *));

    if (flag_verbose) {

      sprintf(event,"Input flat list %s contains %d FITS images",

	      argv[1],nflat);

      reportevt(flag_verbose,STATUS,1,event);

    }

    for(i = 0;i < nflat;i++)

      flatnames[i] = (char *)calloc(1000,sizeof(char));

    if(populate_list_from_file(flatnames,argv[1])){

      sprintf(event,"mkbpm could not process input flat list: %s",argv[1]);

      reportevt(flag_verbose,STATUS,5,event);

      return(1);;

    }

  }

  else {

    sprintf(event,"mkbpm requires an input FITS flatcor image or list of images: %s",argv[1]);

    reportevt(flag_verbose,STATUS,5,event);

    return(1);;

  }

  

  /* copy input biascor image name if FITS file*/

  if (!strncmp(&(argv[2][strlen(argv[2])-5]),".fits",5)

      || !strncmp(&(argv[2][strlen(argv[2])-8]),".fits.gz",8))  {

    biasnames = (char **)calloc(1,sizeof(char *));

    biasnames[0] = (char *)calloc(1000,sizeof(char));

    sprintf(biasnames[0],"%s",argv[2]);

    imnum=1;

  }

  else if(check_for_fits_list(argv[2],&nbias)){

    if (flag_verbose) {

      sprintf(event,"Input bias list %s contains %d FITS images",

	      argv[2],nbias);

      reportevt(flag_verbose,STATUS,1,event);

    }

    biasnames = (char **)calloc(nflat,sizeof(char *));

    for(i = 0;i < nbias;i++)

      biasnames[i] = (char *)calloc(1000,sizeof(char));

    if(populate_list_from_file(biasnames,argv[2])){

      sprintf(event,"mkbpm could not process input bias list: %s",argv[2]);

      reportevt(flag_verbose,STATUS,5,event);

      return(1);

    }

  }

  else {

    sprintf(event,"mkbpm requires an input biascor image or list of images: %s",argv[2]);

    reportevt(flag_verbose,STATUS,5,event);

    return(1);;

  }

  

  /* prepare output image */

  

  if (!strncmp(&(argv[3][strlen(argv[3])-5]),".fits",5)){

    sprintf(output.name,"!%s",argv[3]);

    output.image=NULL;

    output.varim=NULL;

    output.mask=NULL;

    output.bitpix=USHORT_IMG;

  }

  else {

    sprintf(event,"mkbpm output must be FITS image: %s",argv[3]);

    reportevt(flag_verbose,STATUS,5,event);

    return(1);

  }

  /* ********************************************** */

  /* ********* PROCESS COMMAND LINE *************** */

  /* ********************************************** */

  for (i=4;i<argc;i++) {

    if (!strcmp(argv[i],"-flatmin")) sscanf(argv[i+1],"%f",&min);

    if (!strcmp(argv[i],"-flatmax")) sscanf(argv[i+1],"%f",&max);

    if (!strcmp(argv[i],"-biashot")) sscanf(argv[i+1],"%f",&biashot);

    if (!strcmp(argv[i],"-biaswarm")) sscanf(argv[i+1],"%f",&biaswarm);

    if (!strcmp(argv[i],"-numfmin")) sscanf(argv[i+1],"%i",&num_flat_min);
    if (!strcmp(argv[i],"-numfmax")) sscanf(argv[i+1],"%i",&num_flat_max);
    if (!strcmp(argv[i],"-numbhot")) sscanf(argv[i+1],"%i",&num_bias_hot);
    if (!strcmp(argv[i],"-numbwarm")) sscanf(argv[i+1],"%i",&num_bias_warm);

     if (!strcmp(argv[i],"-doGeneric")) doGeneric=1;
 

    if (!strcmp(argv[i],"-respectmask")) flag_mask=YES;

    if (!strcmp(argv[i],"-ignorebiasimage")) flag_bias=NO;

    if (!strcmp(argv[i],"-mask_edges")) {  flag_edgemask = YES; }

    if (!strcmp(argv[i],"-edgesize")) { sscanf(argv[i+1],"%d",&edgesize);}

    if (!strcmp(argv[i],"-image_compare"))  {

      flag_image_compare=YES;

      sprintf(template.name,"%s",argv[++i]);

      if (!strncmp(&(template.name[strlen(template.name)-5]),".fits",5)

	  && !strncmp(&(template.name[strlen(template.name)-8]),".fits.gz",8))  {

	sprintf(event,"Template image must be FITS: %s",template.name);

	reportevt(flag_verbose,STATUS,5,event);

	return(1);

      }

    }

  }

  /* Issue event that contains the plan */

  if (flag_verbose) {

    sprintf(event,"Using flatcor (value>%.2f and value<%.2f) and biascor (>%.2f) to identify bad pixels",max,min,biashot);

    reportevt(flag_verbose,STATUS,1,event);

  }

  

  

  /* **************************************************** */

  /* *************** READING INPUT IMAGES *************** */

  /* **************************************************** */

  // *TODO: check ccdnum, hdutype

  /* read the flat image */

  for(j=0;j<nflat;j++){

    destroy_desimage(&flat);

    sprintf(flat.name,"%s",flatnames[j]);

    rd_desimage(&flat,READONLY,flag_verbose);

    if (flag_verbose) {

      sprintf(event,"Processing flat: %s",flat.name);

      reportevt(flag_verbose,STATUS,1,event);

    }

    /* check flat image */

    if (fits_movabs_hdu(flat.fptr,1,&hdutype,&status)) {

      sprintf(event,"Moving to HDU=1 failed: %s",flat.name);

      reportevt(flag_verbose,STATUS,5,event);

      printerror(status);

    }

    headercheck(&flat,"NOCHECK",&ccdnum,"DESMKFCR",flag_verbose);

    if(initialize==1){

      output.npixels=flat.npixels;

      output.nfound=flat.nfound;

      for (i=0;i<output.nfound;i++) output.axes[i]=flat.axes[i];

      output.mask=(short *)calloc(output.npixels,sizeof(short));

      mintal = (int *) calloc(output.npixels,sizeof (int));
      maxtal = (int *) calloc(output.npixels,sizeof (int));
      hottal = (int *) calloc(output.npixels,sizeof (int));
      warmtal = (int *) calloc(output.npixels,sizeof (int));

      if (output.mask==NULL) {

	reportevt(flag_verbose,STATUS,5,"Calloc of output.mask failed");

	return(1);

      }

      for (i=0;i<output.npixels;i++) {

	output.mask[i] = 0;

      }

      initialize = 0;

    }

    /* create the BPM */

    for (i=0;i<output.npixels;i++) {

      if (flat.image[i]>max) ++maxtal[i];
      if (flat.image[i]<min) ++mintal[i];

    }

    if(flat.mask && flag_mask)
      {
	for (i=0;i<output.npixels;i++)  output.mask[i] |= BADPIX_FLAT_MASK;
      }

  }

  /* read the bias image */

  if(flag_bias){

    for(j=0;j<nbias;j++){

      destroy_desimage(&bias);

      sprintf(bias.name,"%s",biasnames[j]);

      rd_desimage(&bias,READONLY,flag_verbose);

      if (flag_verbose) {

	sprintf(event,"Processing bias: %s",bias.name);

	reportevt(flag_verbose,STATUS,1,event);

      }

    /* check flat image */

      if (fits_movabs_hdu(bias.fptr,1,&hdutype,&status)) {

	sprintf(event,"Moving to HDU=1 failed: %s",bias.name);

	reportevt(flag_verbose,STATUS,5,event);

	printerror(status);

      }

      headercheck(&bias,"NOCHECK",&ccdnum,"DESMKBCR",flag_verbose);

      if(initialize==1){

	output.npixels=bias.npixels;

	output.nfound=bias.nfound;

	for (i=0;i<output.nfound;i++) output.axes[i]=bias.axes[i];

	output.mask=(short *)calloc(output.npixels,sizeof(short));

	if (output.mask==NULL) {

	  reportevt(flag_verbose,STATUS,5,"Calloc of output.mask failed");

	  return(1);

	}

	for (i=0;i<output.npixels;i++) {
	  output.mask[i] = 0;
	}

	initialize = 0;

      }

      /* create the BPM */

      for (i=0;i<output.npixels;i++) {
	if (bias.image[i]>biashot)
	  {
	    nhotbias = 1;
	    io = i - output.axes[0];
	    if (io>=0)
	      {
		if (bias.image[i]<0.2*bias.image[io]) nhotbias=0;
	      }
	    io = i + output.axes[0];
	    if (io<output.npixels)
	      {
		if (bias.image[i]<0.2*bias.image[io]) nhotbias=0;
	      }
	    if (nhotbias) ++hottal[i];
	  }
	if (bias.image[i]>biaswarm) ++warmtal[i];
      }

      if(bias.mask && flag_mask){
	
	for (i=0;i<output.npixels;i++) {
	  output.mask[i] |= BADPIX_BIAS_MASK;
	}
	
      }
      
    }
    
  }


  if (flag_edgemask)

    {	

      for (i=0;i<output.npixels;i++) {

	xpos = i%output.axes[0] ;

	ypos = i/output.axes[0] ;

	if (xpos<edgesize || ypos<edgesize) output.mask[i] |= BADPIX_EDGE;

	if (abs(xpos-output.axes[0])<=edgesize || abs(ypos-output.axes[1])<=edgesize) output.mask[i] |= BADPIX_EDGE;

      }

    }

  //Flag bad pixels based on majority logic
  for (i=0;i<output.npixels;++i)
    {
      if (mintal[i]>=num_flat_min && (output.mask[i]&BADPIX_EDGE)==0) 
	output.mask[i]|=BADPIX_FLAT_MIN;
      if (maxtal[i]>=num_flat_max) output.mask[i]|=BADPIX_FLAT_MAX;
      if (warmtal[i]>=num_bias_warm) output.mask[i]|=BADPIX_BIAS_WARM;
      if (hottal[i]>=num_bias_hot) output.mask[i]|=BADPIX_BIAS_HOT;
    }
  for (ix=0;ix<output.axes[0];++ix)
    {
      nhotbias = 0;
      for (iy=0;iy<output.axes[1];++iy)
	{
	  io = iy*output.axes[0] + ix;
	  if (output.mask[io]&BADPIX_BIAS_HOT) 
	    {
	      printf("hotpixel x=%i y=%i \n",ix,iy);
	      ++nhotbias;
	      iyhot = iy;
	      for (iyp=iyhot-5;iyp<=iyhot+5;++iyp)
		{
		  io = iyp*output.axes[0] + ix;
		  if (io>=0 && io<output.npixels) output.mask[io]|=BADPIX_BIAS_WARM;
		  io = iyp*output.axes[0] + ix - 1;
		  if (io>=0 && io<output.npixels) output.mask[io]|=BADPIX_BIAS_WARM;
		  io = iyp*output.axes[0] + ix + 1;
		  if (io>=0 && io<output.npixels) output.mask[io]|=BADPIX_BIAS_WARM;
		}
	    }
	}
      //If there are hot bias pixels we want to flag surrounding pixels
      //as bad or correctable
      if (nhotbias==0) continue;
      printf("column=%i row=%i  nhotbias=%i \n",ix,iyhot,nhotbias);
      //Correctable column has only 1 hot pixel
      if (nhotbias==1)
	{
	  printf("Column=%i is correctable.\n",ix);
	  if (ccdnum<=31)
	    {
	      for (iy=iyhot+60;iy<output.axes[1];iy+=60) 
		{
		  io = iy*output.axes[0] + ix;
		  output.mask[io] &= ~BADPIX_BIAS_WARM;
		  output.mask[io] |= BADPIX_CORR;
		}
	      for (iy=0;iy<iyhot-5;++iy)
		{
		  io= iy*output.axes[0] + ix;
		  output.mask[io] &= ~BADPIX_BIAS_WARM;
		  output.mask[io] |= BADPIX_CORR;
		}
	    }
	  else
	    {
	      for (iy=iyhot-60;iy>=0;iy-=60)
		{
		  io = iy*output.axes[0] + ix;
		  output.mask[io] &= ~BADPIX_BIAS_WARM;
		  output.mask[io] |= BADPIX_CORR;
		}
	      for (iy=output.axes[1];iy>iyhot-5;--iy)
		{
		  io= iy*output.axes[0] + ix;
		  output.mask[io] &= ~BADPIX_BIAS_WARM;
		  output.mask[io] |= BADPIX_CORR;
		}
	    }
	}
      //Not correctable. Flag whole column as bad
      else
	{
	  printf("Column=%i is NOT correctable.\n",ix);
	  for (iy=0;iy<output.axes[1];++iy) 
	    {
	      io= iy*output.axes[0] + ix;
	      output.mask[io] |= BADPIX_BIAS_COL;
	    }
	}	    
    }
	  
  printf("CCDNUM=%i\n",ccdnum);
  nbadcol = badcol("funky_column.lis",ccdnum,output);
  //Count bad pixels by type
  for (j=0;j<16;++j) kount[j] = 0;
  badpix = 0;
  corrpix = 0;
  int col_bias_hot, col_bias_warm, col_flat_max, col_flat_min;
  int col_sum_bias_hot=0;
  int col_sum_bias_warm=0;
  int col_sum_flat_max=0;
  int col_sum_flat_min=0;
  int col_sum_corr=0;
  int col_sum_hot_corr=0;
  int count_edge = 0;
  int col_bad, col_corr;
  int count_bias_hot = 0;
  int count_bias_warm = 0;
  int count_flat_min = 0;
  int count_flat_max = 0;
  int correct;

  for (ix=0;ix<output.axes[0];++ix)
    {
      col_bias_hot = 0;
      col_bias_warm = 0;
      col_flat_min = 0;
      col_flat_max = 0;
      col_bad = 0;
      col_corr = 0;

      for (iy=0;iy<output.axes[1];++iy)
	{
	  io = iy*output.axes[0] + ix;
	  // Clear warm bias pixels in bad columns
	  if (output.mask[io]&BADPIX_CORR) output.mask[io] &= ~BADPIX_BIAS_WARM;
	  // Also at edges
	  if (output.mask[io]&BADPIX_EDGE) output.mask[io] &= ~BADPIX_BIAS_WARM;
	  for (j=0;j<16;++j) if (output.mask[io] & 1<<j) ++kount[j]; 
	  if (output.mask[io]&(~BADPIX_CORR)) ++badpix;
	  if (output.mask[io]&BADPIX_CORR) ++corrpix;
	  //If enabled -- statistics will be screwed up (maybe)
	  if (doGeneric && (output.mask[io]!=0)) output.mask[io] = BADPIX;
	  if (output.mask[io]&BADPIX_EDGE)
	    {
	      ++count_edge;
	      continue;
	    }
	  if (output.mask[io]&BADPIX_BIAS_HOT)
	    {
	      ++count_bias_hot;
	      col_bias_hot = 1;
	      continue;
	    }
	  if (output.mask[io]&BADPIX_FLAT_MIN)
	    {
	      ++count_flat_min;
	      col_flat_min = 1;
	      continue;
	    }
	  if (output.mask[io]&BADPIX_FLAT_MAX)
	    {
	      ++count_flat_max;
	      col_flat_max = 1;
	      continue;
	    }
	  if (output.mask[io]&BADPIX_BIAS_WARM)
	    {
	      ++count_bias_warm;
	      col_bias_warm = 1;
	      continue;
	    }
	}
      if (col_flat_max) ++col_sum_flat_max;
      if (col_flat_min) ++col_sum_flat_min;
      if (col_bias_warm) ++col_sum_bias_warm;
      if (col_bias_hot) ++col_sum_bias_hot;
      correct = (output.mask[ix]&BADPIX_CORR) 
	|| (output.mask[ix+output.axes[1]-1]&BADPIX_CORR);
      if (correct) ++col_sum_corr;
      if (col_bias_hot && correct) ++col_sum_hot_corr;
    }	  

  //This diagnostic print assumes particular values for BADPIX_*
  printf("Flat low=%i (%i < %6.3f)\n",kount[0],num_flat_min,min);
  printf("Flat high=%i (%i > %6.3f \n",kount[1],num_flat_max,max);
  printf("Flat mask=%i \n",kount[2]);
  printf("Bias hot=%i (%i > %7.2f)\n",kount[3],num_bias_hot,biashot);
  printf("Bias warm=%i (%i > %7.2f)\n",kount[4],num_bias_warm,biaswarm);
  printf("Bias mask=%i \n",kount[5]);
  printf("Bias column=%i \n",kount[6]);
  printf("Edge=%i (%i)\n",kount[7],edgesize);
  printf("Correctable=%i \n",kount[8]);
  badfract = (100.0*badpix)/output.npixels;
  corrfract = (100.0*corrpix)/output.npixels;
  if (flag_verbose)
    {
      printf("Total Number of bad pixels masked = %d (%6.2f\%)\n",
	     badpix,badfract);
      printf("Total Number of correctable pixels masked = %d (%6.2f\%)\n",
	     corrpix,corrfract);
    }
  printf("pix edge=%i hot=%i min=%i max=%i warm=%i \n",count_edge,count_bias_hot,count_flat_min,count_flat_max,count_bias_warm);
  printf("col hot=%i min=%i max=%i warm=%i corr=%i hot_corr=%i\n",col_sum_bias_hot,col_sum_flat_min,col_sum_flat_max,col_sum_bias_warm,col_sum_corr,col_sum_hot_corr);

  /*
  // ADW: Comment out creation of summary file
  char filnam[100];
  FILE *fout;
  sprintf(filnam,"Summary_%i.txt",ccdnum);
  fout=fopen(filnam,"w");
  fprintf(fout,"%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n",
	  count_edge,count_bias_hot,count_flat_min,count_flat_max,count_bias_warm,
	  col_sum_bias_hot,col_sum_flat_min,col_sum_flat_max,col_sum_bias_warm,
	  col_sum_corr,col_sum_hot_corr,badpix,corrpix,nbadcol);
  fclose(fout);
  */

  /* ************************************************************ */

  /* *************** Test bpm against template ****************** */

  /* ************************************************************ */

  if (flag_image_compare) {

    /*  Read template image */

    rd_desimage(&template,READONLY,flag_verbose);

    /* check image for ccdnumber */

    headercheck(&template,"NOCHECK",&ccdnum,"DESMKBPM",flag_verbose);

    template.image=template.varim=NULL; 

    /* make sure that image pointers are NULL */

    if (output.image!=NULL || template.image!=NULL) {

      sprintf(event,"Image and weight map pointers must be NULLed to compare masks");

      reportevt(flag_verbose,STATUS,5,event);

      return(1);

    }

    rms=maxdev=offset=0.0;

    ncompare=0;

    image_compare(&output,&template,&offset,&rms,&maxdev,&ncompare,

		  flag_verbose);

    /* issue STATUS events according to differences measured */

    

    /* ***************************************************** */

    /* ***************************************************** */

    /* ***************************************************** */

  }

  



  /* *********************** */	

  /* **** SAVE RESULTS ***** */

  /* *********************** */

  

  /* make sure path exists for new image */

  if (mkpath(output.name,flag_verbose)) {

    sprintf(event,"Failed to create path to file: %s",output.name+1);

    reportevt(flag_verbose,STATUS,5,event);

    return(1);

  }

  else {

    sprintf(event,"Created path to file: %s",output.name+1);

    reportevt(flag_verbose,STATUS,1,event);

  }

  

  /* create the file */

  if (fits_create_file(&output.fptr,output.name,&status)) {

    sprintf(event,"Creating file failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  

  /* create image HDU */

  if (fits_create_img(output.fptr,USHORT_IMG,2,output.axes,&status)) {

    sprintf(event,"Creating image failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  

  /* write the corrected image*/

  if (fits_write_img(output.fptr,TUSHORT,1,output.npixels,output.mask,

		     &status))  {

    sprintf(event,"Writing image failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  

  if (fits_write_key_lng(output.fptr,"CCDNUM",ccdnum,

			 "CCD Number",&status)) {

    sprintf(event,"Setting CCDNUM=%d failed: %s",ccdnum,output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  

  /* Write processing history into the header */

  /* get system time */

  tm=time(NULL);

  sprintf(comment,"%s",asctime(localtime(&tm)));

  comment[strlen(comment)-1]=0;

  if (fits_write_key_str(output.fptr,"DESMKBPM",comment,"bad pixel map created",&status)) {

    sprintf(event,"Writing DESMKBPM keyword failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  sprintf(longcomment,"DESDM:");

  for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);

  if (flag_verbose) reportevt(flag_verbose,STATUS,1,longcomment);

  if (fits_write_comment(output.fptr,longcomment,&status)) {

    sprintf(event,"Writing longcomment failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  if (fits_write_key_str(output.fptr,"DES_EXT","MASK","Extension type",&status)) {

    sprintf(event,"Writing DES_EXT=MASK failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  

  /* close the corrected image */

  if (fits_close_file(output.fptr,&status)) {

    sprintf(event,"File close failed: %s",output.name);

    reportevt(flag_verbose,STATUS,5,event);

    printerror(status);

  }

  if (flag_verbose) {

    sprintf(event,"Closed %s: 2D ( %ld X %ld )",

	    &(output.name[flag_output]),output.axes[0],output.axes[1]);	  

    reportevt(flag_verbose,STATUS,1,event);

  }

  return(0);

}

int badcol(char *fname,int ccdnum,desimage out)
{
  FILE* fin;
  int i;
  char line[128];
  int ccd, col;
  float  d;
  float dev[2049];
  int bcol[2049];
  int seq=0;
  int nc;
  int naxis1, naxis2;
  int start;
  double size, sum;
  int icol, ncol;

  const float MINTHRES = -5.0;
  for (i=1;i<=2048;++i) dev[i]=0.0;
  fin = fopen("funky_column.lst","r");
  if (fin > 0) 
    {
      // If successful, read the funky columns
      while (fgets(line,127,fin))
        {
          sscanf(line,"%i %i %f",&ccd,&col,&d);
          if (ccd != ccdnum) continue;
          printf("ccd=%i col=%i dev=%f\n",ccd,col,d);
          if (col>=1 && col<=2048) dev[col] = d;
        }
    }
  seq=0;
  ncol = 0;
  for (i=1;i<1024;++i)
    {
      //See if part of a previous sequence
      if (seq!=0 && i<=seq+2)
	{
	  if ((dev[i]>0.0 && (dev[i-1]!=0.0 || dev[i-2]<0.0)) || 
	      (dev[i]<0.0 && dev[i-1]<0.0)  || 
	      (dev[i]==0.0 && dev[i-1]<0.0))
	    {
	      //      printf(" seq %i %f \n",i,dev[i]);
	      seq = i;
	      sum += dev[i];
	      bcol[icol] = i;
	      ++icol;
	      continue;
	    }
	}
      if (seq!=0 && fabs(sum)/size<1.0) 
	{
	  printf("start=%i size=%f sum=%f \n",start,size,sum);
	  ncol = icol;
	}
      //See if start of a new sequence
      seq = 0;
      icol = ncol;
      if (dev[i]>MINTHRES) continue;
      //      printf(" Start %i %f \n",i,dev[i]);
      seq = i;
      start = seq;
      size = dev[i];
      sum = size;
      bcol[icol] = i;
      ++icol;
    }
  if (seq!=0 && fabs(sum)/size<1.0) 
    {
      printf("start=%i size=%f sum=%f \n",start,size,sum);
      ncol = icol;
    }

  seq = 0;
  for (i=2048;i>=1025;--i)
    {
      //See if part of a previous sequence
      if (seq!=0 && i>=seq-2)
	{
	  if ((dev[i]>0.0 && (dev[i+1]!=0.0 || dev[i+2]<0.0)) || 
	      (dev[i]<0.0 && dev[i+1]<0.0)  || 
	      (dev[i]==0.0 && dev[i+1]<0.0))
	    {
	      //    printf(" seq %i %f \n",i,dev[i]);
	      seq = i;
	      sum += dev[i];
	      bcol[icol] = i;
	      ++icol;
	      continue;
	    }
	}
      if (seq!=0 && fabs(sum)/size<1.0) 
	{
	  printf("start=%i size=%f sum=%f \n",start,size,sum);
	  ncol = icol;
	}
      //See if start of a new sequence
      seq = 0;
      if (dev[i]>MINTHRES) continue;
      //     printf(" Start %i %f \n",i,dev[i]);
      seq = i;
      start = seq;
      size = dev[i];
      sum = size;
      bcol[icol] = i;
      ++icol;
    }

  if (seq!=0 && fabs(sum)/size<1.0) 
    {
      printf("Bad column start=%i size=%f sum=%f \n",start,size,sum);
      ncol = icol;
    }
  naxis1 = out.axes[0];
  naxis2 = out.axes[1];
  printf("Found %i bad columns \n",ncol);
  for (nc=0;nc<ncol;++nc)
    {
      icol = bcol[nc];
      printf("Masking bad column=%i \n",icol);
      for (i=0;i<naxis2;++i)
	{
	  out.mask[i*naxis1+icol] |= BADPIX_CORR;
	}
    }
     
  return(ncol);
}

int main(int argc, char *argv[])

{
  return(mkbpm(argc,argv));
};

