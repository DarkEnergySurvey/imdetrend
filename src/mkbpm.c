/* Creates bad pixel mask from biascor and flatcor */
/* $Id$ */
#include "imageproc.h"

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
  int	i,j,hdutype,hdunum,flag_verbose=1,flag_bias=YES,edgesize=15,
    flag_output=YES,flag_edgemask=NO,flag_list=NO,imnum,imoutnum,im,
    mkpath(),flag_image_compare=NO,ncompare=0,ccdnum=0,xpos,ypos,badpix=0,flag_mask=NO;
  char	comment[200],imagename[500],longcomment[10000],
    event[10000];
  float	min=0.1,max=3.0,biasmax=1000.0,maxdev,rms,offset;
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
  
  init_desimage(&flat);
  init_desimage(&bias);
  init_desimage(&output);
  init_desimage(&template);

  if (argc<3) {
    printf("mkbpm <flatcor file | list> <biascor file | list> <bpm output> <options>\n");
    printf("  -flatmax <#,3.0>\n");
    printf("  -flatmin <#,0.1>\n");
    printf("  -biasmax <#,1000>\n");
    printf("  -image_compare <template>\n");
    printf("  -verbose <0-3>\n");
    printf("  -ignorebiasimage \n");
    printf("  -mask_edges\n");	
    printf("  -edgesize <#,15> \n");
    printf("  -respectmask\n");
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
    if (!strcmp(argv[i],"-biasmax")) sscanf(argv[i+1],"%f",&biasmax);
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
    sprintf(event,"Using flatcor (value>%.2f and value<%.2f) and biascor (>%.2f) to identify bad pixels",max,min,biasmax);
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
      if (flat.image[i]>max || flat.image[i]<min) output.mask[i] |= BADPIX_BPM;
    }
    if(flat.mask && flag_mask){
      for (i=0;i<output.npixels;i++) {
	output.mask[i] |= flat.mask[i];
      }
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
	if (bias.image[i]>biasmax) output.mask[i] |= BADPIX_BPM;
      }
      if(bias.mask && flag_mask){
	for (i=0;i<output.npixels;i++) {
	  output.mask[i] |= bias.mask[i];
	}
      }
    }
  }

  if (flag_edgemask)
    {	
      for (i=0;i<output.npixels;i++) {
	xpos = i%output.axes[0] ;
	ypos = i/output.axes[0] ;
	if (xpos<edgesize || ypos<edgesize) output.mask[i] |= BADPIX_BPM;
	if (abs(xpos-output.axes[0])<=edgesize || abs(ypos-output.axes[1])<=edgesize) output.mask[i] |= BADPIX_BPM;
      }
    }

  for (i=0;i<output.npixels;i++) {
    if(output.mask[i]&BADPIX_BPM) badpix++;			 
  }	  

  if (flag_verbose)printf("Total Number of Bad Pixels masked = %d\n",badpix);
  
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

int main(int argc, char *argv[])
{
  return(mkbpm(argc,argv));
};
