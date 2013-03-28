/*--main

Basic syntax: mkbiascor <input_list> <output_biascor_image> <options>
  Input data:
    <input_list> = ASCII file containing a list of bias FITS files 
                   to be combined
    <output_biascor_image> = output combined biascor FITS image
  Options:
    -average
    -avsigclip <Sigma>
    -median (default)
    -avminmaxclip <NCLIP>
    -variancetype <SIGMA, WEIGHT>
    -image_compare <template>
    -verbose <0-3>
    -version (print version and exit)

Summary:
  This routine combines a list of bias image frames from a single CCD to form 
  a master bias calibration frame. Images to be combined are provided through 
  a list (ASCII text file) with full paths to each frame.  The combination 
  option used is provided by a command line switch.  

Detailed Description:
  The are currently four different options that control how individual images 
  are combined:
    -median (default)
    -average
    -avminmaxclip <Nclip>
    -avsigclip <NSigma>

  The median and average combine options perform a simple pixel-by-pixel median 
  or average for all of the input images.  The avminmaxclip and avsigclip 
  variants also perform an average on a pixel-by-pixel basis but prior to this 
  calculation attempt to reject outlier values.  Specifically, the avminmaxclip
  option rejects the Nclip high and Nclip low values for each pixel prior to 
  calculating the average.  The avgsigclip option determines the local variance 
  and then rejects pixels with values outside the range of +/-(NSigma*variance) 
  before proceeding to calculate the average of the remaining values.

  -variancetype <SIGMA,WEIGHT>
    An associated uncertainty (or WEIGHT) image is calculated and output along 
    with the biascor frame.  The variancetype option is used to indicate the 
    type of statistic to report in the WEIGHT image plane of the output 
    (combined) image.  The current options are SIGMA and WEIGHT which are 
    evaluated relative to the value output in the image plane.  The SIGMA 
    options writes a variance type value (i.e., sqrt(sum((value-combined)^2)/N),
    while the WEIGHT options writes the value as 1/(variance^2).

 -image_compare <template>
    The image_compare option compares the resulting combined bias with a 
    previously calculated image (presumably a bias) and reports the number 
    of pixels with significant deviation.  (THIS NEEDS TO BE VERIFIED)

 -verbose (0-3)
    Sets the verbosity level for reporting actions being taken as the 
    processing progresses.

Known "Features":
  - If a single input image is given this program will provide a warning (once 
    for each pixel) but will output the image as a bias image.
  - There is a wait/sleep command present in most routines that write images.  
    The result is that if you write to a directory other than the current 
    directory the subroutine attempts to verify the presence of the target 
    directory and if not present builds the path.  If errors are detected the 
    routine sleeps prior to a second attempt (and does so for each level of 
    subdirectory that must be checked and/or created).

*/

/* Programmer Notes follow */
/*                         */
/* Combine bias frames into a composite bias correction that is mean to correct */
/* for residual structure perpendicular to the overscan column            */
/*  author: unknown                                                       */
/*  modified by: N. Kuropatkin   01/27/2012                               */
/*  Implemented min/max clipped combine, sigma clipped combine            */
/*  Rebuild the variance calculation to make it more efficient.           */
/* Code analysis: Given code has potential memory problem when            */
/* combining many images. The possible solution -- read acceptable        */
/* amount of images and keep results of averaging and variance            */
/* accumulating them. Rebuild method of the variance calculation          */
/* permit to achieve this.                                                */
/* N. Kuropatkin 02/06/2012.  In case of AVERAGE mode, if one image has   */
/* significant shift in signal the average value also will be shifted.    */
/* as a result calculation of variance can fail sigma limit test. In this */
/* case the assigned variance will be calculated from sigma limit and not */
/* be zero                                                                */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <getopt.h>

#include "imsupport.h"
#include "argutils.h"

#define SIGMACUT 8.0

static const char *svn_id = "$Id$";

void print_usage()
{
  printf("mkbiascor <input list> <output biascor image> <options>\n");
  printf("  -average\n");
  printf("  -avsigclip <Sigma>\n");
  printf("  -median (default)\n");
  printf("  -avminmaxclip <NCLIP>\n");
  printf("  -variancetype <SIGMA, WEIGHT>\n");
  printf("  -image_compare <template>\n");
  printf("  -verbose <0-3>\n");
  printf("  -help (print usage and exit)\n");
  printf("  -version (print version and exit)\n");
}

int 
MakeBiasCorrection(int argc,char *argv[])
{
  int	i,im,num,x,y,loca,locb,delta,ncompare=0,
    len,xmin,xmax,ymin,ymax,loc,hdutype,flag_verbose=1,
    mkpath(),flag_combine=MEDIAN,flag_variance=NO,
    flag_image_compare=0,overscantype=0,
    flag_aver = 0,flag_avsigclip=0,flag_median=0,flag_avminmaxclip=0;
  float	sigma1,sigma2,sigma3,sigma4,
    meanval1,meanval2,meanval3,meanval4;
  /* */
  float sigmapv, meanpv, minpv, maxpv,imfactor;
  int sel_pix, tailcut, k;
  /* */
  int	ampregiona[4]={0,1024,0,4096},ampregionb[4]={1024,2048,0,4096};
  static int status=0;
  int	ccdnum=0;
  int   rag_check=0;
  void	printerror(),readimsections(),overscan(),retrievescale(),desimstat();
  char	comment[1000],imagename[1000],command[1000],outname_temp[1000],inname_temp[1000],
    longcomment[10000],obstype[200],filter[100],command_line[1000],
    event[10000],*strip_path(),
    imtypename[6][10]={"","IMAGE","VARIANCE","MASK","SIGMA","WEIGHT"};
  unsigned long imnum;
  float	avsigclip_sigma,*vecsort=NULL,val,sigmalima=0,sigmalimb=0,
    offset,rms,maxdev,*imagepointer, *tempimage=NULL;
  int	minmaxclip_npix;
  double	sum2a,sum2b,sum1a,sum1b;
  void	rd_desimage(),shell(),decodesection(),headercheck(),
    reportevt(),image_compare();
  FILE	*inp=NULL;
  time_t	tm;
  desimage *data=NULL,template,datain,output;
  fitsfile *fptr=NULL;
  
  enum {OPT_AVERAGE=1,OPT_AVSIGCLIP,OPT_MEDIAN,OPT_AVMINMAXCLIP,OPT_VARIANCETYPE,
	OPT_IMAGE_COMPARE,OPT_VERBOSE,OPT_HELP,OPT_VERSION};

  if (argc<2) {
    print_usage();
    exit(0);
  }
  
  if(build_command_line(argc,argv,command_line,1000) <= 0){
    reportevt(2,STATUS,1,"Failed to record full command line.");
  }

  /* RAG: Added to print version of code to standard output (for logs) */
  
  sprintf(event,"%s",svn_id);
  reportevt(2,STATUS,1,event);
  reportevt(2,STATUS,1,command_line);

  /* ****************************************************************** */
  /* ******************* Process Command Line ************************* */
  /* ****************************************************************** */
  
  /* process command line options */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-verbose")) {
      sscanf(argv[++i],"%d",&flag_verbose);
      if (flag_verbose<0 || flag_verbose>6){
        if (flag_verbose < 0){
          sprintf(event,"Verbose level (%d) out of range (0 <= verbose <= 6). Reset to 0.",flag_verbose);
          flag_verbose=0;
        }else{
          sprintf(event,"Verbose level (%d) out of range (0 <= verbose <= 6). Reset to 6.",flag_verbose);
          flag_verbose=6;
        }
        reportevt(2,STATUS,3,event);
      }
    }
  }

  int clop;
  int cloperr = 0;
  int command_line_errors = 0;
  while(1){
    int curind = optind;
    static struct option biascor_options[] =
      {
	{"avsigclip",     required_argument, 0, OPT_AVSIGCLIP},
	{"avminmaxclip",  required_argument, 0, OPT_AVMINMAXCLIP},
	{"variancetype",  required_argument, 0, OPT_VARIANCETYPE},
	{"image_compare", required_argument, 0, OPT_IMAGE_COMPARE},
	{"verbose",       required_argument, 0, OPT_VERBOSE},
	{"average",       no_argument,       0, OPT_AVERAGE},
	{"median",        no_argument,       0, OPT_MEDIAN},
	{"version",       no_argument,       0, OPT_VERSION},
	{"help",          no_argument,       0, OPT_HELP},
	{0,0,0,0}
      };

    int clopx = 0;
    clop = getopt_long_only(argc,argv,"",
			    biascor_options,&clopx);
    if(clop == -1)
      break;
    switch(clop){
    case 0:
      // For straightforward flags
      if(biascor_options[clopx].flag != 0)
	break;
      printf("Option %s is set",biascor_options[clopx].name);
      if(optarg)
	printf(" with %s",optarg);
      printf(".\n");
      break;
    case OPT_AVSIGCLIP: // -avsigclip
      cloperr = 0;
      flag_avsigclip = 1;
      flag_combine=AVSIGCLIP;
      if(optarg)
	sscanf(optarg,"%f",&avsigclip_sigma);
      else{
	cloperr = 1;
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -avsigclip requires an argument.");
	exit(1);
      }
      break;
    case OPT_AVMINMAXCLIP: // -avminmaxclip
      cloperr = 0;
      flag_combine=CLIPPEDAVERAGE;
      flag_avminmaxclip = 1;	      
      if(optarg)
	sscanf(optarg,"%d",&minmaxclip_npix);
      else{
	cloperr = 1;
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -avminmaxclip requires an argument.");
	exit(1);
      }
      break;
    case OPT_VARIANCETYPE: // -variancetype
      cloperr = 0;
      if(optarg){
	if (!strcmp(optarg,"SIGMA")) flag_variance=DES_SIGMA;
	else if (!strcmp(optarg,"WEIGHT")) flag_variance=DES_WEIGHT;
	else {
	  sprintf(event,"Variancetype %s undefined",optarg);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
      }
      else cloperr = 1;
      if(cloperr){
	reportevt(flag_verbose,STATUS,5,
		  "Option -variancetype requires an argument.");
	command_line_errors++;
	exit(1);
      }
      break;
    case OPT_IMAGE_COMPARE: // -image_compare
      cloperr = 0;
      flag_image_compare=YES;
      if(optarg){
	sprintf(template.name,"%s",optarg);
	if (!strncmp(&(template.name[strlen(template.name)-5]),".fits",5)  
	    && !strncmp(&(template.name[strlen(template.name)-8]),".fits.gz",8))  {
	  sprintf(event,"Template image must be FITS: %s",template.name);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
      }
      else cloperr = 1;
      if(cloperr){
	reportevt(flag_verbose,STATUS,5,
		  "Option -image_compare requires a FITS image template argument.");
	command_line_errors++;
	exit(1);
      }
      break;
    case OPT_AVERAGE: // -average
      cloperr = 0;
      flag_combine=AVERAGE;
      flag_aver = 1;
      break;
    case OPT_MEDIAN: // -median
      cloperr = 0;
      flag_combine=MEDIAN;
      flag_median = 1;
      break;
    case OPT_VERBOSE: // -verbose
      // already parsed verbosity
      break;
    case OPT_VERSION: // -version
      // Version has already been printed, just exit!
      //	printf("Version: %s\n",svn_id);
      exit(0);     
      break;
    case OPT_HELP: // -help
      print_usage();
      exit(0);     
      break;
    case '?': // unknown/unrecognized
      sprintf(event,"Unrecognized option detected, ");
      if(optopt){
	sprintf(event,"%c.",optopt);
      }
      else {
	sprintf(event,"%s.",argv[curind]);
      }
      reportevt(flag_verbose,STATUS,5,event);
      command_line_errors++;
      break;
    default:
      // should never get here
      abort();
    }
  }          
    
  /* ********************************************** */
  /* ********* Handle Input Image/File  *********** */
  /* ********************************************** */
  int nnoptargs = 0;
  if(optind < argc){
    /* copy input image name */
    sprintf(inname_temp,"%s",argv[optind]);
  }
  else {
    reportevt(flag_verbose,STATUS,5,"Missing required input list of FITS images.");
    command_line_errors++;
    exit(1);
  }
  optind++;
  if(optind < argc){
    /* copy output image name */
    sprintf(outname_temp,"%s",argv[optind]);
  }
  else {
    reportevt(flag_verbose,STATUS,5,"Missing required output FITS image.");
    command_line_errors++;
    exit(1);
  }
  optind++;
  if(optind != argc){
    reportevt(flag_verbose,STATUS,5,"Extra/Unknown command line arguments encountered:");
    while(optind < argc){
      sprintf(event,"%s",argv[optind++]);
      reportevt(flag_verbose,STATUS,5,event);
    }
    command_line_errors++;
    exit(1);
  }
    
  // Finally, go ahead and fail if any command line errors were encountered - currently
  // this will only cause an exit if unrecognized/unknown options were passed to the
  // program.  
  if(command_line_errors){
    sprintf(event,"%d command line errors were detected, exiting.",command_line_errors);
    reportevt(flag_verbose,STATUS,5,event);
    exit(1);
  }
    
  if((flag_avsigclip && (flag_median || flag_aver || flag_avminmaxclip)) ||
     (flag_median    && (flag_avsigclip || flag_aver || flag_avminmaxclip)) ||
     (flag_aver      && (flag_avsigclip || flag_median || flag_avminmaxclip)) ||
     (flag_avminmaxclip && (flag_avsigclip || flag_median || flag_aver)))
    {
      reportevt(flag_verbose,STATUS,5,"Cannot use more than one combine method.");
      exit(1);
    }
    
  /* ****************************************************************** */
  /* ******************** Test input image list *********************** */
  /* ****************************************************************** */
    
  if (!strncmp(&(inname_temp[strlen(inname_temp)-5]),".fits",5)  
      || !strncmp(&(inname_temp[strlen(inname_temp)-8]),".fits.gz",8))  {
    sprintf(event,"mkbiascor requires an input image list: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  else { /* expect file containing list of bias images */
    imnum=0;
    inp=fopen(inname_temp,"r");
    if (inp==NULL) {
      sprintf(event,"File %s not found",inname_temp);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    /* *********************************************************** */
    /* * cycle through image list checking existence and OBSTYPE * */
    /* *********************************************************** */
    while (fscanf(inp,"%s",imagename)!=EOF) {
      imnum++;
      if (strncmp(&(imagename[strlen(imagename)-5]),".fits",5)
	  || !strncmp(&(imagename[strlen(imagename)-8]),".fits.gz",8)) {
	sprintf(event,"File must contain list of FITS or compressed FITS images: %s",
		inname_temp);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      else { /* open file and check header */
	if (fits_open_file(&fptr,imagename,READONLY,&status))  {
	  sprintf(event,"Input image didn't open: %s",imagename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* confirm OBSTYPE */
	if (fits_read_key_str(fptr,"OBSTYPE",obstype,comment,&status)==
	    KEY_NO_EXIST) {
	  sprintf(event,"OBSTYPE keyword not found in %s",imagename);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (strncmp(obstype,"raw_bias",5) && strncmp(obstype,"RAW_BIAS",5)) {
	  sprintf(event,"Input images doesn't have required  OBSTYPE='raw_bias' in %s",
		  imagename);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (fits_close_file(fptr,&status)) {
	  sprintf(event,"Input image didn't close: %s",imagename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
    }
    if (fclose(inp)) {
      sprintf(event,"Input image list didn't close: %s",inname_temp);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
  }  /* end check of the file list validity */
    
  /* report processing details if requested */
  if (flag_verbose) {
    sprintf(event,"Input list %s contains %ld FITS images",
	    inname_temp,imnum);
    reportevt(flag_verbose,STATUS,1,event);
  }
    
  /* ****************************************************************** */
  /* ********************* Test output image ************************** */
  /* ****************************************************************** */
    
  /* prepare output file - set to overwrite existing image if need be */
  sprintf(output.name,"!%s",outname_temp);
  /* check that it is a fits file */
  if (strncmp(&(output.name[strlen(output.name)-5]),".fits",5)) {
    sprintf(event,"Output file must be FITS file: %s",outname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
    
    
  /* ****************************************************************** */
  /* ********************** Begin processing ************************** */
  /* ****************************************************************** */
    
  /* create an array of image structures to hold the input images */
  data=(desimage *)calloc(imnum,sizeof(desimage));
  if (data==NULL) {
    sprintf(event,"Calloc failed for data");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  if (flag_combine==MEDIAN || flag_combine==CLIPPEDAVERAGE ||
      flag_combine==AVSIGCLIP) {         /* allocate memory for pixel vector */
    vecsort=(float *)calloc(imnum,sizeof(float));
    if (vecsort==NULL) {
      sprintf(event,"Calloc failed for vecsort");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
  }
	
  /* reopen input list for processing */
  inp=fopen(inname_temp,"r");
  if (inp==NULL) {
    sprintf(event,"File %s not found",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }

  /* *********************************************************** */
  /* ******* now read and overscan subtract input images ******* */
  /* *********************************************************** */
  for (im=0;im<imnum;im++) {
    /* get next image name */
    fscanf(inp,"%s",datain.name);
    /* read input image */
    rd_desimage(&datain,READONLY,flag_verbose);
    /* copy the name into the current image */
    sprintf(data[im].name,"%s",datain.name);
    /* check image for ccdnumber */
    headercheck(&datain,"NOCHECK",&ccdnum,"NOCHECK",flag_verbose);
	
    /* retrieve basic image information from header */
    /* first make sure we are at the first extension */
    if (fits_movabs_hdu(datain.fptr,1,&hdutype,&status)) {
      sprintf("Move to hdu=1 failed: %s",data[im].name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }

    /* get the BIASSEC information */
    /* skip biassec for already debiased biases: BY 2012
	     
    if (fits_read_key_str((datain.fptr),"BIASSECA",(datain.biasseca),comment,&status)
    ==KEY_NO_EXIST) {
    sprintf(event,"Keyword BIASSECA not defined in %s",datain.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
    }
    decodesection((datain.biasseca),(datain.biassecan),flag_verbose);
    */
	  
    /* get the AMPSEC information */
    if (fits_read_key_str(datain.fptr,"AMPSECA",datain.ampseca,comment,&status)
	==KEY_NO_EXIST) {
      sprintf(event,"Keyword AMPSECA not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    decodesection(datain.ampseca,datain.ampsecan,flag_verbose);

    /* get the BIASSEC information */
    /*
      if (fits_read_key_str(datain.fptr,"BIASSECB",datain.biassecb,comment,&status)
      ==KEY_NO_EXIST) {
      sprintf(event,"Keyword BIASSECB not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
      }
      decodesection(datain.biassecb,datain.biassecbn,flag_verbose);
    */
    /* get the AMPSEC information */
    if (fits_read_key_str(datain.fptr,"AMPSECB",datain.ampsecb,comment,&status)
	==KEY_NO_EXIST) {
      sprintf(event,"Keyword AMPSECB not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    decodesection(datain.ampsecb,datain.ampsecbn,flag_verbose);

    /* get the TRIMSEC information */
    /*
      if (fits_read_key_str(datain.fptr,"TRIMSEC",datain.trimsec,comment,&status)
      ==KEY_NO_EXIST) {
      sprintf(event,"Keyword TRIMSEC not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
      }
      decodesection(datain.trimsec,datain.trimsecn,flag_verbose);
    */
	  
    /* get the DATASEC information */
    /*
      if (fits_read_key_str(datain.fptr,"DATASEC",datain.datasec,comment,&status)
      ==KEY_NO_EXIST) {
      sprintf(event,"Keyword DATASEC not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
      }
      decodesection(datain.datasec,datain.datasecn,flag_verbose);
    */
	  
    /*readimsections(&(datain),flag_verbose);*/
	  
    if (flag_verbose==3) {
      /*
	sprintf(event,"BIASSECA=%s AMPSECA=%s BIASSECB=%s AMPSECB=%s TRIMSEC=%s DATASEC=%s",
	datain.biasseca,datain.ampseca,datain.biassecb,datain.ampsecb,datain.trimsec,
	datain.datasec);        
      */
      sprintf(event,"AMPSECA=%s AMPSECB=%s",
	      datain.ampseca,datain.ampsecb);
      reportevt(flag_verbose,STATUS,1,event);
    }
	  
    /* ***** OVERSCAN Section ****** */
    /* check to see if DESOSCAN keyword is set */
    /* 	  if (fits_read_keyword(datain.fptr,"DESOSCN",comment, */
    /* 	    comment,&status)==KEY_NO_EXIST) { */
    /* 	    status=0; */
    /* 	    if (flag_verbose==3) { */
    /* 	      sprintf(event,"Overscan correcting: %s",datain.name); */
    /* 	      reportevt(flag_verbose,STATUS,1,event); */
    /* 	    } */
    /* 	    overscan(&(datain),data+im,flag_verbose,overscantype); */
    /* 	  } */
    //	  else {/* image has already been OVERSCAN corrected and trimmed */
    /* simply copy the data into the output array  */
    if (flag_verbose) {
      sprintf(event,"Already overscan corrected: %s",datain.name);
      reportevt(flag_verbose,STATUS,1,event);
    }
	  
    data[im].bitpix=datain.bitpix;
    data[im].npixels=datain.npixels;
    data[im].nfound=datain.nfound;
    // RAG 2012oct12 Moved the orient_section prior to assignment for data[im] structure.  Thus,
    // datain.ampsecan datain.ampsecbn are also oriented into the positive sense 
    // as well as the data[im] equivalents (and can therefore be used when working 
    // on combined images rather than on a specific input image 
    // orients the section in the positive sense  
    orient_section(datain.ampsecan);
    orient_section(datain.ampsecbn); 
    data[im].ampsecan[0]=datain.ampsecan[0];
    data[im].ampsecan[1]=datain.ampsecan[1];
    data[im].ampsecan[2]=datain.ampsecan[2];
    data[im].ampsecan[3]=datain.ampsecan[3];
    data[im].ampsecbn[0]=datain.ampsecbn[0];
    data[im].ampsecbn[1]=datain.ampsecbn[1];
    data[im].ampsecbn[2]=datain.ampsecbn[2];
    data[im].ampsecbn[3]=datain.ampsecbn[3];
    // RAG 2012oct12 moved the orient_sections prior to assign to data[im] since some algorithms 
    // below also use datain (and need the elements sorted
    // orients the section in the positive sense
/*    orient_section(data[im].ampsecan);  */
/*    orient_section(data[im].ampsecbn);  */
    data[im].gainA=datain.gainA;data[im].gainB=datain.gainB;
    data[im].rdnoiseA=datain.rdnoiseA;data[im].rdnoiseB=datain.rdnoiseB;
    sprintf(data[im].name, "%s", datain.name);  /* Y */	    
    for (i=0;i<datain.nfound;i++) data[im].axes[i]=datain.axes[i];
    data[im].image=(float *)calloc(data[im].npixels,sizeof(float)); /* allocate image memory */
    if (data[im].image==NULL) {
      sprintf(event,"Calloc failed for data[im].image");    
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    for (i=0;i<data[im].npixels;i++) data[im].image[i]=datain.image[i];
    //	  }
	  
	  
    if (im==0) { /* prepare output image first time through*/
      output.npixels=data[im].npixels;
      output.nfound=data[im].nfound;
      for (i=0;i<output.nfound;i++) output.axes[i]=data[im].axes[i];
      output.bitpix=FLOAT_IMG;
      output.image=(float *)calloc(output.npixels,sizeof(float));
      if (output.image==NULL) {
	sprintf(event,"Calloc failed for output.image");
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      if (flag_variance) {
	output.varim=(float *)calloc(output.npixels,sizeof(float));
	if (output.varim==NULL) {
	  sprintf(event,"Calloc failed for output.varim");
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	tempimage=(float *)calloc(output.npixels,sizeof(float));
	if (tempimage==NULL) {
	  sprintf(event,"Calloc failed for temporary variance image");
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	output.variancetype=flag_variance;
      }
      if (flag_image_compare) { /* need mask for comparison */
	/* also create an image mask for the output image */
	output.mask=(short *)calloc(output.npixels,sizeof(short));
	if (output.mask==NULL) {
	  sprintf(event,"Calloc failed for output.mask");
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	/* turn off masking on all pixels */
	for (i=0;i<output.npixels;i++) output.mask[i]=0;
      }
    }
    /* close input FITS file before moving on */
    if (fits_close_file(datain.fptr,&status)) {
      sprintf(event,"Closing image failed: %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    free(datain.image);
  }  /* all the input images have been read and overscan subtracted */
	
  /* close input file list */	
  if (fclose(inp)) {
    sprintf(event,"Close file list failed: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
	  
  }
	
  /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
  /*     compute  mean, standard deviation of single CCD node     */
  /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
	
  for (im=0; im<imnum;im++){
    desimstat(data+im,ampregiona,&meanval1,&sigma1); 
    sprintf(event,"Image=%s & AmpAMean=%f & AmpAStandardDeviation =%f", 
	    data[im].name,meanval1,sigma1);
    reportevt(flag_verbose,QA,1,event);		
    desimstat(data+im,ampregionb,&meanval2,&sigma2); 
    sprintf(event,"Image=%s & AmpBMean=%f & AmpBStandardDeviation =%f", 
	    data[im].name,meanval2,sigma2);
    reportevt(flag_verbose,QA,1,event);			
  }
	
  /* *********************************************************** */
  /* ********* now combine input images to make biascor  ******* */
  /* *********************************************************** */
	
  if (flag_verbose) {
    if (flag_combine==AVERAGE) 
      sprintf(event,"Combining images using AVERAGE");
    else if (flag_combine==MEDIAN) 
      sprintf(event,"Combining images using MEDIAN");
    else if (flag_combine==AVSIGCLIP) 
      sprintf(event,"Combining images using AVSIGCLIP with +/-%.1f sigma",
	      avsigclip_sigma) ;
    else if (flag_combine==CLIPPEDAVERAGE)
      sprintf(event,"Combining images using CLIPPEDAVERAGE with clipping %d pixels",
	      minmaxclip_npix) ;
    reportevt(flag_verbose,STATUS,1,event);
  }
  /*We need to calculate limits on the variance here in case we need it */
  /* */
  /* */
  if (flag_variance ) {
    /* instead of calculated sigma of average value */
    /* we will calculate average sigma per plane as */
    /* we will use it to cut variance in a plane    */
    sum2a=sum2b=0.0;num=0;
    for (im=0;im<imnum;im++){
	    
      sum2a += data[im].rdnoiseA/data[im].gainA;
      sum2b += data[im].rdnoiseB/data[im].gainB;
    }
    /* */
    sigmalima=SIGMACUT*sum2a/(float)imnum ;
    /* introduce lower limit for sigmalima */
    if (sigmalima<1.0e-2) sigmalima=1.0e-2;
    /* */
    sigmalimb=SIGMACUT*sum2b/(float)imnum ;
    /* introduce lower limit for sigmalimb */
    if (sigmalimb<1.0e-2) sigmalimb=1.0e-2;
    if (flag_verbose) {
      sprintf(event,"Image=%s & EstimatedVarianceForAmpA/B= %.2f/%.2f",
	      output.name,sigmalima/SIGMACUT,sigmalimb/SIGMACUT);
      reportevt(flag_verbose,QA,1,event);
    }
	  
  }
  imfactor = 1.0/((float)(imnum));
  /*                                              */
  /* Combine the input images to create composite */
  if (flag_combine==AVERAGE) {
    for (i=0;i<output.npixels;i++) {
      x=i%output.axes[0];y=i/output.axes[0];
      output.image[i]=0;
      for (im=0;im<imnum;im++) output.image[i]+=data[im].image[i];
      output.image[i]*=imfactor;
      /* CALCULATE VARIANCE  */
      if (flag_variance) {
	sel_pix = 0;
	sigmapv = 0.0;
	tempimage[i] = 0.0;
	for (im=0;im<imnum;im++){
	  val=fabs(data[im].image[i]-output.image[i]);
	  if((column_in_section(x+1,datain.ampsecan) && (val < sigmalima)) ||
	     (column_in_section(x+1,datain.ampsecbn) && (val < sigmalimb))){
	    sigmapv	+=Squ(val);
	    sel_pix++;
	  }
	}
	if (sel_pix >= 2) {
	  sigmapv = sigmapv/((float)sel_pix-1.)/((float)sel_pix);  /* variance */
	}
	tempimage[i]=sigmapv;
	if(sel_pix == 0) tempimage[i]= Squ(sigmalima);
      } 
    }
  }
	
  if (flag_combine==MEDIAN) {
    for (i=0;i<output.npixels;i++) { /* for each pixel */
      x=i%output.axes[0];y=i/output.axes[0];
      /* copy values into sorting vector */
      for (im=0;im<imnum;im++) vecsort[im]=data[im].image[i];
      shell(imnum,vecsort-1);
      /* odd number of images */
      if (imnum%2) output.image[i]=vecsort[imnum/2]; /* record the median value */  
      else output.image[i]=0.5*(vecsort[imnum/2]+vecsort[imnum/2-1]);
	    
      /* CALCULATE VARIANCE  */
      if (flag_variance) {
	tempimage[i] = 0.0;
	sigmapv = 0.0;
	sel_pix = 0;
	for (im=0 ;im<imnum; im++){
	  val=fabs(data[im].image[i]-output.image[i]);
	  if((column_in_section(x+1,datain.ampsecan) && (val < sigmalima)) ||
	     (column_in_section(x+1,datain.ampsecbn) && (val < sigmalimb))){
	    sigmapv	+=Squ(val);
	    sel_pix++;
	  }
	}
	if (sel_pix >= 2) {
	  sigmapv=sigmapv/((float)sel_pix-1.)/((float)sel_pix);
	}
	tempimage[i] =sigmapv;
	if(sel_pix == 0) tempimage[i]= Squ(sigmalima);
      }
    }
  }
	
  /* ********************************************* */
  /* ************ Implementing clipped average *** */
  /* ******* By N. Kuropatkin   01/20/2012 ******* */
  /* ********************************************* */
	
  if (flag_combine==CLIPPEDAVERAGE) {
    for (i=0;i<output.npixels;i++) {
      x=i%output.axes[0];y=i/output.axes[0];
      output.image[i]=0.;
      for (im=0;im<imnum;im++) vecsort[im]=data[im].image[i];
      shell(imnum,vecsort-1);            /* sort the vector */
      if (2*minmaxclip_npix >=imnum - 2) {
	sprintf(event,"CLIPPEDAVERAGE minmaxclip_npix can not be more than imnum/2 Reset to 1");
	minmaxclip_npix = 1;
	reportevt(flag_verbose,STATUS,1,event);
      }
      for ( k=minmaxclip_npix ;k<imnum-minmaxclip_npix; k++) output.image[i]+=vecsort[k];
      output.image[i]/=(float)(imnum - minmaxclip_npix - minmaxclip_npix );
      /* CALCULATE VARIANCE  */
      if (flag_variance) {
	tempimage[i] = 0.0;
	sigmapv = 0.0;
	sel_pix = 0;
	      
	for (k=minmaxclip_npix ;k<imnum-minmaxclip_npix; k++) {
	  val=fabs(vecsort[k] - output.image[i]);
	  if((column_in_section(x+1,datain.ampsecan) && (val < sigmalima)) ||
	     (column_in_section(x+1,datain.ampsecbn) && (val < sigmalimb))){
	    sigmapv	+=Squ(val);
	    sel_pix++;
	  }
	}
	if (sel_pix >= 2) {
	  sigmapv=sigmapv/((float)sel_pix-1.)/((float)sel_pix);
	}
	      
	tempimage[i] =sigmapv;
	if(sel_pix == 0) tempimage[i]= Squ(sigmalima);
      }
    }
  }
	
  /* ********************************************* */
  /* ****** Implementing sigma clipped average *** */
  /* **** By N. Kuropatkin 01/20/2012   ********** */
  /* sigma calculation can be replaced with sigma  */
  /* calculated from noise and gain                */
  /* ********************************************* */
	
  if (flag_combine==AVSIGCLIP) {
    if ( imnum < 4) {
      sprintf(event,"Sigma clipping require at least 4 images to combine");
      reportevt(flag_verbose,STATUS,5,event);
      exit(1);
    }
    for (i=0;i<output.npixels;i++) { /* for each pixel */
      x=i%output.axes[0];y=i/output.axes[0];
      output.image[i]=0;
      sigmapv = 0.0;
      meanpv = 0.0;
      /* first cut tails of the distribution to estimate sigma */
      tailcut = (int) (0.15*imnum);  /* cut tails 30% to estimate sigma */
      if (tailcut < 1) tailcut = 1;  /* this requires at list 4 images to combine */
      for (im=0;im<imnum;im++)vecsort[im]=data[im].image[i];
      shell(imnum,vecsort-1);
      /* what is left is 2 sigma for Gaussian distribution */
      sigmapv = 0.5*(vecsort[imnum-tailcut-1] - vecsort[tailcut]);
      if (imnum%2) meanpv=vecsort[imnum/2]; /* record the median value */
      else meanpv=0.5*(vecsort[imnum/2]+vecsort[imnum/2-1]);
      /* now calculate average within selected sigma cut */
      minpv = meanpv - avsigclip_sigma*sigmapv;
      maxpv = meanpv + avsigclip_sigma*sigmapv;
      if (minpv < vecsort[0]) minpv = vecsort[0];
      if (maxpv > vecsort[imnum-1]) maxpv = vecsort[imnum - 1];
      sel_pix = 0;
      for (im=0;im<imnum;im++){
	if((data[im].image[i] >= minpv) && (data[im].image[i] <= maxpv)) {
	  output.image[i] += data[im].image[i];
	  sel_pix++;
	}
      }
      output.image[i] /= (float)sel_pix;    /* AVERAGE with sigma clip */
      /* CALCULATE VARIANCE  */
      if (flag_variance) {
	sel_pix = 0;
	sigmapv = 0.0;
	tempimage[i] = 0.0;
	for (im=0;im<imnum;im++){
	  if((data[im].image[i] >= minpv) && (data[im].image[i] <= maxpv)) {
	    val=fabs(vecsort[im] - output.image[i]);
	    if((column_in_section(x+1,datain.ampsecan) && (val < sigmalima)) ||
	       (column_in_section(x+1,datain.ampsecbn) && (val < sigmalimb))){
	      //		  if ((x<1024 && val<sigmalima) || (x>=1024 && val<sigmalimb)) {
	      sigmapv	+=Squ(val);
	      sel_pix++;
	    }
	  }
	}
	if (sel_pix >= 2) {
	  sigmapv = sigmapv/((float)sel_pix-1.)/((float)sel_pix);  /* variance */
	}
	tempimage[i]=sigmapv;
	if(sel_pix == 0) tempimage[i]= Squ(sigmalima);
      }    
    } /* end loop on pixels */
  }
	
	
  /* ******************************************* */
  /* *********  VARIANCE SECTION *************** */
  /* ******************************************* */
  /* N. Kuropatkin Here we will average variance for each pixel */
  /* using 2*delta+1,2*delta+1 matrix                           */
  /*  */
  if (flag_variance ) {
    if (flag_verbose) {
      sprintf(event,"Creating %s image",imtypename[flag_variance]);
      reportevt(flag_verbose,STATUS,1,event);
    }
	  
    sum2a=sum2b=0.0;num=0;
    delta=VARIANCE_DELTAPIXEL;
    if (flag_verbose) {
      sprintf(event,"Extracting variance using %dX%d square centered on each pixel",
	      2*delta+1,2*delta+1);
      reportevt(flag_verbose,STATUS,1,event);
    }
    for (i=0;i<output.npixels;i++) {
      x=i%output.axes[0];y=i/output.axes[0];
      /* define a small square centered on this pixel */
	    
      xmin=x-delta;xmax=x+delta;

      /* RAG for test */
      /* if (y == 2024){rag_check=1;}else{rag_check=0;} */
      /* if (rag_check == 1){printf(" %d %d %d   ",x,xmin,xmax);} */

      if(column_in_section(x+1,datain.ampsecan)){
	if(xmin < datain.ampsecan[0]){
	  xmax += (datain.ampsecan[0] - xmin);
	  xmin = datain.ampsecan[0];
	}
	if(xmax >= datain.ampsecan[1]){
	  xmin -= (xmax - datain.ampsecan[1] + 1);
	  xmax = datain.ampsecan[1]-1;
	} 
      }
      else if(column_in_section(x+1,datain.ampsecbn)){
	if(xmin < datain.ampsecbn[0]){
	  xmax += (datain.ampsecbn[0] - xmin);
	  xmin = datain.ampsecbn[0];
	}
	if(xmax >= datain.ampsecbn[1]){
	  xmin -= (xmax - datain.ampsecbn[1] + 1);
	  xmax = datain.ampsecbn[1]-1;
	} 
      }
	    
      ymin=y-delta;ymax=y+delta;
      /* hanging off bottom edge? */
      if (ymin<0) {ymax+=abs(ymin);ymin=0;}
      /* hanging off top edge? */
      if (ymax>=output.axes[1])
	{ymin-=(ymax-output.axes[1]+1);ymax=output.axes[1]-1;}
      sum2a=0.0;num=0;
      for (x=xmin;x<=xmax;x++) for (y=ymin;y<=ymax;y++) {
	      
	loc=x+y*output.axes[0];
	      

	if ((tempimage[loc] > 0.0) ) {
	  sum2a+=tempimage[loc];
	  num++;
	}
      }
	    
      /* RAG for test */
      /* if (rag_check == 1){printf(" %d %d %d %.3f \n",x,xmin,xmax,sum2a);} */

      if (num>1) val=sum2a/(float)num;
      else {
	sprintf(event,"No pixels meet criterion for inclusion (%ld,%ld)",
		i%output.axes[0],i/output.axes[0]);
	reportevt(flag_verbose,STATUS,3,event);
	      
	if(column_in_section((i%output.axes[0])+1,datain.ampsecan))
	  val=Squ(sigmalima/SIGMACUT)/(float)imnum;
	else val=Squ(sigmalimb/SIGMACUT)/(float)imnum;
      }
      /* IMPOSE VARIANCE LOWER LIMIT */
      if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) {
	if (val>1.0e-6) output.varim[i]=1.0/val;
	else output.varim[i]=1.0e+6;
      }
      else if (flag_variance==DES_SIGMA) {
	if (val>1.0e-10) val=sqrt(val);
	output.varim[i]=val;
      }
    } /* end loop on pixels  */
  }  /* end of flag_variance section */  
	
	
	
  /* ************************************************************ */
  /* ************** Test image against template ***************** */
  /* ************************************************************ */
  if (flag_image_compare) {
    /*  Read template image */
    rd_desimage(&template,READONLY,flag_verbose);
    /* check image for ccdnumber */
    headercheck(&template,"NOCHECK",&ccdnum,"DESMKBCR",flag_verbose);
    /* first compare images */
    rms=maxdev=offset=0.0;
    ncompare=0;
    image_compare(&output,&template,&offset,&rms,&maxdev,&ncompare,
		  flag_verbose);
    /* issue STATUS events according to differences measured */
	  
    /* ***************************************************** */
    /* ***************************************************** */
    /* ***************************************************** */
    /* second compare weight maps */
    if (output.varim!=NULL && template.varim!=NULL) {
      imagepointer=output.image;
      output.image=template.image=NULL;
      image_compare(&output,&template,&offset,&rms,&maxdev,&ncompare,
		    flag_verbose);
      output.image=imagepointer;
      /* issue STATUS events according to differences measured */
	    
      /* ***************************************************** */
      /* ***************************************************** */
      /* ***************************************************** */
    }
  }
	
  /* ************************************************************ */
  /* ******************* Save biascor image ********************* */
  /* ************************************************************ */
	
  if (flag_verbose) {
    sprintf(event,"Writing results to %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,1,event);
    fflush(stdout);
  }
	
  /* make sure path exists for new image */
  if (mkpath(output.name,flag_verbose)) {
    sprintf(event,"Failed to create path to file: %s",output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  else {
    sprintf(event,"Created path to file: %s",output.name+1);
    reportevt(flag_verbose,STATUS,1,event);
  }
	
  /* create the file */
  if (fits_create_file(&output.fptr,output.name,&status)) {
    sprintf(event,"Creating file failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
	
  /* create image extension */
  if (fits_create_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
    sprintf(event,"Creating image failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
	
  /* write the corrected image*/
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.image,&status)) {
    sprintf(event,"Writing image failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
	
  /* write basic information into the header */
  if (fits_write_key_str(output.fptr,"OBSTYPE","biascor",
			 "Observation type",&status)) {
    sprintf(event,"Setting OBSTYPE=biascor failed: %s",&(output.name[1]));
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
  if (fits_write_key_str(output.fptr,"DESMKBCR",comment,
			 "biascor image from mkbiascor",&status)) {
    sprintf(event,"Writing DESMKBCR failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  sprintf(longcomment,"DESDM:");
  sprintf(longcomment,"%s %s",longcomment,command_line);
  if (flag_verbose) {
    sprintf(event,"%s",longcomment);
    reportevt(flag_verbose,STATUS,1,event);
  }
  if (fits_write_history(output.fptr,longcomment,&status)) {
    sprintf(event,"Writing mkbiascor call failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  /* */
  /* Change to output image names, one per line, to better accomodate supercal */
  /* */
  sprintf(longcomment,"DESDM: %ld input files ",imnum);
  sprintf(event,"%s",longcomment);
  reportevt(flag_verbose,STATUS,1,event);
  if (fits_write_history(output.fptr,longcomment,&status)) {
    sprintf(event,"Writing mkbiascor call failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  for (im=0;im<imnum;im++){
    sprintf(longcomment,"  %s",strip_path(data[im].name));
    if (flag_verbose > 3) {
      sprintf(event,"%s",longcomment);
      reportevt(flag_verbose,STATUS,1,event);
    }
    if (fits_write_history(output.fptr,longcomment,&status)) {
      sprintf(event,"Writing mkbiascor input image failed: %s",&(output.name[1]));
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }

  if (fits_write_key_str(output.fptr,"DES_EXT",imtypename[DES_IMAGE],
			 "Image extension",&status)) {
    sprintf(event,"Writing DES_EXT=%s failed: %s",
	    imtypename[DES_IMAGE],&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  if (flag_variance) {
    /* now store the variance image that has been created or updated */
    /* first create a new extension */
    if (fits_create_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
      sprintf(event,"Creating extension for %s image failed: %s",
	      imtypename[flag_variance],&(output.name[1]));
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* write the data */	  
    if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.varim,
		       &status)) {
      sprintf(event,"Writing %s extension failed: %s",
	      imtypename[flag_variance],&(output.name[1]));
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (fits_write_key_str(output.fptr,"DES_EXT",
			   imtypename[flag_variance],"Extension type",&status)) {
      sprintf(event,"Writing DES_EXT=%s failed: %s",
	      imtypename[flag_variance],&(output.name[1]));
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }
	
  /* mean and standard deviation of the output image*/

  desimstat(&output,ampregiona,&meanval3,&sigma3);
  sprintf(event,"Image=%s & AmpAMean=%f & AmpAStandardDeviation =%f", 
	  output.name,meanval3,sigma3);
  reportevt(flag_verbose,QA,1,event);	
  desimstat(&output,ampregionb,&meanval4,&sigma4);
  sprintf(event,"Image=%s & AmpBMean=%f & AmpBStandardDeviation =%f", 
	  output.name,meanval4,sigma4);
  reportevt(flag_verbose,QA,1,event);	

  /* close the corrected image */
  if (fits_close_file(output.fptr,&status)) {
    sprintf(event,"Closing biascor image failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (flag_verbose) {
    sprintf(event,"Closed %s: 2D ( %ld X %ld )",
	    &(output.name[1]),output.axes[0],output.axes[1]);
    reportevt(flag_verbose,STATUS,1,event);
  }
	

  /* ********************************************************** */
  /* ************ free internal arrays ************************ */
  /* ********************************************************** */
  free(output.image);
  if (tempimage != NULL) free(tempimage);
  for (im=0;im<imnum;im++) free(data[im].image);

  return(0);
}


int main(int argc,char *argv[]){
  return(MakeBiasCorrection(argc,argv));
}
#undef SIGMACUT


