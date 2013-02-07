/*--main

Basic syntax: mksupersky <input list> <output image> <options>
  Preprocessing Options:
    -noscale
    -scaleregion <Xmin:Xmax,Ymin:Ymax>
  Filter Options:
    -norejectobjects
    -avsigclip <Sigma,3>
    -srcgrowrad <radius,0>
    -replacedeadpixels
  Combine Options
    -average
    -median (default)
  Output Options
    -outputmasks
    -variancetype <SIGMA, WEIGHT>
    -image_compare <template>
    -verbose <0-3>
    -version (print version and exit)

Summary:

This program combines images on a pixel-by-pixel basis to create a sky
flat (or supersky).  The combination option either a median (default) or
can be set to a mean by setting the flag -average.  A number of other
options are used to reject objects 


Detailed Description:

Scaling (-noscale -scaleregion):
  If the option -noscale is given then images are not scaled prior to 
  combination.  Otherwise, images are scaled by the median value found in
  the area specified by -scaleregion prior to combining them.  The scale
  region is set to default to [500:1500,1500:2500] where the format is: 
  [xmin:xmax,ymin:ymax].  All values are integers and non-numeric characters
  (including ".") serve as delimeters when parsing this argument.  Currently
  the same scale region is used throughout the processing. 

Filter Options:
  Four types of rejection/masking can take place on the individual images 
  prior to combining the images.  The current processing order these options
  are applied are: object masking, avsigclip, srcgrow, replacedeadpixels.

  Object rejection:
    By default object rejection occurs unless the -norejectobjects flag is 
    set.  Currently this is accomplished by setting a threshold that masks 
    all pixels with values greater than 12*sigma (5*FWHM) above the median 
    (background) level.  

  Source-grow:
    If the -srcgrowrad options is set the area around masked pixels are 
    also flagged (in order to grow the area around masked objects).  The
    exception is pixels which were masked but interpolated which are not
    grown.  The argument is in units of pixels (i.e. grow masks by N pixels
    around previously masked pixels.

  Replace Dead Pixels:
    If the -replacedeadpixels option is chosen then the routine attempts 
    to obtain a value for pixels where the BPM flag has been set.  This 
    is currently done by interpolating from the good pixels within the 
    grow radius (-srcgrowrad) of pixels with that were flagged due to the
    BPM flag.

  Sigma-clipping:
    Prior to combination if the -avsigclip option is set then another
    round of flagging will occur.  The algorithm used first calculates
    sigma for non-flagged pixels on a pixel-by-pixel basis, and considers
    the first (currently 10000) pixels to obtain an estimate of the average
    value of sigma, avsig.  Then on a pixel-by-pixel basis the mean is 
    calculated among the images and pixels are rejected when:
       abs( data[im].image[pix]/scale[im] - mean ) > N * avsig 
    Where N is the argument passed for the -avsigclip option.

Output Options:
  -outputmasks
    For debugging, this option writes the image masks (sets masked pixels 
    to a value of 1.0) to images with name <input>_mask.fits. 
  
  -variancetype <SIGMA, WEIGHT>
    An associated uncertainty (or WEIGHT) image is calculated and output along 
    with the supersky image.  The variancetype option is used to indicate the 
    type of statistic to report in the WEIGHT image plane of the output 
    (combined) image.  The current options are SIGMA and WEIGHT which are 
    evaluated relative to the value output in the image plane.  The SIGMA 
    options writes a variance type value (i.e., sqrt(sum((value-combined)^2)/N),
    while the WEIGHT options writes the value as 1/(variance^2).

  -image_compare <template>
    The image_compare option compares the resulting supersky with a previously 
    calculated image (presumably another supersky) and reports the number 
    of pixels with significant deviation.  (THIS NEEDS TO BE VERIFIED)

  -verbose <0-3>

Known "Features":
  -The -replacebadpixels option has been put in place to attempt to deal 
   with pixels that are bad in all images.  The result has not yet been 
   heavily tested.

*/

/* mksupersky 
 *
 * Combine reduced science frames into a composite image 
 * for use in the illumination and fringe corrections 
 *
 * 06/27/2012: modified by V. Kindratenko
 *  - Introduced -replacedeadpixels to replace hot/dead (BADPIX_BPM) pixels
 *    with averages over some neighborhood
 *
 */

#include "imdetrend.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

#include "argutils.h"

#define S2N_THRESH 75   /* minimum signal to noise threshold for pixel */

static const char *svn_id = "$Id$";

void print_usage(char *program_name)
{
  printf("%s <input_list> <output_image> <options>\n",program_name);
  printf("  Preprocessing Options\n");
  printf("    -noscale\n");
  printf("    -scaleregion <Xmin:Xmax,Ymin,Ymax>\n");
  printf("  Filter Options\n");
  printf("    -norejectobjects\n");
  printf("    -srcgrowrad <radius,0>\n");
  printf("    -replacedeadpixels\n");
  printf("    -avsigclip <NSigma,3>\n");
  printf("  Combine Options\n");
  printf("    -average\n");
  printf("    -median (default)\n");
  printf("  Output Options\n");
  printf("    -outputmasks\n");
  printf("    -variancetype <SIGMA, WEIGHT>\n");
  printf("  Testing Options\n");
  printf("    -image_compare <template>\n");
  printf("    -verbose <0-3>\n");
  printf("    -help    (print usage and exit)\n");
  printf("    -version (print version and exit)\n");
  
}

static int flag_scale     = YES;
static int flag_rejectob  = YES;
static int flag_combine   = MEDIAN;
static int flag_dead      = NO;
static int flag_omask     = NO;

int MakeSuperSky(int argc, char *argv[])
{
  int	i,j,im,imnum,x,y,loc,xmin,xmax,ymin,ymax,num,delta,
    xdim=0,ydim=0,ccdnum=0,ncompare,mkpath(),
    nummasked,flag_variance=NO;
  static int status=0;
  char	comment[1000],imagename[1000],command[1000],
    longcomment[10000],scaleregion[100],filter[200]="",
    obstype[200],maskname[1000],temp[1000],event[10000],inname_temp[1000],
    outname_temp[1000],*striparchiveroot(),command_line[1000],
    imtypename[6][10]={"","IMAGE","VARIANCE","MASK","SIGMA","WEIGHT"};
  float	avsigclip_sigma,srcgrowrad=0.0f,srcgrowrad2=0.0f,srcgrowrad3=0.0f,avsig,sigma,
    *vecsort=NULL,*scalefactor,*scalesort=NULL,maxval,*meanimage,
    val,sigmalim,mean,mode,fwhm,threshold,*imagepointer,
    scalefactor_output,rms,offset,maxdev;
  /*static*/ desimage *data,output,template;
  int	flag_verbose=1,flag_filter=NO,
    flag_bpm=YES,hdutype,flag_srcgrowrad=NO,
    scaleregionn[4]={500,1500,1500,2500},scalenum;
  int avsignum,srcgrowrad_int,valnum,dy,dx,newx,newy,
    flag_image_compare=NO;
  void	rd_desimage(),shell(),decodesection(),retrievescale(),
    headercheck(),printerror(),reportevt(),image_compare();
  double	avsigsum,val1sum,val2sum,weightsum,variancesum,
    variance_threshold=Squ(S2N_THRESH);
  short	mask_loc, interppix=4,maskpix=8,growpix=16,rejectpix=32, growpix2=64;
  time_t	tm;
  fitsfile *fptr;
  FILE	*inp,*pip;

  enum {OPT_NOSCALE=1,OPT_SCALEREGION,OPT_NOREJECTOBJECTS,OPT_SRCGROWRAD,
	OPT_REPLACEDEADPIXELS,OPT_AVSIGCLIP,OPT_AVERGAGE,OPT_MEDIAN,OPT_OUTPUTMASKS,
	OPT_VARIANCETYPE,OPT_IMAGE_COMPARE,OPT_VERBOSE,OPT_HELP,OPT_VERSION};
  
  if (argc<2) {
    print_usage(argv[0]);
    exit(0);
  }

  if(build_command_line(argc,argv,command_line,1000) <= 0){
    reportevt(2,STATUS,1,"Failed to record full command line.");
  }
  /* RAG: Added to print version of code to standard output (for logs) */
  sprintf(event,"%s",svn_id);
  reportevt(2,STATUS,1,event);
  reportevt(2,STATUS,1,command_line);

  /* ******************************************* */
  /* ********  PROCESS COMMAND LINE ************ */
  /* ******************************************* */

  /* first cycle through looking for quiet flag */
  for (i=1;i<argc-1;i++) {
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
    static struct option supersky_options[] =
      {
	{"scaleregion",       required_argument, 0,                 OPT_SCALEREGION},
	{"srcgrowrad",        required_argument, 0,                 OPT_SRCGROWRAD},
	{"avsigclip",         required_argument, 0,                 OPT_AVSIGCLIP},
	{"variancetype",     required_argument, 0,                  OPT_VARIANCETYPE},
        {"image_compare",     required_argument, 0,                 OPT_IMAGE_COMPARE},
	{"verbose",           required_argument, 0,                 OPT_VERBOSE},
	{"version",           no_argument,       0,                 OPT_VERSION},
	{"help",              no_argument,       0,                 OPT_HELP},
	{"average",           no_argument,       &flag_combine, AVERAGE},
	{"median",            no_argument,       &flag_combine,  MEDIAN},
        {"noscale",           no_argument,       &flag_scale,        NO},
	{"norejectobjects",   no_argument,       &flag_rejectob,     NO},
	{"replacedeadpixels", no_argument,       &flag_dead,        YES},
	{"outputmasks",       no_argument,       &flag_omask,       YES},
	{0,0,0,0}
      };

    int clopx = 0;

    clop = getopt_long_only(argc,argv,"",supersky_options,&clopx);

    if(clop == -1)
      break;
    switch(clop){
    case 0:
      // For straightforward flags
      if(supersky_options[clopx].flag != 0)
	break;
      printf("Option %s is set",supersky_options[clopx].name);
      if(optarg)
	printf(" with %s",optarg);
      printf(".\n");
      break;
    case OPT_SCALEREGION: // -scaleregion
      sprintf(scaleregion,"%s",optarg);
      decodesection(scaleregion,scaleregionn,flag_verbose); 
      break;
    case OPT_SRCGROWRAD: // -srcgrowrad
      flag_srcgrowrad=YES;
      sscanf(optarg,"%f",&srcgrowrad);
      break;
    case OPT_AVSIGCLIP: // -avsigclip
      flag_filter=AVSIGCLIP;
      sscanf(optarg,"%f",&avsigclip_sigma);
      break;
    case OPT_VARIANCETYPE: // -variancetype
      if (!strcmp(optarg,"SIGMA") || !strcmp(optarg,"sigma")) flag_variance=DES_SIGMA;
      else if (!strcmp(optarg,"WEIGHT") || !strcmp(optarg,"weight")) flag_variance=DES_WEIGHT;
      else {
	sprintf(event,"Variancetype %s undefined",optarg);
	reportevt(flag_verbose,STATUS,5,event);
	exit(1);
      }
      break;
    case OPT_IMAGE_COMPARE: // -image_compare
      flag_image_compare=YES;
      sprintf(template.name,"%s",optarg);
      if (!strncmp(&(template.name[strlen(template.name)-5]),".fits",5)  
	  && !strncmp(&(template.name[strlen(template.name)-8]),".fits.gz",8))  {
	sprintf(event,"Template image must be FITS: %s",template.name);
	reportevt(flag_verbose,STATUS,5,event);
	exit(1);
      }
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
      print_usage(argv[0]);
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
    
  /* ****************************************************************** */
  /* ******************** Test input image list *********************** */
  /* ****************************************************************** */
    
  if (!strncmp(&(inname_temp[strlen(inname_temp)-5]),".fits",5)  
      || !strncmp(&(inname_temp[strlen(inname_temp)-8]),".fits.gz",8))  {
    sprintf(event,"mksupersky requires an input image list: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  else { /* expect file containing list of fits images */
    imnum=0;
    inp=fopen(inname_temp,"r");
    if (inp==NULL) {
      sprintf(event,"File %s not found",inname_temp);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    /* *********************************************************** */
    /* * cycle through image list checking validity              * */
    /* *********************************************************** */
    while (fscanf(inp,"%s",imagename)!=EOF) {
      imnum++;
      if (strncmp(&(imagename[strlen(imagename)-5]),".fits",5)
	  || (!strncmp(&(imagename[strlen(imagename)-8]),".fits.gz",8))){
	sprintf(event,"File must contain list of FITS or compressed FITS images: %s",
		inname_temp);
	reportevt(flag_verbose,STATUS,5,event);
	exit(1);
      }
      else { /* open file and check header */
	if (fits_open_file(&fptr,imagename,READONLY,&status))  {
	  sprintf(event,"Input image didn't open: %s",imagename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
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
      exit(1);
    }
  }  /* end check of the file list validity */
    
  /* report processing details if requested */
  if (flag_verbose) {
    sprintf(event,"Input list %s contains %d FITS images",
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
    exit(1);
  }

  /* ***************************************************************** */
  /* ***************************************************************** */
  /* ********  READ and determine SCALEFACTORS for Input Images ****** */
  /* ********           also, Reject Objects                    ****** */
  /* ***************************************************************** */
  /* ***************************************************************** */

  /* create an array of image structures to hold the input images */
  data=(desimage *)calloc(imnum+1,sizeof(desimage));
  if (data==NULL) {
    sprintf(event,"Calloc for data failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }

	
  /* vector for median combining the images */
  if (flag_combine==MEDIAN) {
    vecsort=(float *)calloc(imnum,sizeof(float));
    if (vecsort==NULL) {
      sprintf(event,"Calloc for vecsort failed");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
  }
  /* vector for determining the scale of an image */
  if (flag_scale) { /* prepare for medianing the scale */
    scalenum=(scaleregionn[1]-scaleregionn[0]+1)*(scaleregionn[3]-scaleregionn[2]+1);
    scalesort=(float *)calloc(scalenum,sizeof(float));
    if (scalesort==NULL) {
      sprintf(event,"Calloc for scalesort failed");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
  }
  /* vector for holding scale factor */
  scalefactor=(float *)calloc(imnum,sizeof(float));	
  if (scalefactor==NULL) {
    sprintf(event,"Calloc for scalefactor failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }

	
  /* now cycle through input images to read and calculate scale factors */
  inp=fopen(inname_temp,"r");
  if (inp==NULL) {
    sprintf(event,"File open failed: %s",inname_temp);
    reportevt(flag_verbose, STATUS,5,event);
    exit(0);
  }
  for (im=0;im<imnum;im++) {

    fscanf(inp,"%s",data[im].name);

    /* read input image */
    rd_desimage(data+im,READONLY,flag_verbose);
    xdim= data[im].axes[0];
    ydim = data[im].axes[1];
    /* go back to the first HDU */
    if (fits_movabs_hdu(data[im].fptr,1,&hdutype,&status)) {
      sprintf("Move to hdu=1 failed: %s",data[im].name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }

    /* ************************************************************* */
    /* **** check image for various header processing parameters *** */
    /* ************************************************************* */
    headercheck(data+im,filter,&ccdnum,"DESOSCN",flag_verbose);
    headercheck(data+im,filter,&ccdnum,"DESBIAS",flag_verbose);
    /* note that we don't pupil correct Mosaic2 data, so we cannot */
    /* require that DESPUPC be present in the image */
    /*headercheck(data+im,filter,&ccdnum,"DESPUPC",flag_verbose);*/
    headercheck(data+im,filter,&ccdnum,"DESFLAT",flag_verbose);

    /* confirm OBSTYPE */
    if (fits_read_key_str(data[im].fptr,"OBSTYPE",obstype,comment,
			  &status)==KEY_NO_EXIST) {
      sprintf(event,"OBSTYPE keyword not found: %s",data[im].name);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    if (strcmp(obstype,"red")) {
      sprintf(event,"OBSTYPE not red but %s:  %s\n",obstype,
	      data[im].name);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    /* confirm image size consistent with last */
    if (data[im].axes[0]!=data[0].axes[0] || 
	data[im].axes[1]!=data[0].axes[1] || 
	data[im].npixels!=data[0].npixels) {
      sprintf(event,
	      "Current image has different size than previous image: %s",
	      data[im].name);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }

    /* close image now */
    if (fits_close_file(data[im].fptr,&status)) {
      sprintf(event,"Image close failed: %s",data[im].name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }

    /* calculate the scalefactor, mode and fwhm in the scaleregion */
    /* also read the mask image */
    if (flag_scale || flag_rejectob) 
      retrievescale(data+im,scaleregionn,scalesort,flag_verbose,
		    scalefactor+im,&mode,&fwhm);

    /* check scalefactor */
    if (flag_scale) {
      /* mask entire image if scalefactor unknown */
      if (fabs(scalefactor[im])<1.0e-4) {
	sprintf(event,"Masking entire image: %s",data[im].name);
	reportevt(flag_verbose,STATUS,3,event);
	for (i=0;i<data[im].npixels;i++) data[im].mask[i]|=maskpix;
	scalefactor[im]=1.0;
      }
    }
    else scalefactor[im]=1.0;

    /* cycle through masking objects */
    if (flag_rejectob) {
      nummasked=0;
      /* set threshold at ~12 sigma above the mode in the sky */
      threshold=mode+5.0*fwhm;
      for (i=0;i<data[im].npixels;i++) 
	if (data[im].image[i]>threshold) {
	  data[im].mask[i]|=rejectpix;
	  nummasked++;
	}
      sprintf(event,"Image= %s & MaskedPixelsAbove %.1f = %.2f percent",
	      data[im].name,threshold,(float)(100.0*nummasked/(data[im].npixels)));
      reportevt(flag_verbose,QA,1,event);
    }
  } 

  /* close input file list */	
  if (fclose(inp)) {
    sprintf(event,"Closing image list failed: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }


  /* *********************************************** */
  /* *********************************************** */
  /* *********  PREPARE OUTPUT IMAGE *************** */
  /* *********************************************** */
  /* *********************************************** */

  im=0;
  output=data[im];
  sprintf(output.name,"!%s",outname_temp);
  for (i=0;i<output.nfound;i++) output.axes[i]=data[im].axes[i];
  output.bitpix=FLOAT_IMG;
  output.npixels=data[im].npixels;
  output.image=(float *)calloc(output.npixels,sizeof(float));
  if (output.image==NULL) {
    sprintf(event,"Calloc for output.image failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  /* prepares variance image */
  if (flag_variance) {
    output.varim=(float *)calloc(output.npixels,sizeof(float));
    if (output.varim==NULL) {
      sprintf(event,"Calloc for output.varim failed");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    for (i=0;i<output.npixels;i++) output.varim[i]=0.0;
  }
  /* prepare mask image */
  output.mask=(short *)calloc(output.npixels,sizeof(short));
  if (output.mask==NULL) {
    sprintf(event,"Calloc for output.mask failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  /* set mask and weight map to zero */
  for (i=0;i<output.npixels;i++) output.mask[i]=0;

  /* ********************************************** */
  /* *********  GROW RADIUS SECTION *************** */
  /* ********************************************** */

  if (flag_srcgrowrad==YES) {
    if (flag_verbose==3) 
      reportevt(flag_verbose,STATUS,1,"Beginning source grow radius filtering"); 

    srcgrowrad_int=srcgrowrad;
    srcgrowrad2=Squ(srcgrowrad);

    for (im=0;im<imnum;im++) {
      for (y=0;y<ydim;y++) 
	for (x=0;x<xdim;x++) {
	  loc=x+y*xdim;
	  /* if bad pixel but not interpolated */
	  if (data[im].mask[loc] && !(data[im].mask[loc]&interppix) && 
	      !(data[im].mask[loc]&growpix))
	    for (dy=-srcgrowrad_int;dy<=srcgrowrad_int;dy++) {
	      newy=y+dy;
	      if (newy>=0 && newy<ydim) {
		xmax=(int)sqrt(srcgrowrad2-Squ(dy));
		for (dx=-xmax;dx<=xmax;dx++) {
		  newx=x+dx;
		  if (newx>=0 && newx<xdim) {
		    data[im].mask[newx+newy*xdim]|=growpix;
		    if (flag_dead && (data[im].mask[loc]&BADPIX_BPM))
		      data[im].mask[newx+newy*xdim]|=growpix2;
		  }
		} /* cycle through x */
	      } /* y within range? */
	    } /* cycle through y */
	} /* cycle through pixels of image */
      /* report on how many pixels are masked */
      nummasked=0;
      for (i=0;i<data[im].npixels;i++) if (data[im].mask[i]) nummasked++;
      sprintf(event,"Image= %s & AfterGrowRadiusPercentagePixelsMasked= %.2f",
	      data[im].name,(float)(100.0*nummasked/(data[im].npixels)));
      reportevt(flag_verbose,QA,1,event);
    } /* cycle through list of images */
  }
	

  /* ***************************************************************** */
  /* *********  FILTER IMAGE using AVSIGCLIP if NEEDED *************** */
  /* ***************************************************************** */

  if (flag_filter==AVSIGCLIP) {
    if (flag_verbose==3) 
      reportevt(flag_verbose,STATUS,1,"Beginning AVSIGCLIP filtering"); 
    /* first determine the average sigma */
    avsigsum=0.0;
    avsignum=0;
    for (y=scaleregionn[2]-1;y<scaleregionn[3];y++) 
      for (x=scaleregionn[0]-1;x<scaleregionn[1];x++) {
	loc=x+y*xdim;
	num=0;val1sum=val2sum=0.0;
	for (im=0;im<imnum;im++)
	  if (!data[im].mask[loc]) {
	    val=(data[im].image[loc]/scalefactor[im]);
	    val1sum+=val;
	    val2sum+=Squ(val);
	    num++;
	  }
	if (num>3) {
	  mean=val1sum/(float)num;
	  sigma=sqrt(val2sum/(float)num-Squ(mean));
	  avsigsum+=sigma;
	  avsignum++;
	  if (avsignum>AVSIGCLIPMAXNUM) break;
	}
      }
    avsig=avsigsum/(float)avsignum;
    if (flag_verbose) {
      sprintf(event,"Image= %s & AverageSigma = %12.3e",data[im].name,avsig); 
      reportevt(flag_verbose,QA,1,event);
    }

    /* now go through clipping with avsig */
    /* scale avsig by the number of sigma requested */
    avsig*=avsigclip_sigma;
    for (y=0;y<ydim;y++) 
      for (x=0;x<xdim;x++) {
	loc=x+y*xdim;
	/* first determine mean value for this pixel */
	val1sum=0.0;valnum=0;
	for (im=0;im<imnum;im++) 
	  /* if pixel not flagged or pixel has been interpolated */
	  if (!data[im].mask[loc] || data[im].mask[loc]&interppix) {
	    val1sum+=(data[im].image[loc]/scalefactor[im]);
	    valnum++;
	  }
	mean=val1sum/(float)valnum;
	/* second, flag deviant pixels */
	for (im=0;im<imnum;im++) 
	  /*if (!data[im].mask[loc] || data[im].mask[loc]&interppix) */
	  if (!data[im].mask[loc]) 
	    if (fabs((data[im].image[loc]/scalefactor[im])-mean)>avsig) 
	      data[im].mask[loc]|=maskpix;
      }
  }
	

  /* ******************************************************** */
  /* *********  RECALCULATE the SCALEFACTORs  *************** */
  /* ******************************************************** */

  if (flag_scale) {
    for (im=0;im<imnum;im++) {
      retrievescale(data+im,scaleregionn,scalesort,flag_verbose,
		    scalefactor+im,&mode,&fwhm);
      /* mask entire image of scale factor cannot be determined */
      if (fabs(scalefactor[im])<1.0e-4) {
	for (i=0;i<data[im].npixels;i++) data[im].mask[i]|=maskpix;
	scalefactor[im]=1.0;
      }
    }
  }

  if (flag_omask){
    for (im=0;im<imnum;im++){
      /* set all masked pixels to 1.0 */
      for (i=0;i<data[im].npixels;i++){
	if (data[im].mask[i]){
	  data[im].image[i]=scalefactor[im];
	}
      }
      /* now write the image */
      sprintf(temp,"%s",data[im].name);
      for (i=strlen(temp);i>=0;i--){
	if (!strncmp(temp+i,".fits",5)){
	  temp[i]=0;
	  break;
	}
      }
      sprintf(maskname,"!%s_mask.fits",temp);
	      
      if (flag_verbose) {
	sprintf(event,"Writing mask file %s",maskname+1);
	reportevt(flag_verbose,STATUS,1,event);
      }
      /* make sure path exists for new image */
      if (mkpath(maskname,flag_verbose)) {
	sprintf(event,"Failed to create path to file: %s",maskname+1);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }else{
	sprintf(event,"Created path to file: %s",maskname+1);
	reportevt(flag_verbose,STATUS,1,event);
      }
      if (fits_create_file(&fptr,maskname,&status)) {
	sprintf(event,"File creation failed: %s",maskname);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (fits_create_img(fptr,FLOAT_IMG,2,data[im].axes,&status)) {
	sprintf(event,"Image creation failed: %s",maskname);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (fits_write_img(fptr,TFLOAT,1,data[im].npixels,data[im].image,
			 &status)) {
	sprintf(event,"Image write failed: %s",maskname);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (fits_close_file(fptr,&status)) {
	sprintf(event,"File close failed: %s",maskname);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
  }
	
  /* ******************************************* */
  /* ******************************************* */
  /* *********  COMBINE SECTION *************** */
  /* ******************************************* */
  /* ******************************************* */
	
  /* Combine the input images to create composite */
  if (flag_combine==AVERAGE) {
    if (flag_verbose==3) 
      reportevt(flag_verbose,STATUS,1,"Combining images using AVERAGE");
    for (i=0;i<output.npixels;i++) {
      val1sum=weightsum=variancesum=0.0;
      valnum=0;
      for (im=0;im<imnum;im++) 
	if (!data[im].mask[i]) {
	  /* carry out a variance weighted sum */
	  val1sum+=(data[im].image[i]);
	  weightsum+=scalefactor[im];
	  valnum++;
	  variancesum+=(1.0/data[im].varim[i]/Squ(scalefactor[im]));
	}
      if (valnum) {
	output.image[i]=val1sum/weightsum;
	variancesum=Squ((float)valnum)/variancesum;
	/* mask pixels that have low S2N */
	if (variancesum<variance_threshold) output.mask[i]=1;
	if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) 
	  output.varim[i]=variancesum;
	else if (flag_variance==DES_SIGMA) {
	  if (variancesum>1.0e-10) variancesum=sqrt(variancesum);
	  output.varim[i]=1.0/variancesum;
	}
      }
      else {
	output.image[i]=1.0f;
	if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) 
	  output.varim[i]=0.0f;
	else if (flag_variance==DES_SIGMA) output.varim[i]=1.0e+6;
	output.mask[i]=1;
      }
    }
  }
	
  if (flag_combine==MEDIAN) {
    if (flag_verbose==3) 
      reportevt(flag_verbose, STATUS,1,"Combining images using MEDIAN");
    for (i=0;i<output.npixels;i++) { /* for each pixel */
      output.image[i]=1.0f;
      /* copy values into sorting vector */
      valnum=0;
      variancesum=0.0;
      for (im=0;im<imnum;im++) {
	if (!data[im].mask[i]) { 
	  vecsort[valnum]=data[im].image[i]/scalefactor[im];
	  variancesum+=(1.0/data[im].varim[i]/Squ(scalefactor[im]));
	  valnum++;
	}
      }
/*      printf("RAG: %d %d \n",i,valnum); */
/*      fflush(stdout);                   */
      if (valnum>0) {
	shell(valnum,vecsort-1);
	/* odd number of images */
	/* record median value */  
	if (valnum%2) output.image[i]=vecsort[valnum/2]; 
	else output.image[i]=0.5*(vecsort[valnum/2]+vecsort[valnum/2-1]);
	variancesum=Squ((float)valnum)/variancesum;
	/* mask pixels that have low S2N */
	if (variancesum<variance_threshold) output.mask[i]=1;
	if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) 
	  output.varim[i]=variancesum;
	else if (flag_variance==DES_SIGMA) {
	  variancesum=1.0/variancesum;
	  if (variancesum>1.0e-10) variancesum=sqrt(variancesum);
	  output.varim[i]=1.0/variancesum;
	}
      }
      else {
	output.image[i]=1.0f;
	if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) 
	  output.varim[i]=0.0f;
	else if (flag_variance==DES_SIGMA) output.varim[i]=1.0e+6;
	output.mask[i]=1;
      }
    }
  }

  /* average values for some masked pixels */ 
  if (flag_dead == YES)
    {
      if (flag_verbose==3) 
	reportevt(flag_verbose, STATUS,1,"Averaging for hot/dead pixels/columns");

      srcgrowrad_int=srcgrowrad*1.025f;
      srcgrowrad3=Squ(srcgrowrad*1.025f);

      for (y=0;y<ydim;y++)
	for (x=0;x<xdim;x++) {
	  loc=x+y*xdim;
	  if (!output.mask[loc]) continue;
	  mask_loc = 0;
	  for (im=0;im<imnum;im++) mask_loc |= data[im].mask[loc];
	  if (mask_loc&growpix2) {
	    valnum = 0;
	    val1sum = 0.0f;
	    for (dy=-srcgrowrad_int;dy<=srcgrowrad_int;dy++) {
	      newy=y+dy;
	      if (newy>=0 && newy<ydim) {
		xmax=(int)sqrt(srcgrowrad3-Squ(dy));
		for (dx=-xmax;dx<=xmax;dx++) {
		  newx=x+dx;
		  if (newx>=0 && newx<xdim) {
		    if (output.mask[newx+newy*xdim]==0) /* good pixel */
		      {
			valnum++;
			val1sum += output.image[newx+newy*xdim];
		      }
		  }
		} /* cycle through x */
	      } /* y within range? */
	    } /* cycle through y */
	    if (valnum > 0)
	      {
		mean = val1sum / valnum;
		output.image[loc] = mean;
	      }
	  }
	} /* cycle through pixels of image */
    }

  /* *********************************** */
  /* *********************************** */
  /* **** RENORMALIZE IMAGE TO 1.0 ***** */
  /* *********************************** */
  /* *********************************** */
  retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		&scalefactor_output,&mode,&fwhm);
			
  /* ************************************************************ */
  /* ************************************************************ */
  /* ************** Test image against template ***************** */
  /* ************************************************************ */
  /* ************************************************************ */
  if (flag_image_compare) {
    /*  Read template image */
    rd_desimage(&template,READONLY,flag_verbose);
    /* check image for ccdnumber */
    headercheck(&template,"NOCHECK",&ccdnum,"DESMKSKY",flag_verbose);
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

      /* ***************************************************** */
      /* ***************************************************** */
      /* ***************************************************** */
    }
  }

  /* ***************************** */
  /* ***************************** */
  /* **** WRITE OUTPUT IMAGE ***** */
  /* ***************************** */
  /* ***************************** */

  if (flag_verbose) {
    sprintf(event,"Writing results to %s",output.name+1);
    reportevt(flag_verbose,STATUS,1,event);
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
    sprintf(event,"File creation failed: %s",output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  /* create image extension */
  if (fits_create_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
    sprintf(event,"Image creation failed: %s",output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  /* write the corrected image*/
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.image,
		     &status)) {
    sprintf(event,"Image creation failed: %s",output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  /* write basic information into the header */
  if (fits_update_key_str(output.fptr,"OBSTYPE","supersky",
			  "Observation type",&status)) {
    sprintf(event,"Keyword OBSTYPE=supersky write failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);  
  }
  if (fits_update_key_str(output.fptr,"FILTER",filter,
			  "Filter name(s)",&status)) {
    sprintf(event,"Keyword FILTER=%s write failed: %s",filter,
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (fits_update_key_lng(output.fptr,"CCDNUM",ccdnum,
			  "CCD number",&status)) {
    sprintf(event,"Keyword CCDNUM=%d write failed: %s",ccdnum,
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  /* Write information into the header describing the processing */
  tm=time(NULL);
  sprintf(comment,"%s",asctime(localtime(&tm)));
  comment[strlen(comment)-1]=0;
  if (fits_write_key_str(output.fptr,"DESMKSKY",comment,
			 "sciencecombine output",&status)) {
    sprintf(event,"Keyword DESMKSKY write failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  sprintf(longcomment,"DESDM:");
  //  for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);
  sprintf(longcomment,"%s %s",longcomment,command_line);
  if (flag_verbose) {
    sprintf(event,"%s",longcomment);
    reportevt(flag_verbose,STATUS,1,event);
  }
  if (fits_write_comment(output.fptr,longcomment,&status)) {
    sprintf(event,"Command line write failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
/* */
/* Change to output image names one per line, to better accomodate large sets of files */
/* */
  sprintf(longcomment,"DESDM: %ld input files ",imnum);
  sprintf(event,"%s",longcomment);
  reportevt(flag_verbose,STATUS,1,event);
  if (fits_write_history(output.fptr,longcomment,&status)){
    sprintf(event,"Writing mksupersky call failed: %s",&(output.name[1]));
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  for (im=0;im<imnum;im++){
    sprintf(longcomment," %s",strip_path(data[im].name));
    if (flag_verbose > 3){
      sprintf(event,"%s",longcomment);
      reportevt(flag_verbose,STATUS,1,event);
    }
    if (fits_write_history(output.fptr,longcomment,&status)){
      sprintf(event,"Writing mksupersky input image list failed: %s",&(output.name[1]));
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }

  if (fits_write_key_str(output.fptr,"DES_EXT",imtypename[DES_IMAGE],
			 "Image extension",&status)) {
    sprintf(event,"Keyword DES_EXT=%s write failed: %s",
	    imtypename[DES_IMAGE],output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  /* now store the variance image that has been created or updated */
  /* first create a new extension */
  if (fits_create_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
    sprintf(event,"Creating weight map image failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  /* write the data */
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.varim,
		     &status)) {
    sprintf(event,"Writing weight map image failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (fits_update_key_str(output.fptr,"DES_EXT",
			  imtypename[flag_variance],"Extension type",&status)) {
    sprintf(event,"Keyword DES_EXT=%s write failed: %s",
	    imtypename[flag_variance],output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
  }

  /* now store the bad pixel mask that has been created or updated */
  /* first create a new extension */
  if (fits_create_img(output.fptr,USHORT_IMG,2,output.axes,&status)) {
    sprintf(event,"Creating image mask failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  /* write the data */
  if (fits_write_img(output.fptr,TUSHORT,1,output.npixels,output.mask,
		     &status)) {
    sprintf(event,"Writing image mask failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (fits_update_key_str(output.fptr,"DES_EXT",imtypename[DES_MASK],
			  "Extension type",&status)) {
    sprintf(event,"Writing keyword DES_EXT=%s failed: %s",
	    imtypename[DES_MASK],output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  /* closing the output image */
  if (fits_close_file(output.fptr,&status)) {
    sprintf(event,"Closing image failed: %s",
	    output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  return(0);

}

int main(int argc,char *argv[])
{
  return(MakeSuperSky(argc,argv));
}

#undef S2N_THRESH

