/*--main
Basic Syntax: imcorrect <input image or list> <options>
  Input Data
    <input image or list> = FITS image (or ASCII file containing a list 
                                        of FITS images)
  Image Corrections:
    -bpm <image>
    -obpm <image>
    -bias <image>
    -pupil <image>
    -flatten <image>
    -darkcor <image>  (see warning below)
    -photflatten <image>
    -illumination <image>
    -fringe <image>
  Other input options:
    -scaleregion <Xmin:Xmax,Ymin,Ymax>
    -overscantype <0-4 default=0>
    -nooverscan
  Output Options
    -output <image or file list>
    -MEF (single file output)
    -variancetype <SIGMA, WEIGHT, INVERSE_VARIANCE>
    -noisemodel <SKYONLY (default), FULL>
    -interpolate_col <variance scaling>
    -minsize <pixels> src removal in WEIGHT maps
    -maxsize <pixels> src removal in WEIGHT maps
    -updatesky
    -verbose <0-3>
    -version (Print version and exit)
    -help    (Print usage w/ options list)

Summary:
  This program applies corrections to remove instrumental effects from images.
  A wide variety of corrections are possible including:
    - overscan subtraction (depricated but still supported)
    - application of a bad pixel mask (BPM)
    - bias subtraction
    - pupil ghost subtraction
    - flat field correction 
    - dark current subtraction
    - pixel area correction
    - illumination correction
    - fringe pattern subtraction
  Most of these rely upon other routines to produce the actual correction that
  is applied.  In all cases the correction images should have the same phsical
  dimensions as the input images (or the trimmed dimensions if overscan 
  subtraction is chosen).

Detailed Description:

The detrending steps within imcorrect occur in the order listed below.
In some cases corrections occur as a group with a null correction
used for corrections not explicitly requested on the command line.
In the descriptions below, each correction is described in isolation
and arithmatic expressions refer to the correction applied in isolation.

Overscan Correction (-overscantype): 
  The first to be applied is an overscan correction.  In the current 
  DESDM framework this step is usually accomplished by DECam_crosstalk. 
  The -overscantype option is still supported within mkflatcor but 
  only functions correctly for -overscantype 0 which removes an estimate
  for the instrumental bias by considering a line-by-line median in the 
  BIASSEC for each amplifier.  After processing the overscan (as well
  as the pre/post-scan sections) of the image are trimmed/removed.
  If the DESOSCN keyword has been previously populated (i.e. an overscan 
  correction has already been applied) then this option is ignored.
  By default overscan correction is attempted unless the -nooverscan 
  option is selected (or the DESOSCN keyword has been populated in the
  input image).

Bias Correction (-bias <biascor>):
  When used, a master bias image will be read and subtracted from each 
  input frame.  The image used is FITS file with dimensions that match
  the post overscan image size (i.e. biascor image should be trimmed).

BPM (-bpm <image> or -obpm <image> ):
  If these option are used then a bad pixel mask (BPM) is propogated into
  the mask plane of the input image.  When -bpm is used this occurs through
  a logical OR (thus any pixels already masked in the input image are 
  respected).  When -obpm is used then the mask in the input image
  is discarded and replaced by the values present in the BPM image.
  Pixels flagged in the BPM are not considered when calculating scale 
  factors.  This mask is later propagated into the output image.

Pupil Correction (-pupil <image>)
  When present an image of the pupil ghost is subtracted from each of the
  input images.  This is accomplished by determining the scale (median) of
  the input image (see the -scale option below) and then subtracting a 
  scaled version of the pupil image:
           output(pix)=input(pix)-scale*pupil(pix)

Flat Field (-flatten <image>)
  When present a flat field correction is applied to the input image(s).
  The correction that occurs has the form:
           output(pix)=input(pix)/flat(pix)   

Dark Current (-dark <image>)
  WARNING: This option has never been tested and the order in which it 
  would occur (i.e., its placement in the code) is suspect.  Currently
  the DESDM/CP frameworks do not apply a dark correction but if it were
  applied, dark subtraction would occur after the (bias/pupil/flat).  
  The form of the correction is: 
           output(pix)=input(pix)-(exptime*darkcor(pix))

Pixel Area Correction (-photflatten <image>)
  The pixel area correction attempts to compensate for the varying area 
  subtended by each pixel on the sky.  This is done by determining the
  sky background level (sky) in the scaleing region and then applying the
  correction to the difference between the flux in each pixel and sky level.
  This correction has the form:
           output(pix)= [ [input(pix) - sky ] * photflat(pix) ] + sky  

Illumination Correction (-illumcor <image>)
  If the -illumcor option is present the an illumination correction is made.
  The correction has a form similar to that used for the flat field.

Fringe Correction (-fringe <image>)
  This option causes a fringe correction to be applied.  Immediately, prior 
  to performing the correction the median value of the image is determined
  and the fringe image is scaled by this value prior to attempting to subtract
  the fringe patter.  The correction has the form:
           output(pix) = input(pix) - [ scale * fringe(pix) ] 


Other Input Options:

Scaling (-scaleregion)
  A number of corrections (pupil, fringe) require a scale factor as they
  are made relative to the background level of the image.  The region over
  which this scale factor is determined can be changed through the 
  -scaleregion option.  The default region used is [500:1500,1500:2500]
  where the format is: [xmin:xmax,ymin:ymax], all values are integers and 
  non-numeric characters (including ".") serve as delimeters when parsing
  this argument.  Currently the same scale region is used throughout the
  code.


Output Options:

Output File Name:
  By default the output file name for each input image (or image within the 
  input list) is also used for output.  Thus input images are overwritten
  unless a new image name/path is given.  The output file name (or lists of
  names) is specified with the -output option.

Output File Format
  By default up to three output files are written for each input image.  Those
  three images have names ending in: "_im.fits" (the image data), "_bpm.fits"
  (the mask data), and "_var.fits" (the weight/uncertainty data).  

  In order to write a standard DESDM image, one using a multi-extension FITS
  (MEF) format, the -MEF flag must be specified.  When selected, no suffixes
  are appended to the output image name, instead an MEF file is written where
  the primary HDU contains the image data, the second HDU the mask data and 
  third HDU the image weight/uncertainty data.  

Weight/Uncertainty Images (-variancetype and -noisemodel):
  The default is not to write a weight image, however, two major options 
  (along with some associated tuning parameters) control the types 
  of weight/uncertainty images which are produced.

  The actual type of weight/uncertainty calculated depends on the 
  -variancetype option.  This option has three possible values: SIGMA, WEIGHT,
  and INVERSE_VARIANCE.  The values in the output weight image from these 
  different types are set after reading the input image as follows:

             SIGMA:  sqrt( image[pix]/gain[pix] + sqr( readnoise/gain ) )
            WEIGHT:  1/sqr(SIGMA)
  INVERSE_VARIANCE: is currently set to the same value as for the option WEIGHT

  These values are also altered based upon the detrending steps that occur.
  Specifically:  
    - flat field correction --> SIGMA/=sqr(flat)
    - illumination correction --> SIGMA/=sqr(illum)
    - bias variance corrections --> SIGMA+=1.0/bias_var, 
    - flat field variance  --> SIGMA+=1.0/flat_var

  The current implementation recomputes the weight image (i.e. does not
  respect any incoming (previously calculated) weight image but rather
  recomputes from scratch.

  The -noisemodel flag determines whether the weight image will reflect
  the presence of sources.  The default, -noisemodel SKYONLY attempts to 
  generate a weight map with sources removed, while -noisemodel FULL will 
  include the effect of increased Poisson noise associated with sources.
  The process to remove sources is more complex than space permits here
  but uses a median filter that samples MAXPIXELS (currently hardcoded to 32) 
  from the area around each pixel.  The samples are acquired using a random
  pattern within a box around each pixel with minumum dimension (-minsize) 
  near the image edges and with maximum dimension set by (-maxsize).

Other Output Options:

  -interpolate_col <variance scaling>

  When this option is chosen, an interpolation of pixels flagged as saturated
  or bad (in the mask image) occurs for features (2 columns wide or smaller).
  The argument controls the output weights which are scaled by the factor
  provided.  This value is currently restricted to be between 1. and 10.

  -updatesky

  When this option is chose the median sky level and uncertainty are 
  calculated based on the scale section (see -scaleregion) and the header 
  keywords SKYBRITE and SKYSIGMA are populated/updated.  Care should be
  used with this option as these values are later used by Mangle to estimate
  the survey depth.  


Known Bugs/Features:
  - The implementation of the dark correction is suspect.  No test data 
    has been presented and the option has never been carefully tested.
  - Input weigth maps are not respected but always recalculated
  - The algorithm used to remove sources from the weight image (-noisemodel
    SKYONLY) currently has the following idiosychrosies:
       - works with a very small window near the chip edges.
       - when the full window is considered (i.e. more than -maxsize from
         an image edge) the selection of values is not truly random but 
         rather has a fixed pattern.
       - masked data are included among the samples used to generate the
         SKYONLY weights
  - the overscan options (soon to be depricated) only work for -overscantype 0
    (i.e. removal of the overscan with a row-by-row median)

*/

/* Carries out basic detrending of CCD data */

#include "imsupport.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include "argutils.h"

#define MAXPIXELS 32

//char *strip_path(const char *path);

static const char *svn_id = "$Id$";

void print_usage(char *program_name)
{
      printf("%s <image or file list> <options>\n",program_name);
      printf("  Input Data\n");
      printf("    <input image or list> = FITS image (or ASCII file containing a list of FITS images)\n");
      printf("  Image Corrections:\n");
      printf("    -bpm <image>\n");
      printf("    -obpm <image>\n");
      printf("    -bias <image>\n");
      printf("    -pupil <image>\n");
      printf("    -flatten <image>\n");
      printf("    -darkcor <image>\n");
      printf("    -photflatten <image>\n");
      printf("    -illumination <image>\n");
      printf("    -fringe <image>\n");
      printf("  Other Input Options: \n");
      printf("    -scaleregion <Xmin:Xmax,Ymin,Ymax>\n");
      printf("  Output Options:\n");
      printf("    -output <image or file list>\n");
      printf("    -MEF (single file output)\n");
      printf("    -variancetype <SIGMA, WEIGHT, INVERSE_VARIANCE>\n");
      printf("    -noisemodel <SKYONLY (default), FULL>\n");
      printf("    -interpolate_col <variance scaling>\n");
      printf("    -minsize <pixels,def=4> src removal in WEIGHT maps\n");
      printf("    -maxsize <pixels,def=128> src removal in WEIGHT maps\n");
      printf("    -updatesky\n");
      printf("    -verbose <0-3>\n");
      printf("    -version (Print version and exit)\n");
      printf("    -help    (Print this list and exit)\n");
}



int ImCorrect(int argc,char *argv[])
{
  
  int	hdunum,x,y,nfound,imtype,c,
    xmin,xmax,anynull,maskedpixels,interppixels,
    saturatepixels,mkpath(),i,keysexist,j,n,nkeys,
    scaleregionn[4]={500,1500,1500,2500},
    miniscaleregionn[4],scalenum,ccdnum=0,
      flag_overscan=YES,flag_bias=NO,flag_flatten=NO,flag_verbose=1,
      flag_output=NO,flag_list=NO,flag_bpm=0,flag_variance=NO,
      flag_noisemodel=DES_SKYONLY,
      imnum,imoutnum,im,hdutype,flag_interpolate_col=NO,
      flag_illumination=NO,flag_fringe=NO,flag_dark=NO,
      flag_pupil=NO,flag_impupil=NO,flag_photflatten=NO,
      flag_imoverscan=NO,flag_imbias=NO,flag_imflatten=NO,flag_imphotflatten,
      flag_imillumination,flag_imfringe,flag_updateskybrite=NO,
      flag_updateskysigma=NO,minsize=4,maxsize=128,
      seed=-15,xlow,xhi,flag_mef=NO,flag_bpm_override=0,
      ymax,ymin,dy,dx,loc,count,ylen,xlen,k,totpix,l,overscantype=0;
    static	int status=0;
    double	dbzero,dbscale,sumval;
    long	axes[2],naxes[2],pixels,npixels,bscale,bzero,bitpix,fpixel,
      ranseed=-564;
    char	comment[1000],longcomment[10000],filter[100]="",imagename[1000],
      bpmname[1000],rootname[1000],varimname[1000],outputlist[1000],input_list_name[1000],
      *striparchiveroot(),event[10000],command_line[1000],
      keycard[100],keyname[10],scaleregion[100],
      imtypename[6][10]={"","IMAGE","VARIANCE","MASK","SIGMA","WEIGHT"};
    float	scale,offset,gasdev(),scale_interpolate,maxsaturate,imval,
      overscanA=0.0,overscanB=0.0,scalefactor,mode,ran1(),
      *randnum=NULL,*vecsort,*scalesort,skybrite,skysigma,thresholdval;
    desimage bias,photflat,flat,darkimage,data,output,bpm,illumination,
      fringe,pupil,nosource;
    int	check_image_name();
#if 0
    void	rd_desimage(),headercheck(),printerror(),decodesection(),
      readimsections(),overscan(),retrievescale(),reportevt(),
      shell();
#endif
    short accept_mask;
    time_t	tm;
    FILE	*inp=NULL;
    FILE        *out=NULL;
    /* list of keywords to be removed from the header after imcorrect */
    char	delkeys[100][10]={"CCDSUM","AMPSECA","AMPSECB","TRIMSEC",
				  "BIASSECA","BIASSECB","GAINA","GAINB","RDNOISEA",
				  "RDNOISEB","SATURATA","SATURATB",
				  ""};

    enum {OPT_BPM=1,OPT_OBPM,OPT_BIAS,OPT_PUPIL,OPT_FLATTEN,OPT_DARKCOR,OPT_PHOTFLATTEN,OPT_ILLUMINATION,
	  OPT_FRINGE, OPT_SCALEREGION,OPT_OUTPUT,OPT_MEF,OPT_VARIANCETYPE,OPT_NOISEMODEL,OPT_INTERPOLATE_COL,
	  OPT_MINSIZE,OPT_MAXSIZE,OPT_UPDATESKY,OPT_VERBOSE,OPT_VERSION,OPT_HELP,OPT_RANSEED};
    
    nosource.image = NULL;
    imoutnum = 0;

    if(build_command_line(argc,argv,command_line,1000) <= 0){
      reportevt(2,STATUS,1,"Failed to record full command line.");
    }

    sprintf(event,"%s",svn_id);
    reportevt(2,STATUS,1,event); 
    reportevt(2,STATUS,1,command_line);

    if (argc<2) {
      print_usage(argv[0]);
      exit(1);
    }


    // ** PROCESS COMMAND LINE **

    // Still need to scan command line for verbosity so we can 
    // properly print status messages during command line parse.
    for (i=1;i<argc;i++) 
      {
	if (!strcmp(argv[i],"-verbose") || !strcmp(argv[i],"-v")) {
	  sscanf(argv[++i],"%d",&flag_verbose);
	  if (flag_verbose<0 || flag_verbose>3) {
	    sprintf(event,"Verbose level out of range %d. Reset to 2",
		    flag_verbose);
	    flag_verbose=2;
	    reportevt(2,STATUS,3,event);
	  }
	}
      }


    int clop;
    int cloperr = 0;
    int command_line_errors = 0;
    while(1){
      int curind = optind;
      static struct option imcorrect_options[] =
	{
	  {"bpm",             required_argument, 0,OPT_BPM},
	  {"obpm",            required_argument, 0,OPT_OBPM},
	  {"bias",            required_argument, 0,OPT_BIAS},
	  {"pupil",           required_argument, 0,OPT_PUPIL},
	  {"flatten",         required_argument, 0,OPT_FLATTEN},
	  {"darkcor",         required_argument, 0,OPT_DARKCOR},
	  {"photflatten",     required_argument, 0,OPT_PHOTFLATTEN},
	  {"fringe",          required_argument, 0,OPT_FRINGE},
	  {"scaleregion",     required_argument, 0,OPT_SCALEREGION},
	  {"illumination",    required_argument, 0,OPT_ILLUMINATION},
	  {"output",          required_argument, 0,OPT_OUTPUT},
	  {"variancetype",    required_argument, 0,OPT_VARIANCETYPE},
	  {"noisemodel",      required_argument, 0,OPT_NOISEMODEL},
	  {"interpolate_col", required_argument, 0,OPT_INTERPOLATE_COL},
	  {"minsize",         required_argument, 0,OPT_MINSIZE},
	  {"maxsize",         required_argument, 0,OPT_MAXSIZE},
	  {"ranseed",         required_argument, 0,OPT_RANSEED},
	  {"verbose",         required_argument, 0,OPT_VERBOSE},
	  {"updatesky",       no_argument,       0,OPT_UPDATESKY},
	  {"MEF",             no_argument,       0,OPT_MEF},
	  {"version",         no_argument,       0,OPT_VERSION},
	  {"help",            no_argument,       0,OPT_HELP},
	  {0,0,0,0}
	};
      int clopx = 0;
      clop = getopt_long_only(argc,argv,"",
			      imcorrect_options,&clopx);
      if(clop == -1)
	break;
      switch(clop){
      case 0:
	// For straightforward flags
	if(imcorrect_options[clopx].flag != 0)
	  break;
	printf("Option %s is set",imcorrect_options[clopx].name);
	if(optarg)
	  printf(" with %s",optarg);
	printf(".\n");
	break;
      case OPT_BPM: // -bpm
	cloperr = 0;
	flag_bpm=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(bpm.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Bad Pixel Mask image must follow -bpm");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_OBPM: // -obpm
	cloperr = 0;
	flag_bpm=YES;
	flag_bpm_override=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(bpm.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Bad Pixel Mask image must follow -obpm");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_BIAS: // -bias
	cloperr = 0;
	flag_bias=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(bias.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Bias correction image must follow -bias");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_PUPIL: // -pupil
	cloperr = 0;
	flag_pupil=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(pupil.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Pupil correction image must follow -pupil");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_FLATTEN: // -flatten
	cloperr = 0;
	flag_flatten=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(flat.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Flat correction image must follow -flatten");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_DARKCOR: // -darkcor
	cloperr = 0;
	flag_dark=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(darkimage.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Master dark correction image must follow -dark");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_PHOTFLATTEN: // -photflatten
	cloperr = 0;
	flag_photflatten=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(photflat.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Photometricflat correction image must follow -photflatten");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_ILLUMINATION: // -illumination
	cloperr = 0;
	flag_illumination=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(illumination.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Illumination correction image must follow -illumination");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_FRINGE: // -fringe
	cloperr = 0;
	flag_fringe=YES;
	if(optarg){
	  if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    cloperr = 1;
	  }	    
	  else sprintf(fringe.name,"%s",optarg);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Fringe correction image must follow -fringe");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_SCALEREGION: // -scaleregion
	cloperr = 0;
	if(optarg){
	  sprintf(scaleregion,"%s",optarg);
	  decodesection(scaleregion,scaleregionn,flag_verbose);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Region specification must follow -scaleregion");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_OUTPUT: // -output
	cloperr = 0;
	flag_output=YES;
	if(optarg){
	  if (check_image_name(optarg,CHECK_FITS,flag_verbose)) {
	    imoutnum = 1;
	    sprintf(output.name,"%s",optarg);
	  }
	  // REMEMBER TO CHECK FOR LIST COMPATABILITY LATER
	  // if(flag_list) DISABLED FOR NOW
	  //	    sprintf(output.name,"%s",argv[i]);
	  //	  else { /* error because input list cannot map to single image */
	  //	    reportevt(flag_verbose,STATUS,5,
	  //		      "File list must follow -output option in multi image mode");
	  //	    exit(1);
	  //	  }
	  //	  }
	  else {
	    // THIS CHECK TO BE DONE LATER
	    //	  if (!flag_list) { /* input image needs single output image */
	    //	  	    reportevt(flag_verbose,STATUS,5,
	    //	  		      "FITS image output name must follow -output option in single image mode");
	    //	  	    exit(1);
	    //	  }
	    sprintf(outputlist,"%s",optarg);
	    imoutnum=0;
	    out=fopen(outputlist,"r");
	    if (out==NULL) {
	      sprintf(event,"Output file list not found: %s",outputlist);
	      reportevt(flag_verbose,STATUS,5,event);
	      command_line_errors++;
	      exit(1);
	    }
	    while (fscanf(out,"%s",imagename)!=EOF) {
	      imoutnum++;
	      if (!check_image_name(imagename,CHECK_FITS,flag_verbose)) {
		reportevt(flag_verbose,STATUS,5,
			  "Output list must contain FITS or compressed FITS images");
		command_line_errors++;
		exit(1);
	      }
	    }
	    if (fclose(out)) {
	      sprintf(event,"Closing file failed: %s",outputlist);
	      reportevt(flag_verbose,STATUS,5,event);
	      command_line_errors++;
	      exit(1);
	    }
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Output FITS image or list must follow -output");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_VARIANCETYPE: // -variancetype
	cloperr = 0;
	if(optarg){
	  if (!strcmp(optarg,"SIGMA")) flag_variance=DES_SIGMA;
	  else if (!strcmp(optarg,"WEIGHT")) flag_variance=DES_WEIGHT;
	  else if (!strcmp(optarg,"INVERSE_VARIANCE")) 
	    flag_variance=DES_VARIANCE;
	  else {
	    sprintf(event,"Variancetype %s undefined",optarg);
	    reportevt(flag_verbose,STATUS,5,event);
	    command_line_errors++;
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Variance type specification must follow -variancetype");
	  command_line_errors++;
	  exit(1);	 
	}
	break;
      case OPT_NOISEMODEL: // -noisemodel
	cloperr = 0;
	if(optarg){
	  if (!strcmp(optarg,"SKYONLY")) flag_noisemodel=DES_SKYONLY;
	  else if (!strcmp(optarg,"FULL")) flag_noisemodel=DES_FULL;
	  else {
	    sprintf(event,"Noisemodel %s undefined",optarg);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Noisemodel specification must follow -noisemodel");
	  command_line_errors++;
	  exit(1);	 	  
	}
	break;
      case OPT_INTERPOLATE_COL: // -interpolate_col
	cloperr = 0;
	flag_interpolate_col=YES;
	if(optarg){
	  sscanf(optarg,"%f",&scale_interpolate);
	  if (scale_interpolate<1.0 || scale_interpolate>10.0) {
	    sprintf(event,"Noise scale factor for interpolated pixels must be between 1 and 10 (%f)",
		    scale_interpolate);
	    reportevt(flag_verbose,STATUS,5,event);
	    command_line_errors++;
	    exit(1);
	  }
	  else {
	    sprintf(event,
		    "Pixel noise scale factor for interpolated pixels %.2f",
		    scale_interpolate);
	    reportevt(flag_verbose,STATUS,1,event);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Pixel noise scale factor must follow -interpolate_col");
	  command_line_errors++;
	  exit(1);	 	  
	}
	break;
      case OPT_MINSIZE: // -minsize
	cloperr = 0;
	if(optarg){
	  sscanf(optarg,"%d",&minsize);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Min size specification must follow -minsize");
	  command_line_errors++;
	  exit(1);	 	  	  
	}
	break;
      case OPT_MAXSIZE: // -maxsize
	cloperr = 0;
	if(optarg){
	  sscanf(optarg,"%d",&maxsize);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Max size specification must follow -maxsize");
	  command_line_errors++;
	  exit(1);	 	  	  
	}
	break;
      case OPT_VERBOSE: // -verbose
	// already parsed verbosity
	break;
      case OPT_RANSEED: // -ranseed
	cloperr = 0;
	if(optarg){
	  sscanf(optarg,"%ld",&ranseed);
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Random seed specification must follow -ranseed");
	  command_line_errors++;
	  exit(1);	 	  	  
	}
	break;
      case OPT_UPDATESKY: // -updatesky
	flag_updateskybrite=YES;
	flag_updateskysigma=YES;
	break;
      case OPT_MEF: // -MEF
	flag_mef=YES;
	break;
      case OPT_VERSION: // -version
	// Actually, the version has already been printed above.  Just exit!
	//	printf("Version: %s\n",svn_id);
	exit(0);     
   	break;
      case OPT_HELP: // -help
	print_usage(argv[0]);
	exit(0);     
	break;
      case '?': // unknown/unrecognized
	// Actually, getopt has already printed an error to stderr.
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
    imnum = 0;
    if(optind < argc){
      /* copy input image name if FITS file*/
      if(check_image_name(argv[optind],CHECK_FITS,flag_verbose)){
	sprintf(data.name,"%s",argv[optind]);
	imnum = 1;
      }
      else { /* expect file containing list of images */
	imnum=0;flag_list=YES;
	sprintf(input_list_name,"%s",argv[optind]);
	inp=fopen(argv[optind],"r");
	if (inp==NULL) {
	  sprintf(event,"File not found: %s",argv[optind]);
	  reportevt(flag_verbose,STATUS,5,event);
	  command_line_errors++;
	  exit(1);
	}
	while (fscanf(inp,"%s",imagename)!=EOF) {
	  imnum++;
	  if (!check_image_name(imagename,CHECK_FITS,flag_verbose)) {
	    reportevt(flag_verbose,STATUS,5,
		      "Input image list must contain FITS or compressed FITS images");
	    command_line_errors++;
	    exit(1);
	  }
	}
	if (fclose(inp)) {
	  reportevt(flag_verbose,STATUS,5,"Closing input image list failed");
	}
      }
    }
    else {
      reportevt(flag_verbose,STATUS,5,"Missing required input FITS image or input list.");
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
    
    // NOW make sure that the input image and output image (or lists) jibe
    if (imoutnum==imnum) {
      if (flag_verbose==3) 
	reportevt(flag_verbose,STATUS,1,
		  "Output and input lists of same length");
    }
    else {
      sprintf(event,"Input(%d) and output(%d) lists must be of same length in multi image mode",
	      imoutnum,imnum);
      reportevt(flag_verbose,STATUS,5,event);
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


    vecsort=NULL;
    /* ******************************************** */
    /* ********** Prepare for Processing  ********* */
    /* ******************************************** */
    if(flag_list)
      sprintf(event,"Input list %s contains %d FITS images",input_list_name,imnum);
    else
      sprintf(event,"Input is single FITS image: %s",data.name);
    reportevt(flag_verbose,STATUS,1,event);

    /* set up scaling vectors */
    scalenum=(scaleregionn[1]-scaleregionn[0]+1)*
      (scaleregionn[3]-scaleregionn[2]+1);
    scalesort=(float *)calloc(scalenum,sizeof(float));
    if (scalesort==NULL) {
      reportevt(flag_verbose,STATUS,5,"Calloc of scalesort failed");
      exit(0);
    }

    /* reopen input/output lists if in multiimage mode */
    if (flag_list) {
      /* reopen file for processing */
      /*inp=fopen(argv[1],"r");
	if (inp==NULL) {
	sprintf(event,"File not found: %s",argv[1]);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
	}
      */
      /* reopen output list for processing */
      out=fopen(outputlist,"r");
      if (out==NULL) {
	sprintf(event,"File not found: %s",outputlist);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
    }


    /* ****************************************** */
    /* ********** READ CALIBRATION DATA ********* */
    /* ****************************************** */
			
    /* read bias image */
    if (flag_bias) {
      rd_desimage(&bias,READONLY,flag_verbose);
      /* confirm that this is actual bias image */
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(bias.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",bias.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is bias image */
      headercheck(&bias,"NOCHECK",&ccdnum,"DESMKBCR",flag_verbose);
      if (fits_close_file(bias.fptr,&status)) {
	sprintf(event,"File close failed: %s",bias.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
    /* read dark image */
    if (flag_dark) {
      rd_desimage(&darkimage,READONLY,flag_verbose);
      /* confirm that this is actual bias image */
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(darkimage.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",bias.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is bias image */
      headercheck(&darkimage,"NOCHECK",&ccdnum,"DESMKDCR",flag_verbose);
      if (fits_close_file(bias.fptr,&status)) {
	sprintf(event,"File close failed: %s",bias.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }


    /* read flat image */	
    if (flag_flatten) {
      rd_desimage(&flat,READONLY,flag_verbose);
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(flat.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",flat.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is flat image */
      headercheck(&flat,filter,&ccdnum,"DESMKFCR",flag_verbose);
      if (fits_close_file(flat.fptr,&status)) {
	sprintf(event,"File close failed: %s",flat.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
    /* read photflat image */	
    if (flag_photflatten) {
      rd_desimage(&photflat,READONLY,flag_verbose);
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(photflat.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",photflat.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is photflat image */
      headercheck(&photflat,"NOCHECK",&ccdnum,"DESMKPFC",flag_verbose);
      if (fits_close_file(photflat.fptr,&status)) {
	sprintf(event,"File close failed: %s",photflat.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
    /* read pupil image */
    if (flag_pupil) {
      rd_desimage(&pupil,READONLY,flag_verbose);
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(pupil.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",pupil.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is pupil ghost image */
      headercheck(&pupil,"NOCHECK",&ccdnum,"DESMKPUP",flag_verbose);
      if (fits_close_file(pupil.fptr,&status)) {
	sprintf(event,"File close failed: %s",pupil.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
    /* read bad pixel mask image */
    if (flag_bpm || flag_bpm_override) {
      rd_desimage(&bpm,READONLY,flag_verbose);
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(bpm.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",bpm.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is bpm image */
      headercheck(&bpm,"NOCHECK",&ccdnum,"DESMKBPM",flag_verbose);
      if (fits_close_file(bpm.fptr,&status)) {
	sprintf(event,"File close failed: %s",bpm.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
    /* create a mock bpm regardless */
    else bpm=bias;	
		
    /* read illumination image */	
    if (flag_illumination) {
      rd_desimage(&illumination,READONLY,flag_verbose);
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(illumination.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",illumination.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is illumination image */
      headercheck(&illumination,filter,&ccdnum,"DESMKICR",flag_verbose);
      if (fits_close_file(illumination.fptr,&status)) {
	sprintf(event,"File close failed: %s",illumination.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
		
    /* read fringe image */	
    if (flag_fringe) {
      rd_desimage(&fringe,READONLY,flag_verbose);
      /* make sure we are in the 1st extension */
      if (fits_movabs_hdu(fringe.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",fringe.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm that this is fringe image */
      headercheck(&fringe,filter,&ccdnum,"DESMKFRG",flag_verbose);
      if (fits_close_file(fringe.fptr,&status)) {
	sprintf(event,"File close failed: %s",fringe.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
    }
		
		
    /********************************************/
    /****** PROCESSING SECTION BEGINS HERE ******/
    /********************************************/
    /* now cycle through input images to process them */
    for (im=0;im<imnum;im++) {
      /* get next image name */
      if (flag_list) {
	inp=fopen(input_list_name,"r");
	if (inp==NULL) {
	  sprintf(event,"Failed to open input file list: %s",input_list_name);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	for (c=0;c<=im;c++) fscanf(inp,"%s",data.name);
	if (fclose(inp)) {
	  sprintf(event,"Failed to close input file list: %s",input_list_name);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
      }
      
      /* ************************************** */
      /* ***** Read Input Image Section  ****** */
      /* ************************************** */

      if (flag_verbose==3) {
	sprintf(event,"Reading INPUT image %s",data.name);
	reportevt(flag_verbose,STATUS,1,event);
      }
      if (flag_output) { /* writing output to new image */
	rd_desimage(&data,READONLY,flag_verbose);  
	/* prepare output image */
	if (out) /* involked with a file list */
	  fscanf(out,"%s",imagename);
	else /* involked with single file */
	  sprintf(imagename,"%s",output.name);
	sprintf(output.name,"!%s",imagename);
      }
      else  /* overwriting input image */
	rd_desimage(&data,READWRITE,flag_verbose);
      /*echo exposure time if applying dark current correction */
      if (flag_dark) {
	sprintf(event, "Using exposure time of %f seconds for dark current correction\n",data.exptime);
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* prepare output image */
      output.npixels=data.npixels;
      sprintf(event,"Input image has %ld pixels.",data.npixels);
      reportevt(flag_verbose,STATUS,1,event);
      output.nfound=data.nfound;
      output.variancetype=data.variancetype;
      for (i=0;i<output.nfound;i++) 
	output.axes[i]=data.axes[i];
      output.bitpix=FLOAT_IMG;
      output.saturateA=data.saturateA;
      output.saturateB=data.saturateB;
      output.gainA=data.gainA;
      output.gainB=data.gainB;
      output.rdnoiseA=data.rdnoiseA;
      output.rdnoiseB=data.rdnoiseB;
      if (output.saturateA>output.saturateB) 
	maxsaturate=output.saturateA;
      else 
	maxsaturate=output.saturateB;

		  
      /* ************************************** */
      /* ***** READING IMAGE SECTIONS    ****** */
      /* ************************************** */
	  
      /* retrieve basic image information from header */
      if (fits_movabs_hdu(data.fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to HDU=1 failed: %s",data.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* confirm correct filter and ccdnum */
      headercheck(&data,filter,&ccdnum,"NOCHECK",flag_verbose);

      /* Why is this sitting here and then appearing 50 lines down again? */
      /* get the DATASEC information */
      if (fits_read_key_str(data.fptr,"DATASEC",data.datasec,
			    comment,&status) ==KEY_NO_EXIST) {
	sprintf(event,"Keyword DATASEC not found: %s",data.name);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      decodesection(data.datasec,data.datasecn,flag_verbose);
  
      /* get the AMPSEC information */
      if (fits_read_key_str(data.fptr,"AMPSECA",data.ampseca,
			    comment,&status)==KEY_NO_EXIST) {
	sprintf(event,"Keyword AMPSECA not defined: %s",data.name);
	reportevt(flag_verbose,STATUS,3,event);
	//	printerror(status);
	status = 0;
	data.ampsecan[0] = 1;
	data.ampsecan[1] = 1024;
	data.ampsecan[2] = 1;
	data.ampsecan[3] = 4096;
      }
      else
	decodesection(data.ampseca,data.ampsecan,flag_verbose);
      
      /* get the AMPSEC information */
      if (fits_read_key_str(data.fptr,"AMPSECB",data.ampsecb,
			    comment,&status) ==KEY_NO_EXIST) {
	sprintf(event,"Keyword AMPSECB not defined: %s",data.name);
	reportevt(flag_verbose,STATUS,3,event);
	//	printerror(status);
	status = 0;
	data.ampsecbn[0] = 1;
	data.ampsecbn[1] = 1024;
	data.ampsecbn[2] = 1;
	data.ampsecbn[3] = 4096;
      }
      else
	decodesection(data.ampsecb,data.ampsecbn,flag_verbose);

      output.ampsecan[0] = data.ampsecan[0];
      output.ampsecan[1] = data.ampsecan[1];
      output.ampsecan[2] = data.ampsecan[2];
      output.ampsecan[3] = data.ampsecan[3];
      output.ampsecbn[0] = data.ampsecbn[0];
      output.ampsecbn[1] = data.ampsecbn[1];
      output.ampsecbn[2] = data.ampsecbn[2];
      output.ampsecbn[3] = data.ampsecbn[3];

      /* ************************************** */
      /* ***** OVERSCAN and TRIM Section ****** */
      /* ************************************** */

/*       if (flag_overscan) { */
/* 	if (fits_read_keyword(data.fptr,"DESOSCN",comment,comment,&status)== */
/* 	    KEY_NO_EXIST) { */
/* 	  flag_imoverscan=1; */
/* 	  status=0; /\* reset flag *\/ */
/* 	  /\* get the required header information *\/ */
/* 	  /\* get the BIASSEC information *\/           */
/* 	  if (fits_read_key_str((data.fptr),"BIASSECA",(data.biasseca), */
/* 				comment,&status)==KEY_NO_EXIST) { */
/* 	    sprintf(event,"Keyword BIASSECA not defined: %s",data.name); */
/* 	    reportevt(flag_verbose,STATUS,5,event); */
/* 	    printerror(status); */
/* 	  } */
/* 	  decodesection((data.biasseca),(data.biassecan),flag_verbose); */
    
/* 	  /\* get the BIASSEC information *\/ */
/* 	  if (fits_read_key_str(data.fptr,"BIASSECB",data.biassecb, */
/* 				comment,&status) ==KEY_NO_EXIST) { */
/* 	    sprintf(event,"Keyword BIASSECB not defined: %s",data.name); */
/* 	    reportevt(flag_verbose,STATUS,5,event); */
/* 	    printerror(status); */
/* 	  } */
/* 	  decodesection(data.biassecb,data.biassecbn,flag_verbose); */

/* 	  /\* get the TRIMSEC information *\/ */
/* 	  if (fits_read_key_str(data.fptr,"TRIMSEC",data.trimsec, */
/* 				comment,&status) ==KEY_NO_EXIST) { */
/* 	    sprintf(event,"Keyword TRIMSEC not defined: %s",data.name); */
/* 	    reportevt(flag_verbose,STATUS,5,event); */
/* 	    printerror(status); */
/* 	  } */
/* 	  decodesection(data.trimsec,data.trimsecn,flag_verbose); */
/* 	  /\* report header parameters read *\/ */
/* 	  if (flag_verbose==3) { */
/* 	    sprintf(event,"BIASSECA=%s AMPSECA=%s BIASSECB=%s AMPSECB=%s TRIMSEC=%s DATASEC=%s", */
/* 		    data.biasseca,data.ampseca,data.biassecb,data.ampsecb, */
/* 		    data.trimsec,data.datasec); */
/* 	    reportevt(flag_verbose,STATUS,1,event); */
/* 	  } */
/* 	  reportevt(flag_verbose,STATUS,1,"OVERSCAN correcting",overscantype); */
/* 	  /\* carry out overscan *\/ */
/* 	  overscan(&data,&output,flag_verbose,overscantype); */
/* 	  status=0; */
/* 	  /\* update dataset within the output image *\/ */
/* 	  sprintf(output.datasec,"[%d:%ld,%d:%ld]",1,output.axes[0], */
/* 		  1,output.axes[1]); */
/* 	} */
//	else {
//	  reportevt(flag_verbose,STATUS,1,"Already OVERSCAN corrected");
//	  flag_imoverscan=0; 
      //	}
      //      }
      //      if(!flag_imoverscan) {
      
      // OVERSCAN is disabled, so:
      /* simply copy the data into the output array  */
      /* note that this assumes input image has already been
	 overscan/trimmed*/
      /* need to make this robust-- check for trimming */
      sprintf(event,"Copying %ld pixels of image array",output.npixels);
      output.image=(float *)calloc(output.npixels,sizeof(float));
      if (output.image==NULL) {
	reportevt(flag_verbose,STATUS,5,"Calloc of output.image failed");
	exit(0);
      }
      for (i=0;i<output.npixels;i++) output.image[i]=data.image[i];
      if (data.mask!=NULL) {
	sprintf(event,"%s and mask array",event);
	output.mask=(short *)calloc(output.npixels,sizeof(short));
	if (output.mask==NULL) {
	  reportevt(flag_verbose,STATUS,5,"Calloc of output.mask failed");
	  exit(0);
	}
	for (i=0;i<output.npixels;i++) output.mask[i]=data.mask[i];
	if(flag_bpm){
	  if(flag_bpm_override){
	    sprintf(event,"%s by setting to input BPM",event);
	    for (i=0;i<output.npixels;i++) output.mask[i] =  bpm.mask[i];
	  }
	  else{
	    sprintf(event,"%s by merging image mask with input BPM",event);
	    for (i=0;i<output.npixels;i++) output.mask[i] |= bpm.mask[i];
	  }		  
	}
	else{
	  sprintf(event,"%s from input image mask",event);
	}
      }
      /* MUST FIX THIS!!!!! */
      if (data.varim!=NULL) {
	/* set variance flag */	
	flag_variance=YES;  
	sprintf(event,"%s and variance array",event);
	output.varim=(float *)calloc(output.npixels,sizeof(float));
	if (output.varim==NULL) {
	  reportevt(flag_verbose,STATUS,5,"Calloc of output.varim failed");
	  exit(0);
	}
	for (i=0;i<output.npixels;i++) output.varim[i]=data.varim[i];
      }
      reportevt(flag_verbose,STATUS,1,event);
    
    
    
      /* prepare space for variance image and initialize */
      if (flag_variance && data.varim==NULL) {
	output.varim=(float *)calloc(output.npixels,sizeof(float));
	if (output.varim==NULL) {
	  reportevt(flag_verbose,STATUS,5,"Calloc of output.varim failed");
	  exit(0);
	}
	for (i=0;i<output.npixels;i++) output.varim[i]=0.0;
      }
      /* prepare space for the bad pixel mask and initialize */
      if (data.mask==NULL) {
	reportevt(flag_verbose,STATUS,1,"Creating new output mask.");
	/* prepare image for mask == assumes mask is not yet present*/
	output.mask=(short *)calloc(output.npixels,sizeof(short));
	if (output.mask==NULL) {
	  reportevt(flag_verbose,STATUS,5,"Calloc of output.mask failed");
	  exit(0);
	}
	if (flag_bpm){ /* copy mask over */
	  reportevt(flag_verbose,STATUS,1,"Coping BPM into output mask.");
	  for (i=0;i<output.npixels;i++) output.mask[i]=bpm.mask[i];
	}
	else{  /* start with clean slate */
	  reportevt(flag_verbose,STATUS,1,"Zeroing output mask.");
	  for (i=0;i<output.npixels;i++) output.mask[i]=0;
	}
      }
      
      
      
      /* *************************************************/
      /* ******** Confirm Correction Image Sizes  ********/
      /* *************************************************/	  
      
      /* need to make sure bias and flat are same size as image */
      for (i=0;i<output.nfound;i++) {
	if (output.axes[i]!=bpm.axes[i] && flag_bpm) {
	  sprintf(event,"BPM image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  bpm.name,bpm.axes[0],bpm.axes[1],
		  output.name,output.axes[0],output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (output.axes[i]!=bias.axes[i] && flag_bias) {
	  sprintf(event,"BIAS image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  bias.name,bias.axes[0],bias.axes[1],
		  output.name,output.axes[0],output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (output.axes[i]!=flat.axes[i] && flag_flatten) {
	  sprintf(event,"FLAT image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  flat.name,flat.axes[0],flat.axes[1],
		  output.name,output.axes[0],output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (output.axes[i]!=fringe.axes[i] && flag_fringe) {
	  sprintf(event,"FRINGE image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  fringe.name,fringe.axes[0],fringe.axes[1],
		  output.name,output.axes[0],output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (output.axes[i]!=pupil.axes[i] && flag_pupil) {
	  sprintf(event,"PUPIL GHOST image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  pupil.name,pupil.axes[0],pupil.axes[1],
		  output.name,output.axes[0],output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (output.axes[i]!=photflat.axes[i] && flag_photflatten) {
	  sprintf(event,"PHOTFLAT image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  photflat.name,photflat.axes[0],photflat.axes[1],
		  output.name,output.axes[0], output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	if (output.axes[i]!=illumination.axes[i] && flag_illumination) {
	  sprintf(event,"ILLUM image %s (%ldX%ld) different size than science image %s       (%ldX%ld)",
		  illumination.name,illumination.axes[0],illumination.axes[1],
		  output.name,output.axes[0], output.axes[1]);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
      }

      /* ************************************************/
      /* ******** Report on Upcoming Processing  ********/
      /* ************************************************/	  

      if (flag_bias) {
	if (fits_read_keyword(data.fptr,"DESBIAS",comment,comment,&status)==
	    KEY_NO_EXIST) {
	  reportevt(flag_verbose,STATUS,1,"BIAS correcting");
	  flag_imbias=1;
	  status=0;
	}
	else {
	  reportevt(flag_verbose,STATUS,1,"Already BIAS corrected");
	  flag_imbias=0;
	}
      }
      if (flag_flatten) {
	if (fits_read_keyword(data.fptr,"DESFLAT",comment,comment,&status)==
	    KEY_NO_EXIST) {
	  flag_imflatten=1;
	  status=0;
	  reportevt(flag_verbose,STATUS,1,"FLAT correcting");
	}
	else {
	  reportevt(flag_verbose,STATUS,1,"Already FLAT corrected");
	  flag_imflatten=0;
	}
      }
      if (flag_pupil) {
	if (fits_read_keyword(data.fptr,"DESPUPC",comment,comment,
			      &status)==KEY_NO_EXIST) {
	  flag_impupil=1;
	  status=0;
	  reportevt(flag_verbose,STATUS,1,"PUPIL GHOST correcting");
	}
	else {
	  reportevt(flag_verbose,STATUS,1,"Already PUPIL GHOST corrected");
	  flag_impupil=0;
	}
      }
	  
      /* ************************************************/
      /* ********* Mask Missing Image Sections   ********/
      /* *********  and mark Saturated Pixels    ********/
      /* ************************************************/	  

      /* get scale value for image */
      /* note that this scalefactor is pre-bias correction */
      /* but post overscan correction */
      retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		    &scalefactor,&mode,&skysigma);
      /* cycle through image masking all pixels below threshold */
//      thresholdval=scalefactor*BADPIX_THRESHOLD;
//      maskedpixels=0;
//      for (i=0;i<output.npixels;i++) {
//      if ((output.image[i]<thresholdval || output.image[i]<0.0) &&
//          !output.mask[i]) {
//        output.mask[i] |= BADPIX_LOW; /* set the masked flag */
//	  maskedpixels++;
//      }
//    }
//    sprintf(event,"%d pixels masked with value below %.2f",
//            maskedpixels,thresholdval);
//    reportevt(flag_verbose,QA,1,event);
      /* recalculate scale factor if many pixels masked */
      if (maskedpixels>0.1*output.npixels) {
	retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		      &scalefactor,&mode,&skysigma);
      }

      /* ************************************************/
      /* ******** BIAS, Pupil Ghost and FLATTEN  ********/
      /* ************************************************/	  

      if ((flag_flatten && flag_imflatten) || (flag_bias && flag_imbias)
	  || (flag_pupil && flag_impupil)) {
	/* retrieve sky brightness in image */
	scale=1.0;offset=0.0;
	saturatepixels=0;
	accept_mask = BADPIX_SATURATE|BADPIX_INTERP|BADPIX_STAR|BADPIX_TRAIL;
	for (i=0;i<output.npixels;i++) {
	  /* grab flat field scalefactor */
	  if (flag_flatten && flag_imflatten) scale=flat.image[i];
	  /* grab bias offset */
	  if(scale <= 0){
	    output.mask[i] |= BADPIX_BPM;
	    output.varim[i] = 0;
	    output.image[i] = 0;
	  }
	  else{
	    if (flag_bias && flag_imbias) offset=bias.image[i];
	    else offset=0;
	    /* grab pupil offset */
	    if (flag_pupil && flag_impupil) 
	      offset+=pupil.image[i]*scalefactor;
	    /* do math using mask  */
	    if (!(output.mask[i])||(output.mask[i]&accept_mask)) {
	      if(column_in_section((i%output.axes[0])+1,output.ampsecan)){
		//	    if (i%output.axes[0]<=1024) { /* amplifier A */
		if (output.image[i]>=output.saturateA) {
		  output.mask[i] |= BADPIX_SATURATE;
		  saturatepixels++;
		  /* boost all saturated pixels above maxsaturate */
		  if (output.image[i]<1.5*maxsaturate) 
		    output.image[i]=1.5*(scale*maxsaturate+offset);
		}
	      }		
	      else {/* amplifier B */
		if (output.image[i]>=output.saturateB) {
		  output.mask[i] |= BADPIX_SATURATE;
		  saturatepixels++;
		  /* boost all saturated pixels above maxsaturate */
		  if (output.image[i]<1.5*maxsaturate) 
		    output.image[i]=1.5*(scale*maxsaturate+offset);
		}
	      }
	      output.image[i]=(output.image[i]-offset)/scale;
	    }
	    /* set pixel to 0.0 if it's in a bad pixel */
	    else //if(!(output.mask[i]&BADPIX_INTERP))
	      output.image[i]=0.0;
	  }
	}
	sprintf(event,"Masked %d saturated pixels",saturatepixels);
	reportevt(flag_verbose,STATUS,1,event);
      } /* end BIAS and FLATTEN */

      /* ************************************************/
      /* *************** DARK Correction ****************/
      /* ************************************************/	  

      if (flag_dark) 
	for (i=0;i<output.npixels;i++)
	  output.image[i]-=data.exptime*darkimage.image[i];

      /* end of DARK Correction */
	    

      /* ********************************************* */
      /* ********** VARIANCE Image Sections *********** */
      /* ********************************************* */


      /* NOTE:  must alter to use rdnoiseA/gainA to determine whether we are */
      /* building weight image for the first time or simply adjusting it during */
      /* an illumination correction */

      if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) {
	/* mark image structure as containing WEIGHT image */
	/* note that VARIANCE and WEIGHT are treated the same */
	output.variancetype=DES_WEIGHT;
	sprintf(event,"Creating CCD statical %s image", 
		imtypename[flag_variance]);
	reportevt(flag_verbose,STATUS,1,event);
	for (i=0;i<output.npixels;i++) 
	  //	      if ((flag_bpm && !output.mask[i]) || !flag_bpm) {
	  if (!output.mask[i]) {
	    imval = output.image[i];
	    if(imval < 0) imval = 0.0;
	    if(column_in_section((i%output.axes[0])+1,output.ampsecan)){
	      //	    if (i%output.axes[0]<=1024) { /* in AMP A section */	      
	      /* each image contributes Poisson noise, RdNoise */
	      /*  and BIAS image noise */
	      sumval=(imval/output.gainA+ 
		      Squ(output.rdnoiseA/output.gainA));
	    }
	    else { /* in AMP B */
	      /* each image contributes Poisson noise, RdNoise */
	      /*  and BIAS image noise */
	      sumval=(imval/output.gainB+
		      Squ(output.rdnoiseB/output.gainB));
	    }
	    /* apply flat field correction */
	    if (flag_flatten) sumval/=Squ(flat.image[i]);
	    if (flag_illumination) sumval/=Squ(illumination.image[i]);
	    /* add uncertainties in the bias, flat and illum corrections */
	    if (flag_bias && bias.varim[i]>0.0) sumval+=1.0/bias.varim[i];
	    if (flag_flatten && flat.varim[i]>0.0) 
	      sumval+=1.0/flat.varim[i];
	    /*
	      if (flag_illumination && illumination.varim[i]>0.0) 
	      sumval+=1.0/illumination.varim[i];
	    */
	    output.varim[i]=1.0/sumval;
	  }
	  else output.varim[i]=0.0;
      }
      else if (flag_variance==DES_SIGMA) {
	/* mark image structure as containing SIGMA image */
	output.variancetype=DES_SIGMA;
	sprintf(event,"Creating CCD statical %s image",
		imtypename[output.variancetype]);
	reportevt(flag_verbose,STATUS,1,event);
	for (i=0;i<output.npixels;i++) 
	  //	      if ((flag_bpm && !output.mask[i]) || !flag_bpm) {
	  if (!output.mask[i]) {
	    if(column_in_section((i%output.axes[0])+1,output.ampsecan)){
	      //    if (i%output.axes[0]<=1024) { /* in AMP A section */	      
	      /* each image contributes Poisson noise, RdNoise */
	      /*  and BIAS image noise */
	      sumval=(output.image[i]/output.gainA+ 
		      Squ(output.rdnoiseA/output.gainA));
	    }
	    else { /* in AMP B */
	      /* each image contributes Poisson noise, RdNoise */
	      /*  and BIAS image noise */
	      sumval=(output.image[i]/output.gainB+
		      Squ(output.rdnoiseB/output.gainB));
	    }
	    /* apply flat field correction */
	    if (flag_flatten) sumval/=Squ(flat.image[i]);
	    if (flag_illumination) sumval/=Squ(illumination.image[i]);
	    /* add uncertainties in the bias, flat and illum corrections */
	    if (flag_bias && bias.varim[i]>0.0) sumval+=1.0/bias.varim[i];
	    if (flag_flatten && flat.varim[i]>0.0) 
	      sumval+=1.0/flat.varim[i];
	    /*
	      if (flag_illumination && illumination.varim[i]>0.0) 
	      sumval+=1.0/illumination.varim[i];
	    */
	    if (sumval>1.0e-10) sumval=sqrt(sumval);
	    output.varim[i]=sumval;
	  }
	  else output.varim[i]=1.0e+10;
      }

      /* ***************************************************** */
      /* ********** remove source flux from image  *********** */
      /* ***************************************************** */

      if (flag_noisemodel==DES_SKYONLY) {
	/* create an image that has no sources */
	/* carry out setup tasks first time through */
	if (vecsort==NULL) {
	  nosource.image=(float *)calloc(output.npixels,sizeof(float));
	  if (nosource.image==NULL) {
	    reportevt(flag_verbose,STATUS,5,
		      "Calloc of nosource.image failed");
	    exit(0);
	  }
	  count=Squ(2*maxsize);
	  vecsort=(float *)calloc(count,sizeof(float));
	  if (vecsort==NULL) {
	    reportevt(flag_verbose,STATUS,5,"Calloc of vecsort failed");
	    exit(0);
	  }
	  randnum=(float *)calloc(count,sizeof(float));
	  if (randnum==NULL) {
	    reportevt(flag_verbose,STATUS,5,"Calloc of randnum failed");
	    exit(0);
	  }
	  for (i=0;i<count;i++) randnum[i]=ran1(&ranseed);
	}

	/* ************************************************************ */
	/* ********** smooth the input image to create sky image ****** */
	/* ************************************************************ */
	for (y=0;y<output.axes[1];y++) {
	  if (y%200==1 && flag_verbose) {printf(".");fflush(stdout);}
	  dy=maxsize;
	  if (y-dy<0)  dy=y;
	  if (y+dy>=output.axes[1])  dy=output.axes[1]-y-1;
	  if (dy<minsize) dy=minsize;
	  ymin=y-dy;if (ymin<0) ymin=0;
	  ymax=y+dy;if (ymax>output.axes[1]) ymax=output.axes[1];
	  for (x=0;x<output.axes[0];x++) {
	    dx=maxsize;
	    if (x-dx<0)  dx=x;
	    if (x+dx>=output.axes[0])  dx=output.axes[0]-x-1;
	    if (dx<minsize) dx=minsize;
	    xmin=x-dx;if (xmin<0) xmin=0;
	    xmax=x+dx;if (xmax>output.axes[0]) xmax=output.axes[0];
	    loc=x+y*output.axes[0];
	    /* extract median */
	    count=0;
	    ylen=ymax-ymin;
	    xlen=xmax-xmin;
	    totpix=ylen*xlen;
	    /* use all the pixels */
	    if (ylen*xlen<MAXPIXELS)
	      for (k=ymin;k<ymax;k++) for (l=xmin;l<xmax;l++)
					vecsort[count++]=output.varim[l+k*output.axes[0]];
	    else /* randomly choose the pixels */
	      while (count<MAXPIXELS) {
		k=ymin+(int)(totpix*randnum[count])/xlen;
		if (k>=ymax) k=ymax-1;
		l=xmin+(int)(totpix*randnum[count])%xlen;
		if (l>=xmax) l=xmax-1;
		vecsort[count++]=output.varim[l+k*output.axes[0]];
	      }
	    /* sort */
	    shell(count,vecsort-1);
	    /* odd or even number of pixels */
	    if (count%2) nosource.image[loc]=vecsort[count/2];
	    else nosource.image[loc]=0.5*(vecsort[count/2]+
					  vecsort[count/2-1]);
	  }
	}
	/* copy source removed image over */
	/* take care only to replace non-zero weights */
	if (output.variancetype==DES_VARIANCE || 
	    output.variancetype==DES_WEIGHT) {
	  for (i=0;i<output.npixels;i++) 
	    if (output.varim[i]>1.0e-10) output.varim[i]=nosource.image[i];
	}
	else if (output.variancetype==DES_SIGMA) {
	  for (i=0;i<output.npixels;i++) 
	    if (output.varim[i]<1.0e+10) output.varim[i]=nosource.image[i];
	}
	printf("\n");
      }
	

      /* ********************************************************** */
      /* ******** INTERPOLATE Image and Variance Section ********** */
      /* ********************************************************** */
      if (flag_interpolate_col){ // && flag_bpm) {
	if (flag_verbose==3) 
	  sprintf(event,"Interpolating over bad pixels");
	interppixels=0;
	for (i=0;i<output.npixels;i++) {
	  //	      if ((output.mask[i] == BADPIX_BPM)) {/* this is bad pixel */
	  if ((output.mask[i]&BADPIX_BPM)&&(!(output.mask[i]&BADPIX_INTERP))) {/* this is bad pixel */
	    //	  if (output.mask[i]) {/* this is bad pixel */
	    xlow=i-1;
	    while (output.mask[xlow]>0) xlow--;
	    xhi=i+1;
	    while (output.mask[xhi]>0) xhi++;
	    /* only proceed with interpolation if no more than two bad */
	    /* pixels in a row are present */
	    if (xhi-xlow<4) {
	      output.image[i]=0.5*(output.image[xlow]+output.image[xhi]);
	      /* could add noise to counteract averaging of 2 pixels */
	      /*if (flag_variance) output.image[i]+=gasdev(&seed)*
		sqrt(0.25*(1.0/output.varim[xlow]+1.0/output.varim[xhi]));*/
	      //		    output.mask[i]+=BADPIX_INTERP; /* set the interpolate flag */
	      output.mask[i] |= BADPIX_INTERP; /* set the interpolate flag */
	      interppixels++;
	      /* now calculate weight for this pixel */
	      /* allows for non-zero weight if interpolation */
	      /* over 2 or fewer pixels */
	      if (flag_variance && (xhi-xlow)<=3) { 
		if (output.varim[xlow]>1.0e-10 && output.varim[xhi]>1.0e-10)
		  output.varim[i]=0.50*(output.varim[xlow]+
					output.varim[xhi]);
		/* scale weight because it is interpolated */
		if (flag_variance==DES_WEIGHT || 
		    flag_variance==DES_VARIANCE) 
		  output.varim[i]/=Squ(scale_interpolate);
		// BUG FIX: logical check was performing assignment (mtc)
		//		else if (flag_variance=DES_SIGMA) output.varim[i]*=
		else if (flag_variance==DES_SIGMA) output.varim[i]*=
						    scale_interpolate;
	      }
	    }
	  }
	}
	if (flag_verbose==3) {
	  sprintf(event,"%s: %d pixels",event,interppixels);
	  reportevt(flag_verbose,STATUS,1,event);
	}
      }

	
      /* **************************** */
      /* **** FRINGE Correction ***** */
      /* **************************** */

      if (flag_fringe) {
	flag_imfringe=YES;
	reportevt(flag_verbose,STATUS,1,"Applying Fringe correction");
	    
	/* retrieve overall scaling from image */
	retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		      &scalefactor,&mode,&skysigma);
	for (i=0;i<output.npixels;i++) 
	  output.image[i]-=scalefactor*fringe.image[i];
      }	
	
      /* ********************************** */
      /* **** ILLUMINATION Correction ***** */
      /* ********************************** */

      if (flag_illumination) {
	flag_imillumination=YES;
	if (flag_verbose && flag_variance) 
	  reportevt(flag_verbose,STATUS,1,
		    "Applying Illumination correction to image and variance image");
	else if (flag_verbose) 
	  reportevt(flag_verbose,STATUS,1,
		    "Applying Illumination correction");
	for (i=0;i<output.npixels;i++) 
	  output.image[i]/=illumination.image[i];
        // FIXME FIXME  - need to figure out variance type here.
	if (flag_variance) {
	  for (i=0;i<output.npixels;i++) 
	    output.varim[i]/=Squ(illumination.image[i]);
	}
      }	
      
      /* *******************`******************* */
      /* ******* PhotFlatten Correction ******** */
      /* ***********************`*************** */

      if (flag_photflatten) {
	flag_imphotflatten=YES;
	if (flag_verbose) 
	  reportevt(flag_verbose,STATUS,1,
		    "Applying Photflatten correction");
	retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		      &skybrite,&mode,&skysigma);
	for (i=0;i<output.npixels;i++) 
	  output.image[i]=(output.image[i]-skybrite)*photflat.image[i]+
	    skybrite; 
      }	


      /* *************************************** */
      /* ********* CALCULATE SKYBRITE ********** */
      /* *************************************** */
      if (flag_updateskybrite || flag_updateskysigma) {
	/* retrieve overall scaling from image */
	retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		      &skybrite,&mode,&skysigma);
	sprintf(event,"Image SKYBRITE = %10.4e & SKYSIGMA = %10.4e",
		skybrite,skysigma);
	reportevt(flag_verbose,QA,1,event);
      }


      /* *********************** */
      /* **** SAVE RESULTS ***** */
      /* *********************** */

      if (flag_output) {
	/* if writing individual images for each component */
	/* then get these names now */
	if (!flag_mef) {
	  /* strip off the ".fits" */
	  sprintf(rootname,"%s",imagename);
	  for (i=strlen(rootname);i>=0;i--) 
	    if (!strncmp(rootname+i,".fits",5)) break;
	  rootname[i]=0;
	  /* create the three names */
	  sprintf(imagename,"!%s_im.fits",rootname);
	  sprintf(output.name,"%s",imagename);
	  sprintf(bpmname,"!%s_bpm.fits",rootname);
	  sprintf(varimname,"!%s_var.fits",rootname);
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
	
	/* create the (image) file */
	if (fits_create_file(&output.fptr,output.name,&status)) {
	  sprintf(event,"File creation failed: %s",output.name);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	
	/* create the image */
	if (fits_create_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
	  sprintf(event,"Image creation failed: %s",output.name);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	
	/* first read the number of keywords */
	if (fits_get_hdrspace(data.fptr,&keysexist,&j,&status)) {
	  sprintf(event,"Retrieving header failed: %s",data.name);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	sprintf(event,"Copying %d keywords into output image",
		keysexist);
	reportevt(flag_verbose,STATUS,1,event);
	/* reset to beginning of output header space */
	if (fits_read_record(output.fptr,0,keycard,&status)) {
	  reportevt(flag_verbose,STATUS,5,
		    "Reset to start of header failed");
	  printerror(status);
	}
	if (fits_read_record(data.fptr,0,keycard,&status)) {
	  reportevt(flag_verbose,STATUS,5,"Reading header record failed");
	  printerror(status);
	}
	for (j=1;j<=keysexist;j++) {
	  if (fits_read_record(data.fptr,j,keycard,&status)) {
	    sprintf(event,"Reading header record %d failed",j);
	    reportevt(flag_verbose,STATUS,5,event);
	    printerror(status);
	  }
	  if (strncmp(keycard,"BITPIX  ",8) && 
	      strncmp(keycard,"NAXIS",5)    &&
	      strncmp(keycard,"PCOUNT  ",8) &&
	      strncmp(keycard,"EXTEND  ",8) &&
	      strncmp(keycard,"GCOUNT  ",8) &&
	      strncmp(keycard,"COMMENT   FITS (Flexible Image",30) &&
	      strncmp(keycard,"COMMENT   and Astrophysics', v",30) &&
	      strncmp(keycard,"EXTNAME ",8) &&
	      strncmp(keycard,"BSCALE  ",8) &&
	      strncmp(keycard,"BZERO   ",8) 
	      ) { 
	    if (fits_write_record(output.fptr,keycard,&status)) {
	      sprintf(event,"Writing record failed: %s",keycard);
	      reportevt(flag_verbose,STATUS,5,event);
	      printerror(status);
	    }
	  }
	}
      }
      else {  /* writing on top of input image */
	output.fptr = data.fptr; 
	//	resize image
	sprintf(event,"Resizing output image %s to %ld X %ld",
		output.name+1,output.axes[0],output.axes[1]);
	reportevt(flag_verbose,STATUS,1,event);
	if (fits_resize_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
	  sprintf(event,"Resizing of image failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* change BZERO */
	if (fits_set_bscale(output.fptr,1.0,0.0,&status)) {
	  sprintf(event,"Reading BSCALE failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	output.bzero=0;
	if (fits_update_key_lng(output.fptr,"BZERO",output.bzero,comment,
				&status)) {
	  sprintf(event,"Reading BZERO failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      /* write the corrected image*/
      if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,
			 output.image,&status)) {
	sprintf(event,"Writing image failed: %s",output.name+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /* free image array */ 
      free(data.image);
      free(output.image); 

      /* Write processing history into the header */
      /* get system time */
      tm=time(NULL);
      sprintf(comment,"%s",asctime(localtime(&tm)));
      comment[strlen(comment)-1]=0;
      //      if (flag_overscan && flag_imoverscan) {
      //	if (fits_write_key_str(output.fptr,"DESOSCN",comment,
      //			       "overscan corrected",&status)) {
      //	  sprintf(event,"Writing DESOSCN failed: %s",output.name+1);
      //	  reportevt(flag_verbose,STATUS,5,event);
      //	  printerror(status);
      //	}
      /* update the DATASEC keyword */
      //      if (fits_update_key_str(output.fptr,"DATASEC",output.datasec,
      //			      "Data section within image",&status)) {
      //	sprintf(event,"Updating DATASEC failed: %s",output.name+1);
      //	reportevt(flag_verbose,STATUS,5,event);
      //	printerror(status);
      //      }
      //    }
      if (flag_bias && flag_imbias){
	if (fits_write_key_str(output.fptr,"DESBIAS",comment,
			       "bias subtracted",&status)) {
	  sprintf(event,"Writing DESBIAS failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_write_key_str(output.fptr,"BIASFIL",strip_path(bias.name),
			       "Bias image",
			       &status)) {
	  sprintf(event,"Keyword BIASFIL insert in %s failed",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      if (flag_pupil && flag_impupil){
	if (fits_write_key_str(output.fptr,"DESPUPC",comment,
			       "pupil corrected",&status)) {
	  sprintf(event,"Writing DESPUPC failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_write_key_str(output.fptr,"PUPFIL",strip_path(pupil.name),
			       "Pupil image",
			       &status)) {
	  sprintf(event,"Keyword PUPFIL insert in %s failed",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      if (flag_flatten && flag_imflatten){
	if (fits_write_key_str(output.fptr,"DESFLAT",comment,
			       "flat fielded",&status)) {
	  sprintf(event,"Writing DESFLAT failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_write_key_str(output.fptr,"FLATFIL",strip_path(flat.name),
			       "Flat image",
			       &status)) {
	  sprintf(event,"Keyword FLATFIL insert in %s failed",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      if (flag_dark){
	if (fits_write_key_str(output.fptr,"DESDARK",comment,
			       "dark corrected",&status)) {
	  sprintf(event,"Writing DESDARK failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_write_key_str(output.fptr,"DARKFIL",strip_path(darkimage.name),
			       "Dark image",
			       &status)) {
	  sprintf(event,"Keyword DARKFIL insert in %s failed",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	
      }
      if (flag_illumination && flag_imillumination){
	if (fits_write_key_str(output.fptr,"DESILLUM",comment,
			       "Illumination Corrected",&status)) {
	  sprintf(event,"Writing DESILLUM failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_write_key_str(output.fptr,"ILLUMFIL",strip_path(illumination.name),
			       "Illumination image",
			       &status)) {
	  sprintf(event,"Keyword ILLUMFIL insert in %s failed",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      if (flag_photflatten && flag_imphotflatten)
	if (fits_write_key_str(output.fptr,"DESPHOTF",comment,
			       "Photometric Flat Corrected",&status)) {
	  sprintf(event,"Writing DESPHOTF failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      if (flag_fringe && flag_imfringe){
	if (fits_write_key_str(output.fptr,"DESFRING",comment,
			       "Fringe Corrected",&status)) {
	  sprintf(event,"Writing DESFRING failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_write_key_str(output.fptr,"FRINGFIL",strip_path(fringe.name),
			       "Fringe image",
			       &status)) {
	  sprintf(event,"Keyword FRINGFIL insert in %s failed",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      if (flag_updateskysigma) 
	if (fits_update_key_flt(output.fptr,"SKYSIGMA",skysigma,7,
				"Sky Sigma (ADU)",&status)) {
	  sprintf(event,"Writing SKYSIGMA failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      if (flag_updateskybrite) 
	if (fits_update_key_flt(output.fptr,"SKYBRITE",skybrite,7,
				"Sky Brightness (ADU)",&status)) {
	  sprintf(event,"Writing SKYBRITE failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      sprintf(longcomment,"DESDM:");
      //      for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);
      sprintf(longcomment,"%s %s",longcomment,command_line);
      if (flag_verbose) reportevt(flag_verbose,STATUS,1,longcomment);
      if (fits_write_history(output.fptr,longcomment,&status)) {
	sprintf(event,"Writing longcomment failed: %s",output.name+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (fits_update_key_str(output.fptr,"DES_EXT",imtypename[DES_IMAGE],
			      "Image extension",&status)) {
	sprintf(event,"Writing DES_EXT=%s failed: %s",
		imtypename[DES_IMAGE],output.name+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (fits_update_key_lng(output.fptr,"EXTVER",1,
			      "Extension version",&status)) {
	reportevt(flag_verbose,STATUS,5,"Setting EXTVER=1 failed");
	printerror(status);
      }
	    

      /* remove unneeded information from the header */
      nkeys=0;
      while (strlen(delkeys[nkeys])) {
	if (fits_read_keyword(output.fptr,delkeys[nkeys],
			      comment,comment,&status)==KEY_NO_EXIST) status=0;
	else {
	  if (fits_delete_key(output.fptr,delkeys[nkeys],&status)) {
	    if (flag_verbose) {
	      sprintf(event,"Keyword %s not deleted from image %s",
		      delkeys[nkeys],output.name+1);
	      reportevt(flag_verbose,STATUS,3,event);
	    }
	    status=0;
	  }
	}
	nkeys++;
      }

      /* write new keywords if SATURATA/SATURATB were present */
      if (maxsaturate>0.0) 
	if (fits_write_key(output.fptr,TFLOAT,"SATURATE",&maxsaturate,
			   "Global Saturate Value",&status)) {
	  sprintf(event,"Writing keyword SATURATE failed: %s",
		  output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	  
      /* update existing keywords */
      if (fits_update_key_str(output.fptr,"OBSTYPE","red",
			      "Observation type",&status)) {
	sprintf(event,"Updating keyword OBSTYPE failed: %s",output.name+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
	    
      /* close and reopen if outputting individual images */
      if (!flag_mef) {
	/* close the image file */
	if (fits_close_file(output.fptr,&status)) {
	  sprintf(event,"Closing file failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	sprintf(event,"Opening bpm image %s",bpmname+1);
	reportevt(flag_verbose,STATUS,1,event);
	/* open the bpm file */
	if (fits_create_file(&output.fptr,bpmname,&status)) {
	  sprintf(event,"Creating image mask file failed: %s",bpmname+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }

      /* now store the bad pixel mask that has been created or updated */
      /* first create a new extension */
      if (fits_create_img(output.fptr,USHORT_IMG,2,output.axes,&status)) {
	sprintf(event,"Creating image mask failed: %s",bpmname+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* write the data */	  
      if (fits_write_img(output.fptr,TUSHORT,1,output.npixels,output.mask,
			 &status)) {
	sprintf(event,"Writing image mask failed: %s",bpmname+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* free mask image array */ 
      free(data.mask);   
      free(output.mask);   

	  
      if (fits_update_key_str(output.fptr,"DES_EXT",imtypename[DES_MASK],
			      "Extension type",&status)) {

	reportevt(flag_verbose,STATUS,5,"Setting DES_EXT=%s failed");
	printerror(status);
      }
	  
      if (fits_update_key_lng(output.fptr,"EXTVER",2,
			      "Extension version",&status)) {
	reportevt(flag_verbose,STATUS,5,"Setting EXTVER=2 failed");
	printerror(status); 
      }
	  

	  

      /* close and reopen if outputting individual images */
      if (!flag_mef) {
	/* close the image file */
	if (fits_close_file(output.fptr,&status)) {
	  sprintf(event,"Closing image failed: %s",output.name+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (flag_verbose) {
	  sprintf(event,"Opening %s image %s",
		  imtypename[output.variancetype],varimname+1);
	  reportevt(flag_verbose,STATUS,1,event);
	}
	/* open the variance image */
	if (fits_create_file(&output.fptr,varimname,&status)) {
	  sprintf(event,"Creating %s image failed: %s",
		  imtypename[output.variancetype],varimname+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }

      /* now store the variance image that has been created or updated */
      /* first create a new extension */
      if (flag_variance) {
	if (fits_create_img(output.fptr,FLOAT_IMG,2,output.axes,&status)) {
	  sprintf(event,"Creating %s image extention failed: %s",
		  imtypename[output.variancetype],varimname+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* write the data */	  
	if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,
			   output.varim,&status)) {
	  sprintf(event,"Writing %s image failed: %s",
		  imtypename[output.variancetype],varimname+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* free variance array */ 
	free(data.varim);   
	free(output.varim);   

	if (fits_update_key_str(output.fptr,"DES_EXT",
				imtypename[output.variancetype],"Extension type",&status)) {
	  sprintf(event,"Writing DES_EXT=%s failed: %s",
		  imtypename[output.variancetype],varimname+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	    
	if (fits_update_key_lng(output.fptr,"EXTVER",3,
				"Extension version",&status)) {
	  reportevt(flag_verbose,STATUS,5,"Setting EXTVER=3 failed");
	  printerror(status);
	}
	    
      }
	  	    
      /* close the corrected image */
      if (fits_close_file(output.fptr,&status)) {
	sprintf(event,"Closing image failed: %s",output.name+1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* close the input image if needed */
      if (flag_output)  
	if (fits_close_file(data.fptr,&status)) {
	  sprintf(event,"Closing input image failed: %s",data.name);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      sprintf(event,"Closed %s: 2D ( %ld X %ld )",
	      &(output.name[flag_output]),output.axes[0],output.axes[1]);
      reportevt(flag_verbose,STATUS,1,event);
    } /* end of image processing cycle with variable im */

    if (flag_list) {
      if (flag_output) 
	if (fclose(out)) {
	  sprintf(event,"Failed to close output file list: %s",outputlist);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
    }

    return (0);
}

int main(int argc,char *argv[])
{
  return(ImCorrect(argc,argv));
}
#undef MAXPIXELS


