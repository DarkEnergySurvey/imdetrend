/*--main

Basic syntax: mkflatcor <input_list> <output_flatcor_image> <options>
  Input data:
    <input_list> = ASCII file containing a list of bias FITS files 
                   to be combined
    <output_flatcor_image> = output combined flatcor FITS image
  Options:
    -average
    -avsigclip <Sigma>
    -avminmaxclip <NCLIP>
    -median (default)
    -scale <scale region>
    -noscale
    -bias <image>
    -linear <lut>
    -pupil <image>
    -bpm <image>
    -flattype <flatcor or tflatcor> 
    -variancetype <SIGMA, WEIGHT>
    -image_compare <template>
    -verbose <0-3>
    -version (print version and exit)
  Performance Options:
    -fast


Summary:
  This routine combines a list of flat field image frames from a single CCD 
  to form a master flat calibration frame. Images to be combined are provided 
  through a list (ASCII text file) with full paths to each frame.  The 
  combination option to be used is provided by a command line switch.  In 
  addition a number of other calibrations (bias, pupil, bpm, scaling) can 
  be applied prior to the actual combination of the individual frames.

Detailed Description:

  Pre-processing:
    A number of pre-processing steps can occur prior to image combination.

    Bias Correction (-bias <biascor>):
      When used, a master bias image will be read and subtracted from each 
      input frame.  The image used is FITS file with dimensions that match
      the post overscan image size (i.e. biascor image should be trimmed).

    Linearity Correction (-linear <lut>):
      This options causes a lookup table (lut) to be read which describes a
      linearity correction.  The format of this table is currently an 
      multi-extension FITS file that contains one table per CCD.  Each
      table is comprised of a list of pixel values and their associated
      value to be corrected to in order to linearize that CCDs response.
      Separate columns describe the correction for amplifiers A and B.

    BPM (-bpm <image>):
      If this options is used then a bad pixel mask (BPM) is propogated into
      the mask plane of the input image.  Note this will overwrite any mask
      that was already present within an input image.  Pixels flagged in the
      BPM are not considered when calculating scale factors.  This mask is NOT
      propogated into the output image.
     
    Pupil Correction (-pupil <image>)
      When present an image of the pupil ghost is subtracted from each of the
      input images.  This is accomplished by determining the scale (median) of
      the input image (see the -scale option below) and then subtracting a 
      scaled version of the pupil image:
                 input(pix)=input(pix)-scale*pupil(pix)

    Scaling (-scale, -noscale)
      By default, prior to combination, all images are scaled by the median 
      value in a region of the image (after all corrections have been applied)
      unless the -noscale option is present.  The -scale option allows the user 
      to specify a region to be used.  The default is to scale using the region
      [500:1500,1500:2500] where the format is: [xmin:xmax,ymin:ymax], all 
      values are integers and non-numeric characters (including ".") serve as 
      delimeters (including ".").  The same region is used when assessing the 
      scale used when applying a pupil correction.

      Since the processing flags are assessed serially at startup, the scale 
      region could in principle be updated for the application of the pupil 
      while scaling for the combined image could be turned off by specifying 
      -noscale someplace on the command line after the -scale option.

  Image Combination:
    There are currently four different options that direct how the image 
    combination is accomplished.  These are:
      -median (default)
      -average
      -avminmaxclip <Nclip>
      -avsigclip <NSigma>

    The median and average combine options perform a simple pixel-by-pixel 
    median or average for all of the input images.  The avminmaxclip and 
    avsigclip variants also perform an average on a pixel-by-pixel basis but 
    prior to this calculation, attempt to reject outlier values.  Specifically,
    the avminmaxclip option rejects the Nclip high and Nclip low values for 
    each pixel prior to calculating the average.  The avgsigclip option 
    determines the variance among the pixel values and then rejects pixels with 
    values outside the range of +/-(NSigma*variance) before proceeding to 
    calculate the average of the remaining values.

  Other Options:

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
    for each pixel) but will output the image as a flat image.
  - There is a wait/sleep command present in most routines that write images.  
    The result is that if you write to a directory other than the current 
    directory the subroutine attempts to verify the presence of the target 
    directory and if not present builds the path.  If errors are detected the 
    routine sleeps prior to a second attempt (and does so for each level of 
    subdirectory that must be checked and/or created).

*/


/* Programmer Notes follow */
/*                         */
/* Combine (dome) flat frames into a composite flat field correction */
/* author: unknown
 * modified by: Nikolay Kuropatkin 02/01/2012
 * Introduced min/max clipped image combine, implemented sigma clipped image combine,
 * rebuild variance calculation, corrected some if statements
 * N. Kuropatkin 02/06/2012.  In case of AVERAGE mode, if one image has
 * significant shift in signal the average value also will be shifted.
 * as a result calculation of variance can fail sigma limit test. In this
 * case the assigned variance will be calculated from sigma limit and not
 * be zero

 * 05/28/2013: modified by V. Kindratenko*
 *  - Added -fast option
 */

#include <float.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <getopt.h>

#include "imsupport.h"
#include "argutils.h"


//void reportevt(int flag_verbose,int type,int level,char *event);
static const char *svn_id = "$Id$";

void print_usage()
{
  printf("mkflatcor <input list> <output flatcor image> <options>\n");
  printf("  -average\n");
  printf("  -avsigclip <Sigma>\n");
  printf("  -avminmaxclip <NCLIP>\n");
  printf("  -median (default)\n");
  printf("  -scale <scale region>\n");
  printf("  -noscale\n");
  printf("  -bias <image>\n");
  printf("  -linear <lut>\n"); 
  printf("  -bpm <image>\n");
  printf("  -pupil <image>\n");
  printf("  -flattype <flatcor or tflatcor> \n");
  printf("  -variancetype <SIGMA, WEIGHT>\n");
  printf("  -image_compare <template>\n");
  printf("  -verbose <0-3>\n");
  printf("  -help (print usage and exit)\n");
  printf("  -version (print version and exit)\n");
  printf("  Performance Options\n");
  printf("    -fast\n");
}

static int flag_scale = YES;

int MakeFlatCorrection(int argc,char *argv[])
{
  int	xp,yp,i,im,imnum,x,y,loc,xmin,xmax,ymin,ymax,num,delta,ncompare=0,
    flag_combine=MEDIAN,flag_bias=NO,flag_pupil=NO,
    flag_variance=NO,flag_bpm=NO,hdutype,
    flag_variance_method=VARIANCE_DIRECT,flag_twilight=0,
    scaleregionn[4]={500,1500,1500,2500};
  int scalenum,flag_verbose=1,mkpath(),flag_image_compare=NO,overscantype=0,
    flag_aver = 0,flag_median=0,flag_avsigclip = 0,flag_avminmaxclip = 0;
  int flag_linear=NO,flag_lutinterp=YES;
  static int flag_fast=NO;
  float	sigma1,sigma2,sigma3,sigma4,
    meanval1,meanval2,meanval3,meanval4;
  int	ampregiona[4]={0,1024,0,4096},ampregionb[4]={1024,2048,0,4096}; 
  static int status=0;
  char	comment[1000],imagename[1000],command[1000],event[10000],
    longcomment[2500],scaleregion[100],filter[200]="",command_line[1000],
    linear_tab[1000],firstimage[1000],retry_firstimage[1000],
    obstype[200],*strip_path(),inname_temp[1000],outname_temp[1000],
    imtypename[6][10]={"","IMAGE","VARIANCE","MASK","SIGMA","WEIGHT"};
  fitsfile *fptr,*tmp_fptr;
  int	ccdnum=0;
  double	sumvariance,sumval,sum1,sum2a, sum2b, *flatstats;
  FILE	*inp;
  time_t	tm;
  float	gasdev(),avsigclip_sigma,*vecsort=NULL,*scaleval,*scalesort=NULL,maxval,
    val,sigmalim=0,mean,mode,fwhm,rms,offset,maxdev,
    *imagepointer, *tempimage=NULL;
    desimage *data,datain,output,bias,bpm,pupil,template;
  float image_val;
  float *lutx;
  double *luta,*lutb;
  void	rd_desimage(),shell(),decodesection(),headercheck(),
    printerror(),readimsections(),overscan(),retrievescale(),
    image_compare(),desimstat(),DataStats();
  long	seed=-15;
  /*  */
  float sigmapv, meanpv, minpv, maxpv, entries;
  int sel_pix, tailcut, loca, locb;
  int	minmaxclip_npix;
  int  ihdu,nhdunum;
  long found_ccdnum;

  /* variables for CAMSYM */

  char camsym_card[80];
  int  flag_camsym_fnd=0;

  /* variables to keep track of min/max date from nite header */
  char nite[80], mindate[80]="99999", maxdate[80]="0000", sawdate = 0;

  enum {OPT_AVERAGE=1,OPT_AVSIGCLIP,OPT_AVMINMAXCLIP,OPT_MEDIAN,OPT_SCALE,OPT_NOSCALE,
	OPT_BIAS,OPT_LINEAR,OPT_BPM,OPT_PUPIL,OPT_FLATTYPE,OPT_VARIANCETYPE,OPT_IMAGE_COMPARE,
	OPT_VERBOSE,OPT_HELP,OPT_VERSION};
  
  /*  */
  if (argc<2) {
    print_usage();
    //	  printf("  -overscantype <0-4 default=0>\n");
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

  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-verbose") || !strcmp(argv[i],"-v")) {
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
    static struct option flatcor_options[] =
      {
	{"avsigclip",     required_argument, 0,         OPT_AVSIGCLIP},
	{"avminmaxclip",  required_argument, 0,         OPT_AVMINMAXCLIP},
	{"variancetype",  required_argument, 0,         OPT_VARIANCETYPE},
	{"image_compare", required_argument, 0,         OPT_IMAGE_COMPARE},
        {"bias",          required_argument, 0,         OPT_BIAS},
        {"linear",        required_argument, 0,         OPT_LINEAR},
	{"verbose",       required_argument, 0,         OPT_VERBOSE},
	{"scale",         required_argument, 0,         OPT_SCALE},
	{"flattype",      required_argument, 0,         OPT_FLATTYPE},
	{"pupil",         required_argument, 0,         OPT_PUPIL},
        {"bpm",           required_argument, 0,         OPT_BPM},
	{"average",       no_argument,       0,         OPT_AVERAGE},
	{"median",        no_argument,       0,         OPT_MEDIAN},
	{"version",       no_argument,       0,         OPT_VERSION},
	{"help",          no_argument,       0,         OPT_HELP},
        {"noscale",       no_argument,       &flag_scale, 0},
	{"fast",          no_argument,       &flag_fast,YES},
	{0,0,0,0}
      };


    int clopx = 0;
    clop = getopt_long_only(argc,argv,"",
			    flatcor_options,&clopx);
    if(clop == -1)
      break;
    switch(clop){
    case 0:
      // For straightforward flags
      if(flatcor_options[clopx].flag != 0)
	break;
      printf("Option %s is set",flatcor_options[clopx].name);
      if(optarg)
	printf(" with %s",optarg);
      printf(".\n");
      break;
    case OPT_BPM: // -bpm
      flag_bpm=1;
      if(optarg){
	if (strncmp(&(optarg[strlen(optarg)-5]),".fits",5)
	    && strncmp(&(optarg[strlen(optarg)-8]),".fits.gz",8)) {
	  sprintf(event,"Bitmap image name required for -bpm");
	  reportevt(2,STATUS,5,event);
	  exit(0);
	}	    
	else sprintf(bpm.name,"%s",optarg);
      }
      else {
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -pupil requires an argument.");
	exit(1);
      }
      break;
    case OPT_PUPIL: // -pupil
      flag_pupil=1;
      if(optarg){
	if (strncmp(&(optarg[strlen(optarg)-5]),".fits",5)
	    && strncmp(&(optarg[strlen(optarg)-8]),".fits.gz",8)) {
	  sprintf(event,"Pupil image name required for -pupil");
	  reportevt(2,STATUS,5,event);
	  exit(0);
	}	    
	else sprintf(pupil.name,"%s",optarg);
      }
      else {
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -pupil requires an argument.");
	exit(1);
      }
      break;
    case OPT_SCALE: // -scale
      flag_scale=1;
      if(optarg){
	sprintf(scaleregion,"%s",optarg);
	decodesection(scaleregion,scaleregionn,flag_verbose); 
      }
      else {
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -scale requires an argument.");
	exit(1);
      }
      break;
    case OPT_BIAS: // -bias
      flag_bias=1;
      if(optarg){
	if (strncmp(&(optarg[strlen(optarg)-5]),".fits",5)
	    && strncmp(&(optarg[strlen(optarg)-8]),".fits.gz",8) && 
	    strncmp(&(optarg[strlen(optarg)-8]),".fits.fz",8)) {
	  sprintf(event,"Biascor image required after -bias");
	  reportevt(2,STATUS,5,event);
	  exit(0);
	}	    
	else sprintf(bias.name,"%s",optarg);
      }
      else{
	cloperr = 1;
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -bias requires an argument.");
	exit(1);
      }
      break;
    case OPT_LINEAR: // -linear
       cloperr = 0;
       flag_linear=YES;
       if(optarg){
          if (!check_image_name(optarg,CHECK_FITS,flag_verbose)) {
             cloperr = 1;
          }
          else sprintf(linear_tab,"%s",optarg);
       }
       else cloperr = 1;
       if(cloperr){
          reportevt(flag_verbose,STATUS,5,
                    "Linearity correction look-up table must follow -linear");
          command_line_errors++;
          exit(1);
       }
       break;
    case OPT_FLATTYPE: // -flattype
      if(optarg){
	if (!strcmp(optarg,"tflatcor"))  { flag_twilight=1;
	  sprintf(event,"Creating a master twilight flats");
	  reportevt(2,STATUS,1,event);
	}
      }
      else{
	command_line_errors++;
	reportevt(flag_verbose,STATUS,5,"Option -flattype requires an argument.");
	exit(1);
      }
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
    sprintf(event,"mkflatcor requires an input image list: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  else { /* expect file containing list of flats */
    imnum=0;
    inp=fopen(inname_temp,"r");
    if (inp==NULL) {
      sprintf(event,"File %s not found",inname_temp);
      reportevt(flag_verbose,STATUS,5,event);
      exit(1);
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
	exit(1);
      }
      else { /* open file and check header */
        if (imnum == 1){
           /* If this is the first image then save the name in the event that
              no other corrections besides linearity correction are requested
              (so that a file is present for the CCDNUM to be probed for)     */
           sprintf(firstimage,"%s",imagename);
        }
	if (fits_open_file(&fptr,imagename,READONLY,&status))  {
	  sprintf(event,"Input image didn't open: %s",imagename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}

        /* collect min/max nite keys*/
	if (fits_read_key_str(fptr,"NITE",nite,comment,&status) == KEY_NO_EXIST) {
          sprintf(event, "NITE keyword not found in %s, we may not get a complete {MAX,MIN}DATE range", imagename);
	  reportevt(flag_verbose,STATUS,3,event);
          status=0;
        } else {
          sawdate = 1;
          if (0 > strcmp(nite,mindate)) {
             strncpy(mindate,nite,80);
          }
          if (0 < strcmp(nite,maxdate)) {
             strncpy(maxdate,nite,80);
          }
        }

        /* check headers for each input until a card containing CAMSYM keyword is found */
        /* When found record the value for output later */
        /* Note cards are used so that this is generic of the keyword type that must be propogated */
        /* Would be nice if this operated as a loop over a set of keywords but currently this is unncessary */

        if (!flag_camsym_fnd){
           if (fits_read_card(fptr,"CAMSYM",camsym_card,&status) == KEY_NO_EXIST){
              sprintf(event, "CAMSYM keyword not found in %s", imagename);
              reportevt(flag_verbose,STATUS,3,event);
              status=0;
           }else{
              flag_camsym_fnd=1;
           }
        }

	/* confirm OBSTYPE */
	/* 	if (fits_read_key_str(fptr,"OBSTYPE",obstype,comment,&status)== */
	/* 	    KEY_NO_EXIST) { */
	/* 	  sprintf(event,"OBSTYPE keyword not found in %s",imagename); */
	/* 	  reportevt(flag_verbose,STATUS,5,event); */
	/* 	  exit(0); */
	/* 	} */
	/* 	if (strncmp(obstype,"raw_flat",5) && strncmp(obstype,"RAW_FLAT",5)) { */
	/* 	  sprintf(event,"Input images doesn't have required  OBSTYPE='raw_flat' in %s", */
	/* 		  imagename); */
	/* 	  reportevt(flag_verbose,STATUS,5,event); */
	/* 	  exit(0); */
	/* 	} */
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
    exit(0);
  }
    
    
	
  if (flag_twilight) {sprintf(obstype,"twilight flat");}
  else{sprintf(obstype,"flatcor");}
	
  /* ******************************************* */
  /* ********  READ CALIBRATION IMAGES ********* */
  /* ******************************************* */
  
  /* reopen file for processing */
  inp=fopen(inname_temp,"r");
  if (inp==NULL) {
    sprintf(event,"File not found: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
  }
  /* read in and check bias image */
  if (flag_bias) {
    rd_desimage(&bias,READONLY,flag_verbose);
    /* confirm that this is actual bias image */
    /* make sure we are in the 1st extension */
    if (fits_movabs_hdu(bias.fptr,1,&hdutype,&status)) {
      printerror(status);
      sprintf("Move to hdu=1 failed: %s",bias.name);
      reportevt(flag_verbose,STATUS,5,event);
    }
    /* check the image */
    headercheck(&bias,"NOCHECK",&ccdnum,"DESMKBCR",flag_verbose);
  }
  
  /* read in and check pupil image */
  if (flag_pupil) {
    rd_desimage(&pupil,READONLY,flag_verbose);
    /* confirm that this is actual pupil image */
    /* make sure we are in the 1st extension */
    if (fits_movabs_hdu(pupil.fptr,1,&hdutype,&status)) {
      sprintf("Move to hdu=1 failed: %s",pupil.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    headercheck(&pupil,"NOCHECK",&ccdnum,"DESMKPUP",flag_verbose);
  }
  
  /* read bad pixel mask image */
  if (flag_bpm) {
    rd_desimage(&bpm,READONLY,flag_verbose);
    /* check the image */
    headercheck(&bpm,"NOCHECK",&ccdnum,"DESMKBPM",flag_verbose);
  }	
  
 
  /* *********************************** */
  /* *****  Linearity correction  ****** */
  /* * read linearity LUT (FITS table) * */
  /* *********************************** */

  if (flag_linear){
/*   If no other corrections have been requested then the CCDNUM is not yet known
     probe the first input image to obtain CCDNUM information */
     if (ccdnum == 0){
        sprintf(event,"No calibration with CCDNUM (needed for linearity correction).");
        reportevt(flag_verbose,STATUS,1,event);
        sprintf(event,"Attempting to open first image to obtain information");
        reportevt(flag_verbose,STATUS,1,event);
        if (fits_open_file(&tmp_fptr,firstimage,mode,&status)){
           sprintf(retry_firstimage,"%s.gz",firstimage);
           status=0;
           if (fits_open_file(&tmp_fptr,retry_firstimage,mode,&status)){
              status=0;
              sprintf(retry_firstimage,"%s.fz",firstimage);
              status=0;
              if (fits_open_file(&tmp_fptr,retry_firstimage,mode,&status)){
                 status=0;
                 sprintf(event,"Failed to open first image (%s) and obtain CCDNUM",firstimage);
                 reportevt(flag_verbose,STATUS,5,event);
                 exit(1);
              }
           }
        }
        /* Determine how many extensions */
        if (fits_get_num_hdus(tmp_fptr,&nhdunum,&status)) {
           sprintf(event,"Reading HDUNUM failed: %s",firstimage);
           reportevt(flag_verbose,STATUS,5,event);
           printerror(status);
        }
        /* Search through extensions and attempt to find a CCDNUM keyword */
        ihdu=1;
        while((ihdu <= nhdunum)&&(ccdnum==0)){
           if (fits_read_key_lng(tmp_fptr,"CCDNUM",&found_ccdnum,comment,&status)==KEY_NO_EXIST){
              ihdu++;
              fits_movabs_hdu(tmp_fptr, ihdu, &hdutype, &status);
           }else{
              ccdnum=(int)found_ccdnum;
           }
        }
        if (fits_close_file(tmp_fptr,&status)) {
           sprintf(event,"Closing input image failed: %s",firstimage);
           reportevt(flag_verbose,STATUS,5,event);
           printerror(status);
        }
        if (ccdnum != 0){
           sprintf(event,"Found CCDNUM=%d",ccdnum);
           reportevt(flag_verbose,STATUS,1,event);
        }else{
           sprintf(event,"All attempts to determine the CCDNUM have failed.  Not possible to determine linearity correction.");
           reportevt(flag_verbose,STATUS,5,event);
           exit(1);
        }
     }

     /* Obtain the linearity correction LUT from the specified table of corrections */

     read_linearity_lut(linear_tab,ccdnum,&lutx,&luta,&lutb);
     sprintf(event,"Successfully read linearity correction for CCDNUM=%d",ccdnum);
     reportevt(flag_verbose,STATUS,1,event);

  }
  
  /* ******************************************* */
  /* ********  READ and PROCESS RAW FLATS ****** */
  /* ******************************************* */
  
  
  /* create an array of image structures to hold the input images */
  data=(desimage *)calloc(imnum,sizeof(desimage));
  if (data==NULL) {
    sprintf(event,"Calloc for data failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  if (flag_combine==MEDIAN || flag_combine==CLIPPEDAVERAGE ||
      flag_combine==AVSIGCLIP) {
    vecsort=(float *)calloc(imnum,sizeof(float));
    if (vecsort==NULL) {
      sprintf(event,"Calloc for vecsort failed");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
  }
  if (flag_scale) { /* prepare for medianing the scale */
    scalenum=(scaleregionn[1]-scaleregionn[0]+1)*
      (scaleregionn[3]-scaleregionn[2]+1);
    scalesort=(float *)calloc(scalenum,sizeof(float));
    if (scalesort==NULL) {
      sprintf(event,"Calloc for scalesort failed");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
  }
  scaleval=(float *)calloc(imnum,sizeof(float));	
  if (scaleval==NULL) {
    sprintf(event,"Calloc for scaleval failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
	
  /* now cycle through input images to read and prepare them */
  for (im=0;im<imnum;im++) {
    /* get next image name */
    fscanf(inp,"%s",datain.name);
    /* copy the name into the current image */
    sprintf(data[im].name,"%s",datain.name);
	  
    /* read input image */
    rd_desimage(&datain,READONLY,flag_verbose);
	
    /* retrieve basic image information from header */
    if (fits_movabs_hdu(datain.fptr,1,&hdutype,&status)) {
      sprintf("Move to hdu=1 failed: %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }

    /* check the image */
    headercheck(&datain,filter,&ccdnum,"NOCHECK",flag_verbose);

    /* get the BIASSEC information */          
    /* remove as unused (crosstalk now takes care of this)
       if (fits_read_key_str((datain.fptr),"BIASSECA",datain.biasseca,
       comment,&status)==KEY_NO_EXIST) {
       sprintf(event,"Keyword BIASSECA not defined in %s",datain.name);
       reportevt(flag_verbose,STATUS,5,event);
       printerror(status);
       }
       decodesection(datain.biasseca,datain.biassecan,flag_verbose);
    */
	  
    /* get the AMPSEC information */
    if (fits_read_key_str(datain.fptr,"AMPSECA",datain.ampseca,
			  comment,&status)==KEY_NO_EXIST) {
      sprintf(event,"Keyword AMPSECA not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    decodesection(datain.ampseca,datain.ampsecan,flag_verbose);

    /* get the BIASSEC information */
    /*
      if (fits_read_key_str(datain.fptr,"BIASSECB",datain.biassecb,
      comment,&status)==KEY_NO_EXIST) {
      sprintf(event,"Keyword BIASSECB not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
      }
      decodesection(datain.biassecb,datain.biassecbn,flag_verbose);
    */
	  
    /* get the AMPSEC information */
    if (fits_read_key_str(datain.fptr,"AMPSECB",datain.ampsecb,
			  comment,&status)==KEY_NO_EXIST) {
      sprintf(event,"Keyword AMPSECB not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    decodesection(datain.ampsecb,datain.ampsecbn,flag_verbose);

    /* get the TRIMSEC information */
    /*
      if (fits_read_key_str(datain.fptr,"TRIMSEC",datain.trimsec,
      comment,&status)==KEY_NO_EXIST) {
      sprintf(event,"Keyword TRIMSEC not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
      }
      decodesection(datain.trimsec,datain.trimsecn,flag_verbose);
    */
	  
    /* get the DATASEC information */
    /*
      if (fits_read_key_str(datain.fptr,"DATASEC",datain.datasec,
      comment,&status)==KEY_NO_EXIST) {
      sprintf(event,"Keyword DATASEC not defined in %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
      }
      decodesection(datain.datasec,datain.datasecn,flag_verbose);
    */
	  
    if (flag_verbose>=3) {
      /*
	sprintf(event,"BIASSECA=%s AMPSECA=%s BIASSECB=%s AMPSECB=%s TRIMSEC=%s DATASEC=%s",
	datain.biasseca,datain.ampseca,datain.biassecb,datain.ampsecb,
	datain.trimsec,datain.datasec);
      */
      sprintf(event,"AMPSECA=%s AMPSECB=%s",
              datain.ampseca,datain.ampsecb);
      reportevt(flag_verbose,STATUS,1,event);

    }


    /* ************************************* */
    /* ******** OVERSCANB Correction ******** */
    /* ************************************* */

    /* check to see if DESOSCN keyword is set */
    /* 	  if (fits_read_keyword(datain.fptr,"DESOSCN",comment,comment,&status) */
    /* 	    ==KEY_NO_EXIST) { */
    /* 	    status=0; */
    /* 	    if (flag_verbose==3) { */
    /*               sprintf(event,"Overscan correcting: %s",datain.name); */
    /*               reportevt(flag_verbose,STATUS,1,event); */
    /* 	    } */
    /* 	    overscan(&datain,data+im,flag_verbose,overscantype); */
    /* 	  } */
    //	  else {/* image has already been OVERSCAN corrected and trimmed */
    if (flag_verbose) {
      sprintf(event,"Already overscan corrected: %s",datain.name);
      reportevt(flag_verbose,STATUS,1,event);
    }

    /* simply copy the data into the output array  */
    sprintf(data[im].name,"%s",datain.name);
    data[im].bitpix=datain.bitpix;
    data[im].npixels=datain.npixels;
    data[im].nfound=datain.nfound;
    for (i=0;i<datain.nfound;i++) data[im].axes[i]=datain.axes[i];
    data[im].saturateA=datain.saturateA;
    data[im].saturateB=datain.saturateB;
    data[im].ampsecan[0] = datain.ampsecan[0];
    data[im].ampsecan[1] = datain.ampsecan[1];
    data[im].ampsecan[2] = datain.ampsecan[2];
    data[im].ampsecan[3] = datain.ampsecan[3];
    data[im].ampsecbn[0] = datain.ampsecbn[0];
    data[im].ampsecbn[1] = datain.ampsecbn[1];
    data[im].ampsecbn[2] = datain.ampsecbn[2];
    data[im].ampsecbn[3] = datain.ampsecbn[3];
    data[im].gainA=datain.gainA;data[im].gainB=datain.gainB;
    data[im].rdnoiseA=datain.rdnoiseA;data[im].rdnoiseB=datain.rdnoiseB;
    data[im].image=(float *)calloc(data[im].npixels,sizeof(float));
    if (data[im].image==NULL) {
      sprintf(event,"Calloc of data[im].image failed");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    for (i=0;i<data[im].npixels;i++) data[im].image[i]=datain.image[i];
    if (flag_verbose) {
      sprintf(event,"Copied input data from input image: %s",datain.name);
      reportevt(flag_verbose,STATUS,1,event);
    }
    //	  }


    /* ************************************************ */
    /* ****** Confirm Correction Image Sizes ********** */
    /* ************************************************ */

    for (i=0;i<data[im].nfound;i++) {
      if (flag_bpm && data[im].axes[i]!=bpm.axes[i] ) {
	sprintf(event,"%s (%ldX%ld) different size than flat %s (%ldX%ld)",
	        bpm.name,bpm.axes[0],bpm.axes[1],data[im].name,
		data[im].axes[0],data[im].axes[1]);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      if ( flag_bias && data[im].axes[i]!=bias.axes[i]) {
	sprintf(event,"%s (%ldX%ld) different size than flat %s (%ldX%ld)",
	        bias.name,bias.axes[0],bias.axes[1],data[im].name,
		data[im].axes[0],data[im].axes[1]);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      if (flag_pupil && data[im].axes[i]!=pupil.axes[i] ) {
	sprintf(event,"%s (%ldX%ld) different size than pupil %s (%ldX%ld)",
	        pupil.name,pupil.axes[0],pupil.axes[1],data[im].name,
		data[im].axes[0],data[im].axes[1]);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
    }


    /* ***************************************** */
    /* ********* Copy BPM into place *********** */
    /* ***************************************** */

    data[im].mask=(short *)calloc(data[im].npixels,sizeof(short));
    if (flag_bpm) { /* copy the bpm into the image mask array */
      if (flag_verbose) {
        sprintf(event,"Copying input bpm from: %s",bpm.name);
        reportevt(flag_verbose,STATUS,1,event);
      }
      if (data[im].mask==NULL) {
        sprintf(event,"Calloc of data[im].mask failed");
        reportevt(flag_verbose,STATUS,5,event);
        exit(0);
      }
      for (i=0;i<data[im].npixels;i++) {
      /*  data[im].mask[i]=bpm.mask[i]; */

      /* BPM mask to supplement (logically or) existing mask) */
        data[im].mask[i]=0;
        if (bpm.mask[i]){
          data[im].mask[i]|=BADPIX_BPM;
        }
      }
      reportevt(flag_verbose,STATUS,1,"Set badpixel mask bit to BADPIX_BPM where bpm != 0");
    }

    /* ****************************************************************** */
    /* ************ prepare output image first time through ************* */
    /* ****************************************************************** */
    if (im==0) { 
      if (flag_verbose) {
	reportevt(flag_verbose,STATUS,1,"Preparing output image.");
      }
      output=data[im];
      for (i=0;i<output.nfound;i++) output.axes[i]=data[im].axes[i];
      output.bitpix=FLOAT_IMG;
      output.npixels=data[im].npixels;
      sprintf(output.name,"!%s",outname_temp);
      output.image=(float *)calloc(output.npixels,sizeof(float));
      if (output.image==NULL) {
	sprintf(event,"Calloc failed for output.image");
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      if (flag_variance) {
	/* prepares variance image */
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
      }
      else output.varim=NULL;
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

	  
    /* ************************************* */
    /* ********* BIAS Correction *********** */
    /* ************************************* */
    if (flag_bias) {
      /* check to see if DESBIAS keyword is set */
      if (fits_read_keyword(datain.fptr,"DESBIAS",comment,comment,
			    &status)==KEY_NO_EXIST) {
        status=0;
        if (flag_verbose>=3) {
          sprintf(event,"Bias subtracting");
          reportevt(flag_verbose,STATUS,1,event);
        }

        /* subtract all pixels, even those flags on tje bpm */
        /* if ( flag_bpm) {
          for (i=0;i<data[im].npixels;i++)
            if (!bpm.mask[i]) data[im].image[i]-=bias.image[i]; 
            } else */

          /* subtract all pixels */

          for (i=0;i<data[im].npixels;i++)
            data[im].image[i]-=bias.image[i];
      }
      /*else */
      if (flag_verbose>=3)
        reportevt(flag_verbose,STATUS,1,"Already bias subtracted");
    } /* end if flag_bias */

    /* ************************************* */
    /* ****** Linearity Correction ********* */
    /* ************************************* */

    if (flag_linear){
       for (i=0;i<data[im].npixels;i++){
          image_val=data[im].image[i];
          if (column_in_section((i%data[im].axes[0])+1,data[im].ampsecan)){
             image_val=lut_srch(image_val,luta,flag_lutinterp);
          }else{
             image_val=lut_srch(image_val,lutb,flag_lutinterp);
          }
          data[im].image[i]=image_val;
       }
    }

    /* ************************************* */
    /* ****** PUPIL GHOST Correction ******* */
    /* ************************************* */

    if (flag_pupil) {
      if (fits_read_keyword(datain.fptr,"DESPUPC",comment,comment,
			    &status)==KEY_NO_EXIST) {
	status=0;
	if (flag_verbose==3) {
	  sprintf(event,"Correcting for Pupil Ghost");
	  reportevt(flag_verbose,STATUS,1,event);
	}
	/* get scale value for image */
       if (flag_fast)
	 retrievescale_fast(data+im,scaleregionn,scalesort,flag_verbose,
		      scaleval+im,&mode,&fwhm); 
       else
	 retrievescale(data+im,scaleregionn,scalesort,flag_verbose,
		      scaleval+im,&mode,&fwhm); 
	/* apply correction */
	for (i=0;i<data[im].npixels;i++)
	  data[im].image[i]-=pupil.image[i]*scaleval[im];
      }
      else if (flag_verbose==3)
	reportevt(flag_verbose,STATUS,1,"Already Pupil Ghost corrected");
    }


    /* ************************************* */
    /* ****** Calculate Image Scale  ******* */
    /* ************************************* */
    /*printf("imagename=%s",data[im].name);*/
    if (flag_scale) { 
      if (flag_fast)
        retrievescale_fast(data+im,scaleregionn,scalesort,flag_verbose,
		    scaleval+im,&mode,&fwhm);
      else
        retrievescale(data+im,scaleregionn,scalesort,flag_verbose,
		    scaleval+im,&mode,&fwhm); 
      if (flag_twilight){
	for (i=1;i<data[im].npixels;i++) {   
	  if (data[im].image[i]>1*fwhm+scaleval[im]/*||data[im].image[i]<-3*fwhm+scaleval[im]*/) {
	    xp = i%data[im].axes[0];
	    yp = i/output.axes[0];
		  
	    /*sprintf(event,"The adjusted,median,sd, x and y values  are %f %lf,and %lf, %d,%d",
	      data[im].image[i],mode,fwhm,xp,yp);
	      reportevt(flag_verbose,STATUS,1,event);*/
	    data[im].image[i]=scaleval[im]+fwhm*gasdev(&seed);
		  
	  }}}
      /*exit(0);*/
      if (flag_verbose) {
	sprintf(event,"Image=%s & Region = [%d:%d,%d:%d] & Scale= %f",
		data[im].name,scaleregionn[0],scaleregionn[1],scaleregionn[2],scaleregionn[3],
		scaleval[im]);
	reportevt(flag_verbose,QA,1,event);
      }
    }
    else scaleval[im]=1.0;
	  
	  		  

    /* ************************************* */
    /* ******** FITS Housekeeping  ********* */
    /* ************************************* */

    /* close input FITS file before moving on */
    if (fits_close_file(datain.fptr,&status)) {
      sprintf(event,"File close failed: %s",datain.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    free(datain.image);
  } /* all input images have been read and overscan+bias subtracted */
  /* close input file list */	
  if (fclose(inp)) {
    sprintf(event,"Closing image list failed: %s",inname_temp);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }

	
  /* ******************************************* */
  /* *********  COMBINE SECTION *************** */
  /* ******************************************* */
	
  if (flag_verbose==3) 
    reportevt(flag_verbose,STATUS,1,"Combining images");

  /* Calculate sigma limits to calculate variance */
  if (flag_variance) {
    sum1=0;
    for (im=0;im<imnum;im++) {
      /* calculate average signal */
      sum2a=0.;  sum2b = 0.; num = 0;
      for (y=1500;y<2500;y++) for (x=400;x<600;x++) {
	loca=y*data[im].axes[0]+x;
	locb=y*data[im].axes[0]+x+1024;
	if(column_in_section(x+1,data[im].ampsecan)){  // loca in Amp A or Amp B
	  val=data[im].image[loca]/data[im].gainA; /* counts in electrons */
	  sum2a+=val;
	  val=data[im].image[locb]/data[im].gainB;
	  sum2b+=val;
	  num++;
	}
	else {
	  val=data[im].image[loca]/data[im].gainB; /* counts in electrons */
	  sum2b+=val;
	  val=data[im].image[locb]/data[im].gainA;
	  sum2a+=val;
	  num++;
	}
		    
      }
      sum1 += 0.5*(sqrt(sum2a/(float)num) + sqrt(sum2b/(float)num))/scaleval[im];
    }
	    
    sigmalim=3.5*sum1/(float)imnum; /* this is average sigma */
	    
    /* introduce lower limit variance for clipping */
    if (sigmalim<3.5*1.0e-02) sigmalim=3.5*1.0e-02;
	    
    sprintf(event,"Calculated sigma= %f sigmalim= %f",sum1/(float)imnum, sigmalim);
    reportevt(flag_verbose,STATUS,1,event);
  } /* end of VARIANCE_CCD section */
	
  /* Combine the input images to create composite */
  if (flag_combine==AVERAGE) {


    for (i=0;i<output.npixels;i++) {
      output.image[i]=0;
      /* */
      /* Change R. Covarrubias: Don't use bpm when combining images. Using the bpm will zero the pixels where a bad pixels exists */
      /* if (!flag_bpm || !bpm.mask[i] ) { */
      for (im=0;im<imnum;im++)
        output.image[i]+=data[im].image[i]/scaleval[im];
      /* } */
      output.image[i]/=(float)imnum;


      /* CALCULATE VARIANCE  */
      if (flag_variance) {
        sel_pix = 0;
        sigmapv = 0.0;
        tempimage[i] = 0.0;

        for (im=0;im<imnum;im++){

          val = data[im].image[i]/scaleval[im]-output.image[i];
          if (fabs(val)<sigmalim ) {

            sigmapv += Squ(val);
            sel_pix++;
          }
        }
        if (sel_pix >= 2) {
          sigmapv = sigmapv/((float)sel_pix-1.)/((float)sel_pix);  /* variance */
        }
        tempimage[i]=sigmapv;
        if(sel_pix == 0) tempimage[i]= Squ(sigmalim);
      }    /* end if variance  */
      /*} else {
        output.image[i] = 0.0;
        if (flag_variance) tempimage[i]= Squ(sigmalim);
        } */
    } /* end loop on pixels */

  } /* end if average */
	
  if (flag_combine==MEDIAN) {
    for (i=0;i<output.npixels;i++) { /* for each pixel */
      output.image[i]=0.0f;

      /* */
      /* Change R. Covarrubias: Don't use bpm when combining images. Using the bpm will zero the pixels where a bad pixels exists */      
      /*if (!flag_bpm || !bpm.mask[i] ) { */
	/* copy values into sorting vector */
        for (im=0;im<imnum;im++) vecsort[im]=data[im].image[i]/scaleval[im];

        if (flag_fast)
          output.image[i] = quick_select(vecsort, imnum);
        else
          {
            shell(imnum,vecsort-1);
            /* odd number of images */
            if (imnum%2) output.image[i]=vecsort[imnum/2]; 
            /* record the median value */  
            else output.image[i]=0.5*(vecsort[imnum/2]+vecsort[imnum/2-1]);
          }

        /* CALCULATE VARIANCE  */
        if (flag_variance) {
          sel_pix = 0;
          sigmapv = 0.0;
          tempimage[i] = 0.0;

          for (im=0;im<imnum;im++){

            val = vecsort[im]-output.image[i];
            if (fabs(val)<sigmalim ) {

              sigmapv += Squ(val);
              sel_pix++;
            }
          }
          if (sel_pix >= 2) {
            sigmapv = sigmapv/((float)sel_pix-1.)/((float)sel_pix);  /* variance */
          }
          tempimage[i]=sigmapv;
          if(sel_pix == 0) tempimage[i]= Squ(sigmalim);
        } /* end if variance */

        /* } End of if (!flag_bpm || !bpm.mask[i] */ 
      /* This else is used when the bpm is taken into account */
      /* else { */ 
        /* output.image[i] = 0.0f; */
        /* if (flag_variance) tempimage[i] = Squ(sigmalim);*/
        /* } */
    }
  }
  /* **************************************************** */
  /* ** here should go implementation of other methods ** */
  /* **************************************************** */

  /* ********************************************* */
  /* ************ Implementing clipped average *** */
  /* ******* By N. Kuropatkin   01/20/2012 ******* */
  /* ********************************************* */

  if (flag_combine==CLIPPEDAVERAGE) {
    for (i=0;i<output.npixels;i++) {
      output.image[i]=0;
      /*if (!flag_bpm || !bpm.mask[i] ) { */
        for (im=0;im<imnum;im++) vecsort[im]=data[im].image[i]/scaleval[im];
        if (flag_fast)
          float_qsort(vecsort, imnum);
        else
          shell(imnum,vecsort-1);            /* sort the vector */
        if (2*minmaxclip_npix >=imnum - 1) {
          sprintf(event,"CLIPPEDAVERAGE minmaxclip_npix can not be more than imnum/2 Reset to 1");
          minmaxclip_npix = 1;
          reportevt(flag_verbose,STATUS,1,event);
        }
        for (im=minmaxclip_npix ;im<imnum-minmaxclip_npix; im++) output.image[i]+=vecsort[im];
        output.image[i]/=(float)(imnum - 2*minmaxclip_npix );
        /* CALCULATE VARIANCE  */
        if (flag_variance) {
          sel_pix = 0;
          sigmapv = 0.0;
          tempimage[i] = 0.0;
          for (im=minmaxclip_npix;im<imnum-minmaxclip_npix;im++){
            val = vecsort[im]-output.image[i];
            if (fabs(val)<sigmalim ) {
              sigmapv += Squ(val);
              sel_pix++;
            }
          }
          if (sel_pix >= 2) {
            sigmapv = sigmapv/((float)sel_pix-1.)/((float)sel_pix);  /* variance */
          }
          tempimage[i]=sigmapv;
          if(sel_pix == 0) tempimage[i]= Squ(sigmalim);
        }
        /* } end of if (!flag_bpm || !bpm.mask[i] ) { */ 
      /* Changing this such that masked bit are not set to zero */
      /* else { */  /* if masked bit let set it to 0 */
      /*output.image[i] = 0; */
      /* if (flag_variance) tempimage[i]= Squ(sigmalim);*/
        /*}*/
    } /* end of for loop */
  } /* end of clipped average */

  /* ********************************************* */
  /* ****** Implementing sigma clipped average *** */
  /* **** By N. Kuropatkin 01/20/2012   ********** */
  /* ********************************************* */

  if (flag_combine==AVSIGCLIP) {
    if ( imnum < 4) {
      sprintf(event,"Sigma clipping require at least 4 images to combine");
      reportevt(flag_verbose,STATUS,5,event);
      exit(1);
    }
    for (i=0;i<output.npixels;i++) { /* for each pixel */
      output.image[i]=0;
      sigmapv = 0.0;
      meanpv = 0.0;
      tempimage[i] = 0.;
      /*if (!flag_bpm || !bpm.mask[i] ) { */ 
      /* first cut tails of the distribution to estimate sigma */
        tailcut = (int) (0.15*imnum);  /* cut tails 30% to estimate sigma */
        if (tailcut < 1) tailcut = 1;  /* this requires at list 4 images to combine */
        for (im=0;im<imnum;im++)vecsort[im]=data[im].image[i]/scaleval[im];
        if (flag_fast)
          float_qsort(vecsort, imnum);
        else
          shell(imnum,vecsort-1);
        sigmapv = 0.5*(vecsort[imnum -1 -tailcut] - vecsort[tailcut]);
        if (imnum%2) meanpv=vecsort[imnum/2]; /* record the median value */
        else meanpv=0.5*(vecsort[imnum/2]+vecsort[imnum/2-1]);
        /* now calculate average within selected sigma cut */
        minpv = meanpv - avsigclip_sigma*sigmapv;
        maxpv = meanpv + avsigclip_sigma*sigmapv;
        if (minpv < vecsort[0]) minpv = vecsort[0];
        if (maxpv > vecsort[imnum-1]) maxpv = vecsort[imnum - 1];
        sel_pix = 0;
        for (im=0;im<imnum;im++){
          if(vecsort[im] >= minpv && vecsort[im] <= maxpv) {
            output.image[i] += vecsort[im];
            sel_pix++;
          }
        }
        output.image[i] /= (float)sel_pix;
      /* CALCULATE VARIANCE  */
        if (flag_variance) {
          sel_pix = 0;
          sigmapv = 0.0;
          tempimage[i] = 0.0;
          for (im=0;im<imnum;im++){
            if((vecsort[im] >= minpv) && (vecsort[im] <= maxpv)) {
              val = vecsort[im]-output.image[i];
              if (fabs(val)<sigmalim ) {
                sigmapv += Squ(val);
                sel_pix++;
              }
            }
          }
          if (sel_pix >= 2) {
            sigmapv = sigmapv/((float)sel_pix-1.)/((float)sel_pix);  /* variance */
          }
          tempimage[i]=sigmapv;
          if(sel_pix == 0) tempimage[i]= Squ(sigmalim);
        }
        /*}  else {end of if (!flag_bpm || !bpm.mask[i] ) { */
        /* Changing this such that masked bit are not set to zero */
        /* output.image[i] = 0; */
        /*if (flag_variance) tempimage[i]=  Squ(sigmalim); */
        /* } */
    } /*end of for loop */
  } /* end of avsigclip */
    

  /* now let's confirm the final scaling */
  if (flag_scale) {
    if (flag_fast)
      retrievescale_fast(&output,scaleregionn,scalesort,flag_verbose,
		  &maxval,&mode,&fwhm);
    else
      retrievescale(&output,scaleregionn,scalesort,flag_verbose,
		  &maxval,&mode,&fwhm); 
    if (flag_verbose) {
      sprintf(event,"Image=%s & Region = [%d:%d,%d:%d] & Scale=%f calcimscale",
	      output.name, scaleregionn[0],scaleregionn[1],scaleregionn[2],scaleregionn[3],
	      maxval);
      reportevt(flag_verbose,QA,1,event);
    }
  }
  /*skip this part for twilight sky flat */
  if (!flag_twilight)
    {  
      if (flag_verbose) {
	sprintf(event,"Rescaling output image in region [%d:%d,%d:%d] of %.5f to 1.0",
		scaleregionn[0],scaleregionn[1],scaleregionn[2],scaleregionn[3],
		maxval);
	reportevt(flag_verbose,STATUS,1,event);
      }
      /* cycle through resulting image to normalize it to 1.0 */
      for (i=0;i<output.npixels;i++) output.image[i]/=maxval;
    }
	
	
  /* ******************************************* */
  /* *********  VARIANCE SECTION *************** */
  /* ******************************************* */
  /* This section was rebuild by N. Kuropatkin   */
  /* as we already have calculated variances we  */
  /* just need to average them                   */
  /* now prepare direct variance image for the flat field */
  if (flag_variance_method==VARIANCE_DIRECT && flag_variance) { 
    if (flag_verbose==3) reportevt(flag_verbose,STATUS,1,
				   "Beginning direct variance calculation");
    /*     */
    /*     */

    delta=VARIANCE_DELTAPIXEL;
	  
	  
    if (flag_verbose) {
      sprintf(event,"Extracting variance using %dX%d square centered on each pixel",
              2*delta+1,2*delta+1);
      reportevt(flag_verbose,STATUS,1,event);
    }
    for (i=0;i<output.npixels;i++) {
      if ( !flag_bpm || !bpm.mask[i] ) {
        x=i%output.axes[0];y=i/output.axes[0];
        /* define a small square centered on this pixel */
	    
        xmin=x-delta;xmax=x+delta;
        if (xmin<0) {xmax+=abs(xmin);xmin=0;} /* hanging off bottom edge? */
        if (xmax>=output.axes[0]) /* hanging off top edge? */
          {xmin-=(xmax-output.axes[0]+1);xmax=output.axes[0]-1;}    
	    
        ymin=y-delta;ymax=y+delta;
        if (ymin<0) {ymax+=abs(ymin);ymin=0;} /* hanging off bottom edge? */
        if (ymax>=output.axes[1])  /* hanging off top edge? */
          {ymin-=(ymax-output.axes[1]+1);ymax=output.axes[1]-1;}

        sum2a=sum1=0.0;num=0;
        for (x=xmin;x<=xmax;x++) for (y=ymin;y<=ymax;y++) {
          loc=x+y*output.axes[0];  
          for (im=0;im<imnum;im++) {
            val=tempimage[loc];  /* variance */
            sum1+=val;
            num++;

          }
        }
        if (num>1) {
          val=sum1/(float)num;
        }
        else {
          sprintf(event,"No pixels meet criterion for inclusion");
          reportevt(flag_verbose,STATUS,3,event);
          val=Squ(sigmalim/3.5/(float)imnum);
        }
        if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) {
          /* introduce lower limit THRESHOLD on VARIANCE */
          if (val>1.0e-6) output.varim[i]=1.0/val;
          else output.varim[i]=1.0e+6;
        }
        else if (flag_variance==DES_SIGMA) {
          if (val>1.0e-10) val=sqrt(val);
          output.varim[i]=val;
        }
      }
      else output.varim[i]=0.0;
    }
  } /* end of VARIANCE_DIRECT section */	
  else if (flag_variance_method==VARIANCE_CCD && flag_variance) {
    if (flag_verbose) 
      reportevt(flag_verbose,STATUS,1,
		"Beginning CCD statical variance calculation");
    for (i=0;i<output.npixels;i++) 
      if ( !flag_bpm || (flag_bpm && !bpm.mask[i])) {
	sumval=sumvariance=0.0;
	for (im=0;im<imnum;im++) {
	  if(column_in_section((i%output.axes[0])+1,data[im].ampsecan)){
	    //	        if (i%output.axes[0]<1024) { /* in AMP A section */	      
	    /* each image contributes Poisson noise, */
	    /* RdNoise and BIAS image noise */
	    sumval=(data[im].image[i]/data[im].gainA+
		    Squ(data[im].rdnoiseA/data[im].gainA)+1.0/bias.varim[i])
	      /Squ(scaleval[im]);
	  }
	  else { /* in AMP B */
	    /* each image contributes Poisson noise, */
	    /* RdNoise and BIAS image noise */
	    sumval=(data[im].image[i]/data[im].gainB+
		    Squ(data[im].rdnoiseB/data[im].gainB)+1.0/bias.varim[i])/
	      Squ(scaleval[im]);		
	  }
	  sumvariance+=(sumval/Squ(maxval));
	}
	sumvariance/=Squ((float)imnum);
	if (flag_variance==DES_VARIANCE || flag_variance==DES_WEIGHT) {
	  if (sumvariance>1.0e-6) output.varim[i]=1.0/sumvariance;
	  else output.varim[i]=1.0e+6;
	}
	else if (flag_variance==DES_SIGMA) {
	  if (sumvariance>1.0e-6) sumvariance=sqrt(sumvariance);
	  output.varim[i]=sumvariance;
	}
      }
      else output.varim[i]=0.0;
  } /* end of VARIANCE_CCD section */

		
  /* ************************************************************ */
  /* ************** Test image against template ***************** */
  /* ************************************************************ */
  if (flag_image_compare) {
    /*  Read template image */
    rd_desimage(&template,READONLY,flag_verbose);
    /* check image for ccdnumber */
    headercheck(&template,"NOCHECK",&ccdnum,"DESMKFCR",flag_verbose);
    /* first compare images */
    offset=rms=maxdev=0.0;
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
      offset=rms=maxdev=0.0;
      ncompare=0;
      image_compare(&output,&template,&offset,&rms,&maxdev,&ncompare,
		    flag_verbose);
      output.image=imagepointer;
      /* issue STATUS events according to differences measured */

      /* ***************************************************** */
      /* ***************************************************** */
      /* ***************************************************** */
    }
  }


  /* ****************************************************************** */
  /* ********************* Save flatcor image ************************* */
  /* ****************************************************************** */
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
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,
		     output.image,&status)) {
    sprintf(event,"Image writing failed: %s",output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  /* output min/max nite date info */
  if (sawdate) {
    if (fits_write_key_str(output.fptr,"MINNITE",mindate,
			   "Earliest NITE covered",&status)) {
      sprintf(event,"Writing MINNITE=%s failed: %s",mindate,output.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);  
    }
    if (fits_write_key_str(output.fptr,"MAXNITE",maxdate,
			   "Latest NITE covered",&status)) {
      sprintf(event,"Writing MAXNITE=%s failed: %s",maxdate,output.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);  
    }
  }

  /* Note cards are used so that this is generic of the keyword type that must be propogated */
  /* Would be nice if this operated as a loop over a set of keywords but currently this is unncessary */
  /* Currently only CAMSYM is needed extra structure not yet necessary */

  if (flag_camsym_fnd){
     if (fits_update_card(output.fptr,"CAMSYM",camsym_card,&status)){
        sprintf(event, "Writing/updating %s failed in :%s",camsym_card,output.name+1);
        reportevt(flag_verbose,STATUS,5,event);
        printerror(status);
     }
  }else{
     sprintf(event, "No CAMSYM card found among inputs.");
     reportevt(flag_verbose,STATUS,3,event);
  }
	  
  /* write basic information into the header */
  if (fits_write_key_str(output.fptr,"OBSTYPE",obstype,
			 "Observation type",&status)) {
    sprintf(event,"Writing OBSTYPE=%s failed: %s",obstype,output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);  
  }
  if (fits_write_key_str(output.fptr,"FILTER",filter,
			 "Filter name(s)",&status)) {
    sprintf(event,"Writing FILTER=%s failed: %s",filter,output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (fits_write_key_lng(output.fptr,"CCDNUM",ccdnum,
			 "CCD Number",&status)) {
    sprintf(event,"Writing CCDNUM=%d failed: %s",ccdnum,output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  /* Write information into the header describing the processing */
  flatstats = (double *)calloc(4,sizeof(double));
  DataStats(scaleval,flatstats,imnum);
  tm=time(NULL);
  sprintf(comment,"%s",asctime(localtime(&tm)));
  comment[strlen(comment)-1]=0;
  if (fits_write_key_str(output.fptr,"DESMKFCR",comment,"mkflatcor output",&status)) {
    sprintf(event,"Writing DESMKFCR %s failed: %s",comment,output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if(fits_write_key_dbl(output.fptr,"SCALMEAN",flatstats[2],-6,"Mean Scale Value",&status)){
    sprintf(event,"Writing SCALMEAN %f failed: %s",flatstats[2],output.name+1);
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
    sprintf(event,"Writing mkflatcor call to history failed: %s",&(output.name[1]));
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
    sprintf(event,"Writing mkflatcor image list failed: %s",&(output.name[1]));
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
      sprintf(event,"Writing mkflatcor input image failed: %s",&(output.name[1]));
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
      sprintf(event,"Creating %s extension failed: %s",
	      imtypename[flag_variance],output.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* write the data */	  
    if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.varim,
		       &status)) {
      sprintf(event,"Writing %s image failed: %s",
	      imtypename[flag_variance],output.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (fits_write_key_str(output.fptr,"DES_EXT",
			   imtypename[flag_variance],"Extension type",&status)) {
      sprintf(event,"Writing DES_EXT=%s failed: %s",
	      imtypename[flag_variance],output.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }


  /* close the corrected image */
  if (fits_close_file(output.fptr,&status)) {
    sprintf(event,"Image close failed: %s",output.name+1);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (flag_verbose) {
    sprintf(event,"Closed %s: 2D ( %ld X %ld )",
            output.name+1,output.axes[0],output.axes[1]);
    reportevt(flag_verbose,STATUS,1,event);
  }

	
  free(output.image);
  for (im=0;im<imnum;im++) free(data[im].image);
  if (flag_bias) {
    free(bias.image);
    if (bias.varim!=NULL) free(bias.varim);
  }
  if (flag_bpm) free(bpm.mask);

  return(0);
}

void DataStats(float *data,double *stats,int ndat)
{
  stats[0] = DBL_MAX;  // min
  stats[1] = -DBL_MAX; // max
  stats[2] = 0.0;      // mean
  stats[3] = 0.0;      // sigma
  int i = 0;
  for(i = 0;i < ndat;i++){
    if(data[i] < stats[0]) stats[0] = data[i];
    if(data[i] > stats[1]) stats[1] = data[i];
    stats[2] += data[i];
    stats[3] += (data[i]*data[i]);
  }
  stats[2] /= ndat;
  stats[3] /= ndat;
  stats[3] = sqrt(stats[3] - (stats[2]*stats[2]));
}
int main(int argc,char *argv[])
{
  return(MakeFlatCorrection(argc,argv));
}
