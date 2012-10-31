/*--main

Basic syntax: DECam_crosstalk <infile.fits> <outfile> <options>
  Input data:
    <infile.fits>  = raw DECam FITS image (i.e. w/ one HDU per chip)
    <outfile> = output file name in the form <pref>_%02d_<suf>[.fits]
                if %02d is present, it will be substituted with CCDN
                if it is absent, CCDN.fits will be attached to the end of the filename
  Options:
    -crosstalk <crosstalk matrix- file>
    -photflag <0 or 1>
    -satmask

    -overscan
    -overscansample <-1 for MEDIAN, 0 for MEAN, 1 for MEAN w/MINMAX>,
    -overscanfunction < -50 < N < -1 for cubic SPLINE, 0 for LINE_BY_LINE, 
                          1 < N < 50 for legendre polynomial>
    -overscanorder <1-6>, default=1 <order for legendre polynomial>
    -overscantrim <ncols>, default=0 

    -maxhdunum <maxhdunum>, default=70 (all images)
    -ccdlist <comma-separated list of CCDs to process, 1-70 range>
    -hdulist <comma-separated list of HDUs to process, 2-71 range>
    -focuschipsout
    -verbose <0-3>

Summary:
  DECam_crosstalk performs up to five functions depending on how it is called.
    - breaks up a raw DECam image (currently referred to as a "src" file) 
      into its individual components (i.e., splits the HDUs from the incoming 
      file into individual FITS images (one per CCD)
    - applies a crosstalk correction between CCDs/AMPs...
    - identifies saturated pixels (prior to any image detrending/manipulation)
    - characterizes and removes overscan (instrumental bias)
    - sets initial value for the photometricity flag

Detailed Description:

  Image Split:
    DECam_crosstalk receives raw DECam data comprised of 70 HDUs (one per ccd).
    Currently the code separates the file into its constituent HDUs keeping 
    track of the images based upon the CCDNUM keyword (necessary because SISPI 
    does not necessarily keep the HDUs in a fixed order). The default behavior
    for the output filename/path is to use the argument from -output as a 
    preamble which is appended with "_<CCDNUM>.fits" for each output file 
    (one per CCD).  

    A more varied output pattern can be created by specifying a filename/path
    that includes the string "%02d" within the output argument.  The behavior 
    for such an argument is to replace the %02d with the CCDNUM.  Therefore, 
    -output filepath/preamble%02dsuffix will write the HDU with CCDNUM=09 as:
    filepath/preamble09suffix.fits .

    Additional options are present that can be used to process specific HDUs 
    and/or CCDNUMs.  These are:

    -focuschipsout 
       When enabled, images with CCDNUM > 62 (the focus chips) are written.

    -ccdlist XX,YY,ZZ
       Here the argument is a comma separated list of CCD numbers that will 
       cause those CCDs to be written (but not others).

    -hdulist AA,BB,CC
       Here the argument is a comma separated list of HDU numbers that will 
       cause those CCDs to be written that were contained within the specified
       HDUs.

    -maxhdunum N
       When this argument is provided then only the first N HDU's are written.

    Finally, note that this program also alters the OBSTYPE header keyword to 
    conform to DESDM standards, and leaves a history in the header.

  Crosstalk Correction:
    If the -crosstalk option is used then the code will read the crosstalk 
    matrix file to determine the coefficients that are applied.  This file 
    currently has a format similar to that used in IRAF for the MOSAIC cameras.
    Specifically the values are:

      victim   source
      ccdXXA   ccdXXA  coefficient (unc,sigma)

    where the names for the victim/source declaration have XX=ccdnum and 
    A=amplifier(A|B).  The values for uncertainty and sigma are included for 
    the user to assess/track the relative fidelity of the measurements but 
    are not used as part of the calculation performed.  The correction 
    performed is:

      corrected = victim - sum (source * coefficient)

    Currently for a correction to occur the absolute value of the coefficient
    must be greater than 1e-6 (note that negative crosstalk corrections are now
    possible).

  Saturation:
    Since the images operated on are unprocssed, DECam_crosstalk represents 
    the only point where the pixel data have had no changes applied (and 
    saturation can be unambiguously identified).  If the -satmask option 
    is given, the SATURATEA and SATURATEB keywords are used to flag pixels 
    with values greater than or equal to the saturation value.  These values 
    are recorded in a mask (as bitplane 2) and the number of saturate pixels 
    is recorded in the image header (NSATPIX).  Note, similar flagging can 
    occur elsewhere in the DESDM processing but this represents the only 
    point where the value are essentially guaranteed to be unaltered (as would 
    occur throughout many of the downstream processing/detrending steps).

  Overscan:
    If the -overscan options are provided then the overscan sections (specified
    by BIASECA and BIASECB) are used to evaluate and remove the instrumental 
    bias.  A number of options are possible which determine how the overscan 
    is evaluated to accomplish this step.  Once the overscan has been subtracted
    for each amplifier, the images are trimmed based on the TRIMSEC keyword 
    prior to output. The keyword DESOSCN is set to prevent subsequent processing
    from repeating this action.

    The additional arguments necessary to invoke overscan (besides the -overscan
    option itself) are:

    -overscan Causes a default version of the overscan correction keywords to
              be set.  Using this option is equivalent to setting:
         -overscansample 0
         -overscanfunction 0
         -overscanorder  1
         -overscantrim   0

    -overscansample <N>
      Specifies how the overscan is initially characterized.  Here the choice 
      of N detemines how the pixel values in each overscan line are combined 
      (prior to further processing).  The options are:
         N=-1 --> will detemine the median value in each overscan line,
         N=0  --> will calculate the average value, and
         N=1  --> will calculate the average after rejecting the highest and 
                  lowest value in each line.

    -overscantrim <ncols>, default=0
      Used to specify the number of start and end columns in the BIASEC regions 
      to ignore. For example -overscantrim 3 would ignore the first and last 
      three pixels in each line when evaluating the overscan.

    -overscanfunction <M>
      Specifies the method to evaluate the overscan (in the columnar direction).
      Choices are determined by the value of N (abs(N) is also used as a 
      binning/sampling factor). The behavior of the routine is as follows.  
      For:
        M<0 --> a cubic spline interpolation (after binning by abs(M) lines),
        M=0 --> causes the overscan to be evaluated on a line-by-line basis 
                (i.e. no fitting),
        M>0 --> will fit the overscan with a set of basis functions (Legendre 
                polynomials), after binning by abs(N) lines).

    -overscanorder <order>
      Is used to specify the order of the fit, order, that will be used for the 
      option -overscanfunction with M>0.  This controls the maximum order of 
      the basis functions (Legendre polynomials) that will be used in the fit.  
      Currently this is restricted to order<7 as higher order would potentially 
      be better characterized by the cubic spline option.

  Photometricity Flag:
    Currently the DESDM framework depends on a photometricity flag to be set 
    before fits to the atmospheric extinction are attempted (near the end of 
    single epoch processing).  In the DESDM pipelines this is accomplished by
    setting the flag at the outset.  The current mode of operation is to 
    assume photometricity until it is shown not to be the case. DECam_crosstalk
    provides the first point where this flag can be set (and subsequently 
    ingested with the image metadata.  The option -photflag is used to set the 
    initial value of the photometricity flag.

Known "Features":
  - Similar to many DESDM codes, the output path for files is checked and 
    created if necessary.  If an error occurs a wait/retry scheme is used 
    because our codes have been known to overwhelm many HPC file systems.  
    These waits can significantly increase the time if a filepath must be 
    checked/built.
  - The current handling of the DATASEC (BIASSEC/POSTSEC/TRIMSEC) have been
    been revamped.  Currently the LTV1 and LTV2 are identically set to 0
    and the CRPIX terms are not altered.  When SISPI begins populating 
    headers with a WCS solution that accruately accounts for the raw image
    format, this change needs to be added.
  - During late testing, a bias jump became apparent in CCDs that shared a 
    backplane with the focus chips (which are smaller and finish reading 
    earlier than the science CCDs populating the image plane.  Currently, 
    the DESDM framework is operating with:
      -overscanfunction 10 -overscanorder 0 
    i.e. fitting a constant value to the overscan.  This assumes that the 
    described bias jumps are stable and can be removed using the biascor... 
    an alternative is to use line-by-line assessment of the overscan but 
    this is known to create a modest correlated line-wise noise in the images.

*/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include "imsupport.h"
#include "argutils.h"

#define DOUBLE 1
#define STRING 2
#define FLOAT  3
#define INTEG  4

#define INCLUDECOMMENTS 0
#define DONTINCLUDECOMMENTS 1

#define NUM_KEYS 200
#define XTDIM 148

/* there are 70 CCDs in currect configuration: 
   1-62, 63-64, (skipping 65-68), 69-74.
   CCDs above 61 are from focus chips */
#define CCDNUM 74
#define CCDNUM1 75
#define CCDNUM2 76

/* there should be 71 headers in the image file
   header 1: applies to all images
   headers 2-71: ccd-specific */
#define HDUNUM 71

/* Amp order */
#define AmpAB 1
#define AmpBA -1

// Prototypes for external functions 
//char *strip_path(const char *path);
int getlistitems(char *strlist, char *indxlist, int maxitems);
static const char *svn_id = "$Id$";
static int flag_sat = 0;
static int flag_overscan = 0;
static int flag_focus = 0;

void print_usage(char *program_name)
{
  printf("%s <infile.fits> <outfile> <options>\n",program_name);
  printf("  -crosstalk <crosstalk matrix- file>\n");
  printf("  -photflag <0 or 1>\n");
  printf("  -satmask\n");
  printf("  -overscan\n");
  printf("  -overscansample <-1 for MEDIAN, 0 for MEAN, 1 for MEAN w/MINMAX>, \n");
  printf("  -overscanfunction < -50 < N < -1 for cubic SPLINE, 0 for LINE_BY_LINE, 1 < N < 50 for legendre polynomial>\n");
  printf("  -overscanorder <1-6>, default=1 <order for legendre polynomial>\n");
  printf("  -overscantrim <ncols>, default=0 <trim ncols at both edges of overscan, default = 0>\n");
  printf("  -maxhdunum <integer>, default=62 <all valid images>\n");
  printf("  -ccdlist <comma-separated list of CCDs to process>\n");
  printf("  -hdulist <comma-separated list of HDUs to process>\n");
  printf("  -focuschipsout\n");
  printf("  -verbose <0-3>\n");
  printf("  -help (print help message and exit)\n");
  printf("  -version (print version and exit)\n");
}

int DECamXTalk(int argc,char *argv[])
{
  /*  ********************************************************  */
  /*  ************  List the Headers Keywords NOT to be ******  */
  /*  ************      copied from 0th Header          ******  */
  /*  ********************************************************  */
		
  char	*zeroheader,**exclkey, *ptr,
    exclist[40][10]={
    "SIMPLE","BITPIX","NAXIS","NAXIS1",
    "NAXIS2","EXTEND","NEXTEND","CHECKSUM",
    "DATASUM","CHECKVER","GCOUNT","PCOUNT",
    "BZERO","BSCALE","INHERIT","XTENSION",
    "TFIELDS","TFORM1","ZIMAGE","ZTILE1",
    "ZTILE2","ZCMPTYPE","ZNAME1","ZVAL1",
    "ZNAME2","ZVAL2","EXTNAME","ZTENSION",
    "ZBITPIX","ZNAXIS","ZNAXIS1","ZNAXIS2",
    "ZPCOUNT","ZGCOUNT","DCREATED","TTYPE1",
    "ZHECKSUM"},
    *nthheader;
    int	numzerokeys,nexc=37,numnthkeys,numwrittenkeys=0,
      keyflag,keynum,totkeynum;
        int     writeCounter, writeImageFailed, retryDelay, nRetry;
	
    /*  ******************************************************** */
    /*  ******** array for storing crosstalk coefficients ****** */
    /*  ******************************************************** */
    float	xtalk[XTDIM][XTDIM];
	
    /*  ******************************************************** */
    /*  ******* other variables needed for the processing ****** */
    /*  ******************************************************** */
    char	filename[500],outname_temp[500],trash[200],nfilename[500],tag1[100],tag2[100],
      newname[500],filetype[50],obstype[80],comment[100],oldcomment[100],command_line[1000],
      longcomment[10000],event[10000],newimagename[500],*xtalk_filename=NULL,
      flag_hdus_or_ccds = 0, ccdlist[CCDNUM2], hdulist[CCDNUM2], AmpOrder[CCDNUM1];
    int	anynull,nfound,i,hdunum,hdutype,j,k,x,y,chdu, ccdnum,locout,
      flag_crosstalk=0,flag_verbose=1,flag_phot=0,bitpix,naxis,
      ext1,ext2,ampextension,ampA,ampB,locA,locB,ampoffset1,nkeys,
      ampoffset2,nsata=0,nsatb=0,flag_osorder=0,
      maxhdunum=HDUNUM, nhdus, needpath=1,
      mkpath(),getheader_flt(),getheader_str();
    static int status=0;
    float	value,uncertainty,significance,nullval,*outdata=NULL,*outdata2=NULL,**indata;
    float       ampAmin = 1e+30,ampAmax = 0,satval;
    float       ampBmin = 1e+30,ampBmax = 0;
    float       new_ltv1,new_ltv2,old_ltv1,old_ltv2;
    double      new_crpix1,new_crpix2,old_crpix1,old_crpix2;
    long	axes[2],naxes[2],taxes[2],pixels,npixels,fpixel,oldpixel=0;
    double	xtalkcorrection=0.0;
    void	printerror(),reportevt(),init_desimage(),decodesection();
    fitsfile *fptr,*nfptr;
    desimage des_image;
    FILE	*inp;
    time_t	tm;
    desimage input_image;
    desimage output_image;
    overscan_config osconfig;
    short *maskdata = NULL, *maskdata2 = NULL;
    char	delkeys[100][10]={"PRESECA","PRESECB","POSTSECA",
				  "POSTSECB","TRIMSECA","TRIMSECB","TRIMSEC",
				  "BIASSECA","BIASSECB",""};

    enum {OPT_CROSSTALK=1,OPT_PHOTFLAG,OPT_SATMASK,OPT_OVERSCAN,OPT_OVERSCANSAMPLE,OPT_OVERSCANFUNCTION,
	  OPT_OVERSCANORDER,OPT_OVERSCANTRIM,OPT_MAXHDUNUM,OPT_CCDLIST,OPT_HDULIST,OPT_FOCUSCHIPSOUT,
	  OPT_VERBOSE,OPT_HELP,OPT_VERSION};
    
    if(build_command_line(argc,argv,command_line,1000) <= 0){
      reportevt(2,STATUS,1,"Failed to record full command line.");
    }

    retryDelay = 150;
    nRetry = 20;

    /* RAG: Added to print version of code to standard output (for logs) */
    sprintf(event,"%s",svn_id);
    reportevt(2,STATUS,1,event);
    reportevt(2,STATUS,1,command_line);

    if (argc<2) {
      print_usage(argv[0]);
      exit(0);
    }

    //    sprintf(filename,"%s",argv[1]);
    /*sprintf(nfilename,"%s",argv[2]);*/

    osconfig.sample   = 0;
    osconfig.function = 0;
    osconfig.order    = 1;
    osconfig.trim     = 0;
    osconfig.debug    = 0;

    for (i = 0; i < CCDNUM2; i++) 
        ccdlist[i] = hdulist[i] = 0;

    /* ****************************************************** */
    /* ****************  Process Command Line *************** */
    /* ****************************************************** */

    // Still need to scan command line for verbosity so we can 
    // properly print status messages during command line parse.
    for (i=1;i<argc;i++) 
      {
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
    flag_sat = 0;
    flag_overscan = 0;
    while(1){
      int curind = optind;
      static struct option xtalk_options[] =
	{
	  {"crosstalk",        required_argument, 0,        OPT_CROSSTALK},
	  {"photflag",         required_argument, 0,         OPT_PHOTFLAG},
	  {"overscansample",   required_argument, 0,   OPT_OVERSCANSAMPLE},
	  {"overscanfunction", required_argument, 0, OPT_OVERSCANFUNCTION},
	  {"overscanorder",    required_argument, 0,    OPT_OVERSCANORDER},
	  {"overscantrim",     required_argument, 0,     OPT_OVERSCANTRIM},
	  {"maxhdunum",        required_argument, 0,        OPT_MAXHDUNUM},
	  {"ccdlist",          required_argument, 0,          OPT_CCDLIST},
	  {"hdulist",          required_argument, 0,          OPT_HDULIST},
	  {"verbose",          required_argument, 0,          OPT_VERBOSE},
	  {"version",          no_argument,       0,          OPT_VERSION},
	  {"help",             no_argument,       0,             OPT_HELP},
	  {"satmask",          no_argument,       &flag_sat,            1},
	  {"overscan",         no_argument,       &flag_overscan,       1},
	  {"focuschipsout",    no_argument,       &flag_focus,          1},
	  {0,0,0,0}
	};
      int clopx = 0;
      clop = getopt_long_only(argc,argv,"",
			      xtalk_options,&clopx);
      if(clop == -1)
	break;
      switch(clop){
      case 0:
	// For straightforward flags
	if(xtalk_options[clopx].flag != 0)
	  break;
	printf("Option %s is set",xtalk_options[clopx].name);
	if(optarg)
	  printf(" with %s",optarg);
	printf(".\n");
	break;
      case OPT_CROSSTALK: // -crosstalk
	cloperr = 0;
	flag_crosstalk=1;
	for (j=0;j<XTDIM;j++) 
	  for (k=0;k<XTDIM;k++) 
	    xtalk[j][k]=0.0;
	if(optarg){
	  xtalk_filename = strip_path(optarg);
	  /* read through crosstalk matrix file */
	  if (flag_verbose>=3) {
	    sprintf(event,"Reading file %s\n",xtalk_filename);
	    reportevt(flag_verbose,STATUS,1,event);
	  }
	  inp=fopen(optarg,"r");
	  if (inp==NULL) {
	    sprintf(event,"Crosstalk file %s not found",xtalk_filename);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	  //	strip_path(xtalk_filename);
	  while (fgets(trash,200,inp)!=NULL) {
	    if (!strncmp(trash,"ccd",3) || !strncmp(trash,"CCD",3)) {
	      /* first replace parentheses and commas with spaces */
	      for (j=0;j<strlen(trash);j++) 
		if (!strncmp(&(trash[j]),"(",1) || 
		    !strncmp(&(trash[j]),")",1)) trash[j]=32; 
	      sscanf(trash,"%s %s %f %f %f",tag1,tag2,
		     &value,&uncertainty,&significance);
	      if (!strncmp(tag1+strlen(tag1)-1,"A",1)) ampoffset1=0;
	      else if (!strncmp(tag1+strlen(tag1)-1,"B",1)) ampoffset1=1;
	      else {
		sprintf(event,"Amp not properly noted as A/B in %s",tag1);
		reportevt(flag_verbose,STATUS,5,event);
		exit(1);
	      }
	      if (!strncmp(tag2+strlen(tag2)-1,"A",1)) ampoffset2=0;
	      else if (!strncmp(tag2+strlen(tag2)-1,"B",1)) ampoffset2=1;
	      else {
		sprintf(event,"Amp not properly noted as A/B in %s",tag2);
		reportevt(flag_verbose,STATUS,5,event);
		exit(1);
	      }
	      /* strip off the A/B suffix */	
	      tag1[strlen(tag1)-1]=0;tag2[strlen(tag2)-1]=0; 
	      /* read in the ccd number */       
	      sscanf(&(tag1[3]),"%d",&ext1);
	      sscanf(&(tag2[3]),"%d",&ext2);
	      ext1=(ext1-1)*2+ampoffset1;
	      if (ext1<0 || ext1>=XTDIM) {
		sprintf(event,"CCD number out of range: %s %d",tag1,ext1);
		reportevt(flag_verbose,STATUS,5,event);
		exit(1);
	      }
	      ext2=(ext2-1)*2+ampoffset2;
	      if (ext2<0 || ext2>=XTDIM) {
		sprintf(event,"CCD number out of range: %s %d",tag2,ext2);
		reportevt(flag_verbose,STATUS,5,event);
		exit(0);
	      }
	      xtalk[ext1][ext2]=value;
	    }
	  }
	  if (fclose(inp)) {
	    sprintf(event,"File close failed: %s",xtalk_filename);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	  if (flag_verbose>3) /* print the crosstalk coefficients */
	    for (j=0;j<XTDIM;j++) {
	      if (j%2==0) printf("  ccd%02dA",j/2+1);
	      else printf("  ccd%02dB",j/2+1);
	      for (k=0;k<XTDIM;k++) {
		printf("  %7.5f",xtalk[j][k]);
	      }
	      printf("\n");
	    }
	}
	else{
	  cloperr = 1;
	  reportevt(flag_verbose,STATUS,5,"Option -crosstalk requires an argument specifying the crosstalk file.");
	  exit(1);
	}
	break;
      case OPT_PHOTFLAG: // -photflag
	cloperr = 0;
	if(optarg)
	  sscanf(optarg,"%d",&flag_phot);
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -photflag requires an argument.");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_OVERSCANSAMPLE: // -overscansample
	cloperr = 0;
 	flag_overscan = 1;
	if(optarg){
	  sscanf(optarg,"%d",&osconfig.sample);
	  if (osconfig.sample<-1 || osconfig.sample>1) {
	    sprintf(event,"Invalid Overscan sample, must be -1, 0 or 1.");
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -overscansample requires an argument.");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_OVERSCANFUNCTION: // -overscanfunction
	cloperr = 0;
 	flag_overscan = 1;
	if(optarg){
	  sscanf(optarg,"%d",&osconfig.function);
	  if (osconfig.function<-51 || osconfig.function>51) {
	    sprintf(event,"Invlaid Overscan function, must be <-50 to -1>, 0, or <1 to 50>.");
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -overscanfunction requires an argument.");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_OVERSCANORDER: // -overscanorder
	cloperr = 0;
 	flag_overscan = 1;
	flag_osorder = 1;
	if(optarg){
	  sscanf(optarg,"%d",&osconfig.order);
	  if (osconfig.order<1 || osconfig.order>6) {
	    sprintf(event,"Invalid Overscan order,  must be <1-6>.");
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	} 
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -overscanorder requires an argument.");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_OVERSCANTRIM: // -overscantrim
	cloperr = 0;
 	flag_overscan = 1;
	if(optarg){
	  sscanf(optarg,"%d",&osconfig.trim);
	  if (osconfig.trim > 20) {
	    sprintf(event,"Invalid Overscan trim, must be < 20 cols.");
	    reportevt(2,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -overscantrim requires an argument.");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_MAXHDUNUM: // -maxhdunum
	cloperr = 0;
	if(optarg){
	  sscanf(optarg,"%d",&maxhdunum);
	  if (maxhdunum<1 || maxhdunum>HDUNUM) {
	    sprintf(event,"Invlaid maxhdunum, must be <1 to %i>.", HDUNUM);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -maxhdunum requires an argument.");
	  exit(1);
	}
	break;
      case OPT_CCDLIST: // -ccdlist
	cloperr = 0;
	flag_hdus_or_ccds = 1;
	if(optarg){
	  if (getlistitems(optarg, ccdlist, CCDNUM2) < 1) {
	    sprintf(event,"Invalid ccdlist.");
	    reportevt(2,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -ccdlist requires an argument.");
	  command_line_errors++;
	  exit(1);
	}
	break;
      case OPT_HDULIST: // -hdulist
	cloperr = 0;
	flag_hdus_or_ccds = 1;
	if(optarg){
	  if (getlistitems(optarg, hdulist, CCDNUM2) < 1) {
	    sprintf(event,"Invalid hdulist.");
	    reportevt(2,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -hdulist requires an argument.");
	  command_line_errors++;
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
      sprintf(filename,"%s",argv[optind]);
    }
    else {
      reportevt(flag_verbose,STATUS,5,"Missing required input FITS image.");
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



    /* process all images (hdus & ccds) unless -ccdlist or -hdulist */
    if (!flag_hdus_or_ccds)
        for (i = 0; i < CCDNUM2; i++) 
            ccdlist[i] = hdulist[i] = 1;
    /*
      for (i = 0; i < CCDNUM2; i++) 
      printf("ccdlist[%i] = %i; hdulist[%i] = %i\n", i, ccdlist[i], i, hdulist[i]);
    */
    // Do some interparameter validation
    if(flag_osorder){ // if order option was given
      if(osconfig.function < 1){
	sprintf(event,"Option ignored: -overscanorder ignored if -overscanfunction < 1.");
	reportevt(flag_verbose,STATUS,3,event);
      }
    } 
    else if (osconfig.function > 0){
      reportevt(flag_verbose,STATUS,3,"Overscan Legendre Polynomial defaulting to order 1.");
    }

    /**********************************************************/
    /***********************  Open Image **********************/
    /**********************************************************/
	
    /* open the file, checked for all supported flavors */
    /* Flavor=Y (uncompressed) */
    if (fits_open_file(&fptr,filename,READONLY,&status)) { 
      status=0;
      sprintf(newname,"%s.fz",filename); 
      //      sprintf(event,"File open failed: %s",filename);
      //      reportevt(flag_verbose,STATUS,3,event);
      /* Flavor=F (fpacked) */
      if (fits_open_file(&fptr,newname,READONLY,&status)) {
	//	sprintf(event,"File open failed: %s",newname);
	//	reportevt(flag_verbose,STATUS,3,event);
	status=0;
	sprintf(newname,"%s.gz",filename); 
	/* Flavor=G (gzipped) */
	if (fits_open_file(&fptr,newname,READONLY,&status)) {
	  sprintf(event,"File open failed: %s{.fz,.gz}",filename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
    }
    if (flag_verbose) {
      sprintf(event,"Opened exposure: %s",newname);	
      reportevt(flag_verbose,STATUS,1,event);
    }
		
    /* get the number of HDUs in the image */
    if (fits_get_num_hdus(fptr,&hdunum,&status)) {
      sprintf(event,"Reading HDUNUM failed: %s",filename);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (flag_verbose) {
      sprintf(event,"%s has %d HDUs",filename,hdunum);
      reportevt(flag_verbose,STATUS,1,event);
    }
    if (hdunum!=HDUNUM) {
      sprintf(event,"%d HDUs not standard DECam format",hdunum);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }

    /* allocate memory for indata */
    indata=(float **)calloc(CCDNUM1, sizeof(float *));
    if (!indata) {
      sprintf(event,"Could not calloc indata");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }


    /**********************************************************/
    /*************  Read Header from 0th Extension ************/
    /**********************************************************/
	
    exclkey=(char **)calloc(nexc,sizeof(char *));
    if (!exclkey) {
      sprintf(event,"Could not calloc exclkey");
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    for (j=0;j<nexc;j++) {
      exclkey[j]=(char *)calloc(10,sizeof(char));
      if (!exclkey[j]) {
	sprintf(event,"Could not calloc exclkey[%d]",j);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      sprintf(exclkey[j],"%s",exclist[j]);
    }
    if (fits_hdr2str(fptr,INCLUDECOMMENTS,exclkey,nexc,
		     &zeroheader,&numzerokeys, &status)) {
      sprintf(event,"Could not read kewords");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* cut off last bogus "END" line */
    zeroheader[strlen(zeroheader)-80]=0;
    if (flag_verbose) {
      printf("  Read %d header keywords from 0th header\n",numzerokeys);
      if (flag_verbose > 4) {
	printf("  *************************************************\n");
	printf("  %s",zeroheader);
	printf("*************************************************\n");
      }
    }

    /* read the OBSTYPE and use it to determine <filetype> */
    if (fits_read_key_str(fptr,"OBSTYPE",obstype,comment,&status)) {
      sprintf(event,"Header Keyword OBSTYPE not found");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    else { /* map OBSTYPE into <filetype> */
      if (!strcmp(obstype,"object") || !strcmp(obstype,"OBJECT"))
	sprintf(filetype,"raw_obj");
      else if (!strcmp(obstype,"bias") || !strcmp(obstype,"BIAS")
	       || !strcmp(obstype,"zero") || !strcmp(obstype,"ZERO"))
	sprintf(filetype,"raw_bias");
      else if (!strcmp(obstype,"sky flat") || !strcmp(obstype,"SKY FLAT"))
	sprintf(filetype,"raw_tflat");
      else if (!strcmp(obstype,"dark") || !strcmp(obstype,"DARK"))
	sprintf(filetype,"raw_dark");
      else if (!strcmp(obstype,"dome flat") || !strcmp(obstype,"DOME FLAT")
	       || !strcmp(obstype,"flat") || !strcmp(obstype,"FLAT"))
	sprintf(filetype,"raw_dflat");
      else {
	sprintf(event,"OBSTYPE=%s not recognized",obstype);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
    }

    /**********************************************************/
    /**********  Read All Image Data into Memory **************/
    /**********    Enables Crosstalk Correction  **************/
    /**********************************************************/
    nhdus = 0;  /* nhdus should always be <= maxhdunum */

    if (flag_verbose) printf("  Reading image data from %s\n", filename);
    for (i=2;i<=hdunum;i++) {

      /* Move to the correct HDU */
      if (fits_movabs_hdu(fptr,i,&hdutype,&status)) {
	sprintf(event,"Move to HDU=%d failed: %s",i,filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	
      }

      /* get CCD number */
      ccdnum = i-1; /* for now */

      /* read current header */
      if (fits_hdr2str(fptr,INCLUDECOMMENTS,exclkey,nexc,&nthheader,
		       &numnthkeys,&status)) {
	sprintf(event,"Header %d not read",i-1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /* Get CCDNUM information from the header */
      if(getheader_str(nthheader,"CCDNUM",trash,flag_verbose)){
	sprintf(event,"CCDNUM not found in header. Using default sequence.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      else {
        sscanf(trash, "%d", &ccdnum);
      }

      /* figure out AmpOrder (AmpAB or AmpBA) */
      if(getheader_str(nthheader,"DATASECA",&input_image.dataseca,flag_verbose)){
	sprintf(event,"DATASECA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      decodesection(input_image.dataseca, input_image.datasecan, flag_verbose);

      if(getheader_str(nthheader,"DATASECB",&input_image.datasecb,flag_verbose)){
	sprintf(event,"DATASECB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      decodesection(input_image.datasecb, input_image.datasecbn, flag_verbose);
      AmpOrder[ccdnum] = (input_image.datasecan[0] < input_image.datasecbn[0]) ? AmpAB : AmpBA;

      free(nthheader); nthheader = NULL;

      /* error checking for ccdnum range */
      if (ccdnum < 1 || ccdnum > CCDNUM) {
	 sprintf(event,"Error: CCDNUM %d is out of range.", ccdnum);
	 reportevt(flag_verbose,STATUS,5,event);
	 printerror(status);
         exit(0);
      }

      if (ccdnum > 62 && !flag_focus) continue;

      /* stop splitting off HDUs after the first maxhdunum */
      if (++nhdus > maxhdunum) break;

      if (fits_get_img_param(fptr,2,&bitpix,&naxis,axes,&status)) {
	 sprintf(event,"Image params not found in %s",filename);
	 reportevt(flag_verbose,STATUS,5,event);
	 printerror(status);
      }

      /* compute image size, etc. */
      pixels  = axes[0] * axes[1];
      fpixel   = 1; nullval  = 0;

      /* double check image size */
      if (ccdnum < 63 ) { /* only for 2x4 chips */
          if (oldpixel) {
              if (pixels!=oldpixel) {
                  sprintf(event,"Image extensions have different sizes:  %ld  vs %ld",
                          pixels,oldpixel);
                  reportevt(flag_verbose,STATUS,5,event);
                  exit(0);
              }
          }
          else {
              oldpixel = pixels;
              naxes[0] = axes[0]; /* remember image size */
              naxes[1] = axes[1];
          }
      }

      indata[ccdnum]=(float *)calloc(pixels, sizeof(float));
      if (!indata[ccdnum]) {
	sprintf(event,"Error: Could not calloc indata[%d] %ld",ccdnum,pixels);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }

      sprintf(event,"HDU=%i, CCDNUM=%i", i, ccdnum);
      reportevt(flag_verbose,STATUS,1,event);

      /* read the CHDU image  */
      if (fits_read_img(fptr,TFLOAT,fpixel,pixels,&nullval,
			indata[ccdnum],&anynull,&status)) {
	sprintf(event,"Reading image extension %d failed: %s",i-1,
		filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /*if (flag_verbose) {printf(".");fflush(stdout);}*/
    }
    if (flag_verbose) { printf("\n"); fflush(stdout); }


    /**********************************************************/
    /**********     Cycle through extensions        ***********/
    /**********    preparing and writing data       ***********/
    /**********************************************************/

    /* naxes[0]=axes[0];naxes[1]=axes[1]; */
    npixels=naxes[0]*naxes[1];

    if (flag_sat) {
      maskdata = (short *)calloc(npixels,sizeof(short));
      if (maskdata == NULL) {
	sprintf(event,"Could not allocate maskdata(in)");
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }

      if (flag_overscan) {
        maskdata2 = (short *)calloc(npixels,sizeof(short));
        if (maskdata2 == NULL) {
	  sprintf(event,"Could not allocate maskdata(out)");
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
        }
      }
    }

    /* jump to 2nd HDU  
    if (fits_movabs_hdu(fptr,2,&hdutype,&status)) {
      sprintf(event,"Move to HDU=2 failed");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);	
    }
*/
    /* now cycle through the HDU's parsing properly and writing */
    /* useful FITS file */
    outdata = (float *)calloc(npixels,sizeof(float));
    if (outdata == NULL) {
	sprintf(event,"Error: Could not calloc outdata %ld",npixels);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
    }

    if (flag_overscan) {
      outdata2 = (float *)calloc(npixels,sizeof(float));
      if (outdata2 == NULL) {
        sprintf(event,"Error: Could not calloc outdata %ld",npixels);
        reportevt(flag_verbose,STATUS,5,event);
        exit(0);
      }
    }

    nhdus = 0;  /* nhdus should always be <= maxhdunum */

    for (i=2;i<=hdunum;i++) {

      /* Move to the correct HDU */
      if (fits_movabs_hdu(fptr,i,&hdutype,&status)) {
	sprintf(event,"Move to HDU=%d failed: %s",i,filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	
      }

      fits_get_hdu_num(fptr,&chdu);
      if (flag_verbose) printf("  Currently located at HDU %d of %d\n",
			       chdu,hdunum);
      if (chdu!=i) {
	sprintf(event,"Not located at correct HDU (%d instead of %d)",
		chdu,i);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }

      /* get CCD number */
      ccdnum = i-1; /* for now */

      /* read current header */
      if (fits_hdr2str(fptr,INCLUDECOMMENTS,exclkey,nexc,&nthheader,
		       &numnthkeys,&status)) {
	sprintf(event,"Header %d not read",i-1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* cut off last bogus "END" line */
      nthheader[strlen(nthheader)-80]=0;

      /* Get CCDNUM information from the header */
      if(getheader_str(nthheader,"CCDNUM",trash,flag_verbose)){
	sprintf(event,"CCDNUM not found in header. Using default sequence.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      else
        sscanf(trash, "%d", &ccdnum); 

     /* error checking for ccdnum range */
     if (ccdnum < 1 || ccdnum > CCDNUM) {
	sprintf(event,"Error: CCDNUM %d is out of range.", ccdnum);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
        exit(0);
     }

     /* does this image need to be processed? */
     if ((ccdnum > 62 && !flag_focus) || (ccdlist[ccdnum] != 1 && hdulist[i] != 1)) {
        if (nthheader != NULL) free(nthheader);
        continue;
     }

     /* stop splitting off HDUs after the first maxhdunum */
     if (++nhdus > maxhdunum) {
        if (nthheader != NULL) free(nthheader);
        break;
     }

      if (flag_verbose) {
	printf("  Read %d header keywords\n",numnthkeys);
	if (flag_verbose > 4) {
	  printf("  *************************************************\n");
	  printf("  %s",nthheader);
	  printf("*************************************************\n");
	}
      }

     if (fits_get_img_param(fptr,2, &bitpix, &naxis, naxes, &status)) {
	sprintf(event,"Image params not found in %s",filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
     }

     npixels  = naxes[0] * naxes[1];
     for(j = 0; j < npixels; j++)
	outdata[j] = indata[ccdnum][j];

      /* 
 	 Use DESIMAGE data structure to store information about
 	 the input image.  This sets input_image image buffers
 	 for all of the processing done in this utility. 
      */
      init_desimage(&input_image);
      input_image.image   = outdata;
      input_image.mask    = maskdata;
      input_image.npixels = npixels;
      input_image.axes[0] = naxes[0];
      input_image.axes[1] = naxes[1];

      init_desimage(&output_image);
      output_image.image  = outdata2;
      output_image.mask   = maskdata2;

      ampAmin = 1e+30;
      ampAmax = 0;
      ampBmin = 1e+30;
      ampBmax = 0;

      /**********************************************************/
      /*************  Read Header from Nth Extension ************/
      /**********************************************************/

      // Get a bunch of useful information from the header
      if(getheader_flt(nthheader,"SATURATA",&input_image.saturateA,flag_verbose)){
	sprintf(event,"SATURATA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }      
      if(getheader_flt(nthheader,"SATURATB",&input_image.saturateB,flag_verbose)){
	sprintf(event,"SATURATB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_flt(nthheader,"GAINA",&input_image.gainA,flag_verbose)){
	sprintf(event,"GAINA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }      
      if(getheader_flt(nthheader,"GAINB",&input_image.gainB,flag_verbose)){
	sprintf(event,"GAINB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_flt(nthheader,"RDNOISEA",&input_image.rdnoiseA,flag_verbose)){
	sprintf(event,"RDNOISEA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }      
      if(getheader_flt(nthheader,"RDNOISEB",&input_image.rdnoiseB,flag_verbose)){
	sprintf(event,"RDNOISEB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }      
      if(getheader_str(nthheader,"BIASSECA",&input_image.biasseca,flag_verbose)){
	sprintf(event,"BIASSECA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"BIASSECB",&input_image.biassecb,flag_verbose)){
	sprintf(event,"BIASSECB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"DATASECA",&input_image.dataseca,flag_verbose)){
	sprintf(event,"DATASECA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"DATASECB",&input_image.datasecb,flag_verbose)){
	sprintf(event,"DATASECB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"AMPSECA",&input_image.ampseca,flag_verbose)){
	sprintf(event,"AMPSECA not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"AMPSECB",&input_image.ampsecb,flag_verbose)){
	sprintf(event,"AMPSECB not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"TRIMSEC",&input_image.trimsec,flag_verbose)){
	sprintf(event,"TRIMSEC not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_str(nthheader,"DATASEC",&input_image.datasec,flag_verbose)){
	sprintf(event,"DATASEC not found in header.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      decodesection(input_image.biasseca,input_image.biassecan,flag_verbose);
      decodesection(input_image.biassecb,input_image.biassecbn,flag_verbose);
      decodesection(input_image.ampseca, input_image.ampsecan, flag_verbose);
      decodesection(input_image.ampsecb, input_image.ampsecbn, flag_verbose);
      decodesection(input_image.dataseca, input_image.datasecan, flag_verbose);
      decodesection(input_image.datasecb, input_image.datasecbn, flag_verbose);
      decodesection(input_image.trimsec, input_image.trimsecn, flag_verbose);
      decodesection(input_image.datasec, input_image.datasecn, flag_verbose);
      if (flag_verbose>=3) {
	sprintf(event,"BIASSECA=%s",input_image.biasseca);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"DATASECA=%s",input_image.dataseca);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"AMPSECA=%s",input_image.ampseca);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"BIASSECB=%s",input_image.biassecb);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"DATASECB=%s",input_image.datasecb);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"AMPSECB=%s",input_image.ampsecb);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"TRIMSEC=%s",input_image.trimsec);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"DATASEC=%s",input_image.datasec);
	reportevt(flag_verbose,STATUS,1,event);
      }
      /****************************************************************** */
      /* ************** apply the crosstalk correction ****************** */
      /* ************** and flag any saturated pixels   ****************** */
      /* **************************************************************** */

      if (ccdnum < 63) {

      /* Reset the mask if enabled */
      if(flag_sat){
	for(y = 0;y < npixels;y++)
	  maskdata[y] = 0;
      }
      nsata = 0;
      nsatb = 0;
      if (flag_crosstalk) {
        char thisAmpOrder = (input_image.datasecan[0] < input_image.datasecbn[0]) ? AmpAB : AmpBA;
	// Loop through *entire* input array.  This includes overscanned regions and 
	// any additional padding. 
	for (y=0;y<naxes[1];y++) 
	for (x=0;x<naxes[0];x++) {
	    locout=y*naxes[0]+x;
	    xtalkcorrection=0.0;
	    if(flag_sat){
	      // Determine whether the current position (wrt input array) is
	      // within the Amp-specific data region by checking DATASEC[A,B]
	      if(column_in_section(x+1,input_image.datasecan)){
		satval = input_image.saturateA;
		if(outdata[locout] > satval){
		  if(outdata[locout] > ampAmax) ampAmax = outdata[locout];
		  if(outdata[locout] < ampAmin) ampAmin = outdata[locout];
		  maskdata[locout] = BADPIX_SATURATE;
		  nsata++;
		}
	      }
	      else if(column_in_section(x+1,input_image.datasecbn)){
		satval = input_image.saturateB;
		if(outdata[locout] > satval){
		  if(outdata[locout] > ampBmax) ampBmax = outdata[locout];
		  if(outdata[locout] < ampBmin) ampBmin = outdata[locout];
		  maskdata[locout] = BADPIX_SATURATE;
		  nsatb++;
		}
	      }
	    }
            if (column_in_section(x+1,input_image.datasecan)) { /* in amp A */
              ampextension=(ccdnum-1)*2; 
              for (j=1;j<63;j++) {
                ampA=(j-1)*2;
                ampB=ampA+1; /* Amp B flipped because readout on far side */
                if (thisAmpOrder == AmpOrder[j]) {
                    locA=locout;
                    locB=y*naxes[0]+(naxes[0]-x-1);
                } else {
                    locB=locout;
                    locA=y*naxes[0]+(naxes[0]-x-1);
                }
                if (fabsf(xtalk[ampextension][ampA])>1.0e-6f) 
                  xtalkcorrection+=xtalk[ampextension][ampA]*indata[j][locA];
                if (fabsf(xtalk[ampextension][ampB])>1.0e-6f) 
                  xtalkcorrection+=xtalk[ampextension][ampB]*indata[j][locB];
              }
              outdata[locout]-=xtalkcorrection;
            }
            else if (column_in_section(x+1,input_image.datasecbn)) { /* in amp B */
              ampextension=(ccdnum-1)*2+1; 
              for (j=1;j<63;j++) {
                ampA=(j-1)*2;
                ampB=ampA+1; /* Amp B flipped because readout on far side */
                if (thisAmpOrder == AmpOrder[j]) {
                    locB=locout;
                    locA=y*naxes[0]+(naxes[0]-x-1);
                } else {
                    locA=locout;
                    locB=y*naxes[0]+(naxes[0]-x-1);
                }
                if (fabsf(xtalk[ampextension][ampA])>1.0e-6f)
                  xtalkcorrection+=xtalk[ampextension][ampA]*indata[j][locA];
                if (fabsf(xtalk[ampextension][ampB])>1.0e-6f) 
                  xtalkcorrection+=xtalk[ampextension][ampB]*indata[j][locB]; 
              }
              outdata[locout]-=xtalkcorrection;
            }
	  }
      }
      if (flag_verbose && flag_sat) {
	sprintf(event,"image=%s,CCD[%d] NSATPIX=%d (NSATA,NSATB)=(%d,%d)",filename,
		ccdnum,nsata+nsatb,nsata,nsatb);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"image=%s,CCD[%d] (SATLEVA,SATLEVB)=(%f,%f)",filename,
		ccdnum,input_image.saturateA,input_image.saturateB);
	reportevt(flag_verbose,STATUS,1,event);
	// Write out min and max of saturated pixels - trying to ferret out 
	// possible remaining A/B mixup.
	sprintf(event,"image=%s,CCD[%d] (SATAMIN,SATAMAX)=(%f,%f)",filename,
		ccdnum,ampAmin,ampAmax);
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"image=%s,CCD[%d] (SATBMIN,SATBMAX)=(%f,%f)",filename,
		ccdnum,ampBmin,ampBmax);
	reportevt(flag_verbose,STATUS,1,event);
      }
      } /* if (ccdnum < 63) */

       /* PERFORM OVERSCAN if not disabled */
       if(flag_overscan && ccdnum < 63) {
 	if(flag_verbose){
 	  sprintf(event,"OVERSCAN: (samp,func,ord,trim) = (%d,%d,%d,%d)",
 		  osconfig.sample,osconfig.function,
 		  osconfig.order,osconfig.trim);
 	  reportevt(flag_verbose,STATUS,1,event);
 	}
 	//  The "input" image to OverScan is the output image from the xtalk 
 	// Upon output , output_image is set up with the overscanned version
 	// of the input_image. 
 	OverScan(&input_image, &output_image, osconfig,flag_verbose);
 	
	output_image.datasecan[0] = (input_image.ampsecan[0] < input_image.ampsecan[1] ?
				    input_image.ampsecan[0] : input_image.ampsecan[1]);
	output_image.datasecan[1] = (input_image.ampsecan[0] < input_image.ampsecan[1] ?
				    input_image.ampsecan[1] : input_image.ampsecan[0]);
	output_image.datasecan[2] = (input_image.ampsecan[2] < input_image.ampsecan[3] ?
				    input_image.ampsecan[2] : input_image.ampsecan[3]);
	output_image.datasecan[3] = (input_image.ampsecan[2] < input_image.ampsecan[3] ?
				    input_image.ampsecan[3] : input_image.ampsecan[2]);
	output_image.datasecbn[0] = (input_image.ampsecbn[0] < input_image.ampsecbn[1] ?
				    input_image.ampsecbn[0] : input_image.ampsecbn[1]);
	output_image.datasecbn[1] = (input_image.ampsecbn[0] < input_image.ampsecbn[1] ?
				    input_image.ampsecbn[1] : input_image.ampsecbn[0]);
	output_image.datasecbn[2] = (input_image.ampsecbn[2] < input_image.ampsecbn[3] ?
				    input_image.ampsecbn[2] : input_image.ampsecbn[3]);
	output_image.datasecbn[3] = (input_image.ampsecbn[2] < input_image.ampsecbn[3] ?
				    input_image.ampsecbn[3] : input_image.ampsecbn[2]);
 	/* update datasec within the output image */
 	sprintf(output_image.datasec,"[%d:%ld,%d:%ld]",1,output_image.axes[0],
 		1,output_image.axes[1]);
 	sprintf(output_image.dataseca,"[%d:%d,%d:%d]",output_image.datasecan[0],
		output_image.datasecan[1],output_image.datasecan[2],
		output_image.datasecan[3]);
 	sprintf(output_image.datasecb,"[%d:%d,%d:%d]",output_image.datasecbn[0],
		output_image.datasecbn[1],output_image.datasecbn[2],
		output_image.datasecbn[3]);
 	if(flag_verbose){
 	  sprintf(event,"OVERSCAN: New DATASEC = %s",output_image.datasec);
 	  reportevt(flag_verbose,STATUS,1,event);	  
 	}
       }
       else { 
 	// If overscan was not used, then just straight copy the input image
 	output_image.image      = input_image.image;
 	output_image.mask       = input_image.mask;
 	output_image.axes[0]    = input_image.axes[0];
 	output_image.axes[1]    = input_image.axes[1];
	output_image.npixels    = input_image.npixels;  
 	strncpy(output_image.datasec,input_image.datasec,100);
	// 	output_image.gainA      = input_image.gainA;
	// 	output_image.gainB      = input_image.gainB;
	// 	output_image.rdnoiseA   = input_image.rdnoiseA;
	// 	output_image.rdnoiseB   = input_image.rdnoiseB;
	// 	output_image.saturateA  = input_image.saturateA;
	// 	output_image.saturateB  = input_image.saturateB;
       }
       //       strncpy(output_image.trimsec,input_image.trimsec,100);
       //       strncpy(output_image.biasseca,input_image.biasseca,100);
       //       strncpy(output_image.biassecb,input_image.biassecb,100);
       //       strncpy(output_image.ampseca,input_image.ampseca,100);
       //       strncpy(output_image.ampsecb,input_image.ampsecb,100);

      /* **************************************************************** */
      /* **************** Write New Single CCD Image ******************** */
      /* ************ Store Processing History in Header **************** */
      /* **************************************************************** */

      /* generate output filename */
      if (strstr(outname_temp, "%02d") != NULL)
      {
          /* new file name format in the form <pref>_%02d_<suf>[.fits] */
          sprintf(&nfilename[1], outname_temp, ccdnum);
      }
      else
      {
          /* old format in the form <pref>_ccdn[.fits] */
          if ((ptr = strstr(outname_temp, ".fits")) != NULL) 
          {
             /* .fits extension present */
             /* assume there is only one '.fits' substring present in the file name */ 
              *ptr = '\0'; /* remove .fits */
          }

          /* generate new file name */
          sprintf(&nfilename[1], "%s_%02d", outname_temp, ccdnum);
      }
      /* add .fits extension, if needed */
      if (strstr(&nfilename[1], ".fits") == NULL)
          strcat(&nfilename[1], ".fits");

      /* add pref '!' */
      nfilename[0] = '!';

      /* make sure path exists first time through */
      if (needpath) {
        needpath = 0;
	/*sprintf(nfilename,"%s_%02d.fits",argv[2],ccdnum);*/
	if (mkpath(nfilename+1,flag_verbose)) {
	  sprintf(event,"Failed to create path to file: %s",nfilename+1);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	else {
	  sprintf(event,"Created path to file: %s",nfilename);
	  reportevt(flag_verbose,STATUS,1,event);
	}
      }

      /*sprintf(nfilename,"!%s_%02d.fits",argv[2],ccdnum);*/
      if (fits_create_file(&output_image.fptr,nfilename,&status)) {
	sprintf(event,"File creation of %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) {
	sprintf(event,"Opened image for CCD %d : %s ",ccdnum,nfilename+1);
	reportevt(flag_verbose,STATUS,1,event);
      }
	  
      /* create image extension  */
      if (fits_create_img(output_image.fptr,FLOAT_IMG,2,output_image.axes,&status)) {
	sprintf(event,"Creating image failed: %s",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      //      /* ************************************************************* */
      //      /* ******************* Problem with resize_img ***************** */
      //      /* ************************************************************* */
      //      /* following approach causes problems with RICE tile compression*/
      //      /* copy the header from the CHDU in the input image */
      //      /*
      //	if (fits_copy_header(fptr,output_image.fptr,&status)) {
      //	sprintf(event,"Header copy into %s failed",nfilename);
      //	reportevt(flag_verbose,STATUS,5,event);
      //	printerror(status);
      //	}
      //	if (flag_verbose) 
      //	printf("  Copied header from %s[%i]\n",filename,i);
      //	/* resize the image */
      //      /*
      //	if (fits_resize_img(output_image.fptr,FLOAT_IMG,2,naxes,&status)) {
      //	sprintf(event,"Resizing image %s failed",nfilename);
      //	reportevt(flag_verbose,STATUS,5,event);
      //	printerror(status);             
      //	}
      //	/* end of problematic section */
      //      /* ************************************************************* */
      //      /* **************** END JUNK SECTION *************************** */
      //      /* ************************************************************* */

      /* first reserve room within the header for the known # of keywords */
      if (fits_set_hdrsize(output_image.fptr,numzerokeys+numnthkeys+15,&status)) {
	sprintf(event,"Reserving header space for %d keywords failed in %s",numzerokeys+numnthkeys+10,nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }


      /* copy zero header information into the new header */
      /* note that last keyword is null keyword and truncates the header! */
      /* only copy the header params that live within zerokeyword */
      if (fits_get_hdrpos(output_image.fptr,&totkeynum,&keynum,&status)) {
	sprintf(event,"Reading header position in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      numwrittenkeys=0;
      for (j=0;j<numzerokeys-1;j++) {
	/* don't write keyword if it is duplicate */
	keyflag=0;
	for (k=0;k<numnthkeys;k++) {
	  if (!strncmp(zeroheader+j*80,nthheader+k*80,8)) {
	    keyflag=1;
	    break;
	  }
	}
	/* copy the record into the header */
	if (!keyflag) {
	  if (fits_insert_record(output_image.fptr,totkeynum+numwrittenkeys,
				 zeroheader+j*80,&status)) {
	    sprintf(event,"Zero Header insert in %s failed",nfilename);
	    reportevt(flag_verbose,STATUS,5,event);
	    printerror(status);
	  }
	  numwrittenkeys++;
	}
      }
      if (flag_verbose) printf("%d keywords from 0th header, ",
			       numwrittenkeys-3);

      /* copy the header information from the current HDU of MEF */
      if (fits_get_hdrpos(output_image.fptr,&totkeynum,&keynum,&status)) {
	sprintf(event,"Reading header position in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose>=3) {
	sprintf(event,"Currently at header position %d with %d keywords",
		keynum,totkeynum);
	reportevt(flag_verbose,STATUS,1,event);
      }
      numwrittenkeys=0;
      for (j=0;j<numnthkeys-1;j++) {
	/* copy the record into the header */
	if (fits_insert_record(output_image.fptr,totkeynum+numwrittenkeys,
			       nthheader+j*80,&status)) {
	  sprintf(event,"Insert of %dth header keyword in %s failed:\n%s",
	          numwrittenkeys,nfilename,nthheader+j*80);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	numwrittenkeys++;
      }
      if (flag_verbose) 
	printf("  Copied %d keywords from current header, ",
	       numwrittenkeys-3);
	    
      
      /* remove unneeded information from the header */
      if (flag_overscan && ccdnum < 63) {
	nkeys=0;
	while (strlen(delkeys[nkeys])) {
	  if (fits_read_keyword(output_image.fptr,delkeys[nkeys],
				comment,comment,&status)==KEY_NO_EXIST) status=0;
	  else {
	    if (fits_delete_key(output_image.fptr,delkeys[nkeys],&status)) {
	      if (flag_verbose) {
		sprintf(event,"Keyword %s not deleted from image %s",
			delkeys[nkeys],output_image.name+1);
		reportevt(flag_verbose,STATUS,3,event);
	      }
	      status=0;
	    }
	  }
	  nkeys++;
	}
      }
      /* write image */
      
      writeCounter = 0;
      while (1) {
	
	writeImageFailed = fits_write_img(output_image.fptr,TFLOAT,1,output_image.npixels,output_image.image,&status);
	
	if(writeImageFailed) {
	  
	  writeCounter++;
	  sprintf(event,"Writing image  %s failed on try %i", nfilename, writeCounter);
	  reportevt(flag_verbose,STATUS,4,event);
	  
	  sleep(retryDelay);
	  
	  if(writeCounter < nRetry) { continue; }
	  
	}
	/* If we get here, we wrote the file,
	   unless writeCounter is up to nRetry
	   break in either case */
	if(writeCounter == nRetry) {
	  sprintf(event,"Writing image %s failed on final try %i", nfilename, writeCounter);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	break;
	
      }
      
      if (flag_verbose) printf("  Wrote image data, ");

      /* want only the image name-- search for first occurence of "/" */
      for (j=strlen(nfilename);j>=0;j--) 
	if (!strncmp(&(nfilename[j]),"/",1)) break;
      sprintf(newimagename,"%s",(nfilename+j+1));
      for (j=strlen(newimagename);j>=0;j--) 
	if (!strncmp(&(newimagename[j]),".fits",5)) break;
      newimagename[j]=0;
      if (!strncmp(newimagename,"!",1)) 
	sprintf(newimagename,"%s",newimagename+1);
      if (fits_modify_key_str(output_image.fptr,"FILENAME",newimagename,"",&status)) {
	sprintf(event,"Keyword FILENAME modify in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) printf("FILENAME = %s, ",nfilename+1);

      if (fits_insert_record(output_image.fptr,totkeynum+(++numwrittenkeys),
			     "COMMENT ------------------------",
			     &status)) {
	sprintf(event,"Comment insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /* add photflag */
      if (fits_update_key_lng(output_image.fptr,"PHOTFLAG",flag_phot,
			      "Night Photometric (1) or not (0)",&status)) {
	sprintf(event,"Keyword PHOTFLAG insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) printf("PHOTFLAG = %d\n",flag_phot);

      /* update OBSTYPE keyword */
      if (fits_update_key_str(output_image.fptr,"OBSTYPE",filetype,
			      "Type of observation",&status)) {
	sprintf(event,"Keyword OBSTYPE update in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) printf("  OBSTYPE = %s\n",filetype);

      /* leave processing history in the header */
      tm=time(NULL);
      sprintf(comment,"%s",asctime(localtime(&tm)));
      comment[strlen(comment)-1]=0;
      if (fits_write_key_str(output_image.fptr,"DESDCXTK",comment,
			     "DECam image conversion and crosstalk correction",
			     &status)) {
	sprintf(event,"Keyword DESDCXTK insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (ccdnum < 63) {
      if (flag_crosstalk) {
      if (fits_write_key_str(output_image.fptr,"XTALKFIL",xtalk_filename,
			     "Crosstalk file",
			     &status)) {
	sprintf(event,"Keyword XTALKFIL insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      }
      if(flag_overscan) {
	if (fits_write_key_str(output_image.fptr,"DESOSCN",comment,
			       "overscan corrected",&status)) {
	  sprintf(event,"Writing DESOSCN failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}	
	/* update the DATASEC keyword if overscan corrected */
	if (fits_update_key_str(output_image.fptr,"DATASEC",output_image.datasec,
				"Data section within image",&status)) {
	  sprintf(event,"Updating DATASEC failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_update_key_str(output_image.fptr,"DATASECA",output_image.dataseca,
				NULL,&status)) {
	  sprintf(event,"Updating DATASECA failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_update_key_str(output_image.fptr,"DATASECB",output_image.datasecb,
				NULL,&status)) {
	  sprintf(event,"Updating DATASECB failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}

/*         printf(" RAG trim_sec: %d:%d,%d:%d \n",input_image.trimsecn[0],input_image.trimsecn[1],input_image.trimsecn[2],input_image.trimsecn[3]);  */

        /* The following unilaterally set the LTV1 and LTV2 keywords to 0.0 */
        /* RAG 09/26/2012: An optional set of code exists below to calculate changes based on trimsec 
           but require that LTV1 and LTV2 are populated correctly and that correctly would be to have positive 
           offsets which is still under contention  */

	if (fits_update_key_flt(output_image.fptr,"LTV1",0.0,2,NULL,&status)) {
	  sprintf(event,"Update LTV1 failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	if (fits_update_key_flt(output_image.fptr,"LTV2",0.0,2,NULL,&status)) {
	  sprintf(event,"Update LTV2 failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}

/*      These are commented out for the time being pending the results of the discussion 
        on how to properly set the LTV keywords in the raw image headers */

#if 0
        if (fits_read_key_flt(output_image.fptr,"LTV1",&old_ltv1,oldcomment,&status)==KEY_NO_EXIST){
          status=0;
          new_ltv1=0.0;
          if (fits_write_key(output_image.fptr,TFLOAT,"LTV1",&new_ltv1,NULL,&status)){
            sprintf(event,"Adding LTV1 failed: %s",output_image.name+1);
            reportevt(flag_verbose,STATUS,5,event);
            printerror(status);
          } 
        }else{
          new_ltv1=old_ltv1-(input_image.trimsecn[0]-1.0);
          if (fits_update_key(output_image.fptr,TFLOAT,"LTV1",&new_ltv1,NULL,&status)){
            sprintf(event,"Update LTV1 failed: %s",output_image.name+1);
            reportevt(flag_verbose,STATUS,5,event);
            printerror(status);
          }
/*          printf(" RAG new LTV1: %f --> %f \n",old_ltv1,new_ltv1);  */
        }  

        if (fits_read_key_flt(output_image.fptr,"LTV2",&old_ltv2,oldcomment,&status)==KEY_NO_EXIST){
          status=0;
          new_ltv2=0.0;
          if (fits_write_key(output_image.fptr,TFLOAT,"LTV2",&new_ltv2,NULL,&status)){
            sprintf(event,"Adding LTV2 failed: %s",output_image.name+1);
            reportevt(flag_verbose,STATUS,5,event);
            printerror(status);
          } 
        }else{
          new_ltv2=old_ltv2-(input_image.trimsecn[2]-1.0);
          if (fits_update_key(output_image.fptr,TFLOAT,"LTV2",&new_ltv2,NULL,&status)){
            sprintf(event,"Update LTV2 failed: %s",output_image.name+1);
            reportevt(flag_verbose,STATUS,5,event);
            printerror(status);
          }
/*          printf(" RAG new LTV2: %f --> %f \n",old_ltv2,new_ltv2); */
        }
#endif

/*      RAG: 09/26/2012 If the raw WCS has been set to properly account for the BIAS and POSTSCAN image 
        sections then the follow code will be needed to update CRPIX1 and CRPIX2 */


        if (fits_read_key_dbl(output_image.fptr,"CRPIX1",&old_crpix1,oldcomment,&status)==KEY_NO_EXIST){
          status=0;
          sprintf(event,"No CRPIX1 keyword detected in image: %s",output_image.name+1);
          reportevt(flag_verbose,STATUS,3,event);
        }else{ 
          new_crpix1=old_crpix1-(input_image.trimsecn[0]-1.0);
          if (fits_update_key(output_image.fptr,TDOUBLE,"CRPIX1",&new_crpix1,NULL,&status)){
            sprintf(event,"Update CRPIX1 failed: %s",output_image.name+1);
            reportevt(flag_verbose,STATUS,5,event);
            printerror(status);
          }else{
            sprintf(event,"Update CRPIX1 due to trim: %lf --> %lf ",old_crpix1,new_crpix1);
            reportevt(flag_verbose,STATUS,1,event);
          }         
/*          printf(" RAG new CRPIX1: %lf --> %lf \n",old_crpix1,new_crpix1);  */
        } 

        if (fits_read_key_dbl(output_image.fptr,"CRPIX2",&old_crpix2,oldcomment,&status)==KEY_NO_EXIST){
          status=0;
          sprintf(event,"No CRPIX2 keyword detected in image: %s",output_image.name+1);
          reportevt(flag_verbose,STATUS,3,event);
        }else{ 
          new_crpix2=old_crpix2-(input_image.trimsecn[2]-1.0);
          if (fits_update_key(output_image.fptr,TDOUBLE,"CRPIX2",&new_crpix2,NULL,&status)){
            sprintf(event,"Update CRPIX2 failed: %s",output_image.name+1);
            reportevt(flag_verbose,STATUS,5,event);
            printerror(status);
          }else{
            sprintf(event,"Update CRPIX2 due to trim: %lf --> %lf ",old_crpix2,new_crpix2);
            reportevt(flag_verbose,STATUS,1,event);
          }
/*         printf(" RAG new CRPIX2: %lf --> %lf \n",old_crpix2,new_crpix2);  */
        } 

      }      
      if(flag_sat) { 
	if (fits_write_key_str(output_image.fptr,"DESSAT",comment,
			       "saturated pixels flagged",&status)) {
	  sprintf(event,"Writing DESSAT failed: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}	
	if (fits_write_key_lng(output_image.fptr,"NSATPIX",(long)(nsata+nsatb),
			       "Number of saturated pixels",&status)) {
	  sprintf(event,"Setting number of saturated pixels failed: %ld",(long)(nsata+nsatb));
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      } /* if (ccdnum < 63) */
      sprintf(longcomment,"DESDM:");
      //      for (j=0;j<argc;j++) sprintf(longcomment,"%s %s",longcomment,argv[j]);
      sprintf(longcomment,"%s %s",longcomment,command_line);
      if (fits_write_history(output_image.fptr,longcomment,&status)) {
	sprintf(event,"DESDM: longcomment insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose>=3) printf("  => %s\n",longcomment);
      /* mark this as an image extension */
      if (fits_write_key_str(output_image.fptr,"DES_EXT","IMAGE","Image extension",
			     &status)) {
	sprintf(event,"Keyword DES_EXT insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      if(flag_sat && ccdnum < 63){
	/* Create extension for the mask */
	if (fits_create_img(output_image.fptr,USHORT_IMG,2,output_image.axes,&status)) {
	  sprintf(event,"Creating image mask failed.");
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* Write the mask data */	  
	if (fits_write_img(output_image.fptr,TUSHORT,1,output_image.npixels,output_image.mask,
			   &status)) {
	  sprintf(event,"Writing image mask failed.");
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* Indicate it's a MASK extension */
	if (fits_update_key_str(output_image.fptr,"DES_EXT","MASK",
				"Extension type",&status)) {
	  reportevt(flag_verbose,STATUS,5,"Setting DES_EXT=MASK failed.");
	  printerror(status);
	}
	/* May be needed */
	if (fits_update_key_lng(output_image.fptr,"EXTVER",2,
				"Extension version",&status)) {
	  reportevt(flag_verbose,STATUS,5,"Setting EXTVER=2 failed");
	  printerror(status); 
	}
      }

      /* close the new image */
      if (fits_close_file(output_image.fptr,&status)) {
	sprintf(event,"Closing %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) {
	sprintf(event,"Closed image %ld X %ld  OBSTYPE = %s PHOTFLAG = %d : %s",
		output_image.axes[0],output_image.axes[1],filetype,flag_phot,nfilename+1);
	reportevt(flag_verbose,STATUS,1,event);
      }


      /*if (flag_verbose) printf("  Closed image %s\n",nfilename+1);*/
      /* Move to the next HDU 
      if (i<hdunum) 
	if (fits_movrel_hdu(fptr,1,&hdutype,&status)) {
	  sprintf(event,"Move to next HDU failed");
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      */

      if (nthheader != NULL) free(nthheader);
    }


    /* ********************************************************* */
    /* ************** Close Input Image and Exit *************** */
    /* ********************************************************* */

    if (fits_close_file(fptr,&status)) {
      sprintf(event,"Closing Exposure failed");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (flag_verbose) {
      sprintf(event,"Closed image %s",filename);
      reportevt(flag_verbose,STATUS,1,event);
    }

    /* free memory */
    if (zeroheader != NULL) free(zeroheader);

    if (indata != NULL) {
       for (i = 0; i < CCDNUM1; i++)
           if (indata[i] != NULL) free(indata[i]);
       free(indata);
    }

    if (exclkey != NULL) {
       for (i = 0; i < nexc; i++)
           if (exclkey[i] != NULL) free(exclkey[i]);
       free(exclkey);
    }

    if (maskdata != NULL) free(maskdata);
    if (outdata != NULL) free(outdata);

    if (outdata2 != NULL) free(outdata2);
    if (maskdata2 != NULL) free(maskdata2);

    return(0);
}

int grabkeywordval(header,numkeys,keyword,keytyp,keyval)
     char header[],keyword[];
     int numkeys,keytyp;
     void *keyval;
{
  int	i,loc,j=0,found=0;
  char 	*val_str;
  float	*val_flt;
  double	*val_dbl;
  int	*val_int;
	
  for (i=0;i<numkeys;i++) {
    loc=i*80;
    /* find the keyword */
    if (!strncmp(keyword,header+loc,strlen(keyword))) {
      loc+=9;
      /* strip of any leading spaces or ' */
      while (!strncmp(header+loc," ",1) || 
	     !strncmp(header+loc,"'",1)) loc++;
      found++;
      if (keytyp==STRING) {
	val_str=(char *)keyval;
	/* copy over until "/" */
	while (strncmp(header+loc,"/",1)) val_str[j++]=header[loc++];
	j--;
	/* strip off any trailing spaces of ' */
	while (!strncmp(val_str+j," ",1) || 
	       !strncmp(val_str+j,"'",1)) j--;
	val_str[j+1]=0;
	/*printf("  KW: %8s %s\n",keyword,val_str);*/
	break;
      }
      if (keytyp==FLOAT)  {
	val_flt=(float *)keyval;
	sscanf(header+loc,"%f",val_flt);
	/*printf("  KW: %8s %f\n",keyword,*val_flt);*/
	break;
      }
      if (keytyp==DOUBLE) {
	val_dbl=(double *)keyval;
	sscanf(header+loc,"%lf",val_dbl);
	/*printf("  KW: %8s %f\n",keyword,*val_dbl);*/
	break;
      }
      if (keytyp==INTEG) {
	val_int=(int *)keyval;
	sscanf(header+loc,"%d",val_int);
	/*printf("  KW: %8s %d\n",keyword,*val_int);*/
	break;
      }
    }
  }
  return(found);
}

int getheader_flt(hdr,key,value,flag_verbose)
     char hdr[],key[];
     float *value;
     int flag_verbose;
{
  int	len,lines,flag=0,keylen,i;
  char	event[500] ,nline[100];
  void	reportevt();
  len=strlen(hdr);
  keylen=strlen(key);
  lines=len/80;
  for (i=0;i<lines;i++) {
    if (!strncmp(key,hdr+i*80,keylen)) {
      flag=1;
      len = strlen(hdr+i*80);
      strncpy(nline,hdr+i*80,80);
      nline[80] = '\0';
      strtok(nline,"=");
      *value = atof(strtok(NULL,"="));
      break;
    }
  }

  if (!flag){
    if(flag_verbose) {
      sprintf(event,"Keyword %s not found.",key);
      reportevt(flag_verbose,STATUS,5,event); 
      /*printf("******\n%s\n*******\n",hdr); */
      exit(0);  
    }
    else
      return(1);
  }
  //  else if (flag_verbose>=3) {
  //    sprintf(event,"Keyword %s found to be %.2f",key,*value);
  //    reportevt(flag_verbose,STATUS,1,event);
  //  }
  return(0);
}

 
int getheader_str(char *hdr,char *key,char *value,int flag_verbose)
{
  int	len,lines,flag=0,keylen,i;
  char	event[500],nline[100];
  void	reportevt();
  char *answer;
  len=strlen(hdr);
  keylen=strlen(key);
  lines=len/80;
  for (i=0;i<lines;i++) {
    if (!strncmp(key,hdr+i*80,keylen)) {
      flag=1;
      len = strlen(hdr+i*80);
      strncpy(nline,hdr+i*80,80);
      nline[80] = '\0';
      strtok(nline,"=");
      answer = strtok(NULL,"=");
      len = strlen(answer);
      strncpy(value,answer,len);
      value[len] = '\0';
      break;
    }
  }

  if (!flag){
    if(flag_verbose) {
      sprintf(event,"Keyword %s not found.",key);
      reportevt(flag_verbose,STATUS,5,event); 
      /*printf("******\n%s\n*******\n",hdr); */
      exit(0);  
    }
    else
      return(1);
  }
  //else if (flag_verbose>=3) {
  //    sprintf(event,"Keyword %s found to be %s",key,value);
  //    reportevt(flag_verbose,STATUS,1,event);
  //}
  return(0);
}

/* expects strlist in comma-separated format, e.g., 1,2,5,65,3
   sets corresponding indxlist entries to 1
   ignores everything >= maxitems
   returns number of set values, or -1 if any list is NULL */
int getlistitems(char *strlist, char *indxlist, int maxitems)
{
    char *ptr = strlist;
    int indx, count = 0;

    if (strlist == NULL || indxlist == NULL) return -1;

    /* clear list */
    for (indx = 0; indx < maxitems; indx++)
        indxlist[indx] = 0;

    while (*ptr != '\0')
    {
        if (sscanf(ptr, "%d", &indx) == 1) {
            if (indx > 0 && indx < maxitems) {
                indxlist[indx] = 1;
                count++;
            }
        }

        while (*ptr != ',' && *ptr != '\0')
            ptr++;

        if (*ptr == ',') ptr++;
    }

    return count;
}


int main(int argc,char *argv[])
{
  return(DECamXTalk(argc,argv));
}

#undef STRING
#undef FLOAT
#undef DOUBLE
#undef INTEG
