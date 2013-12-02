/* mkillumcor
 *
 * Takes an input image (likely produced using mksupersky)
 * and produces a high signal to noise (i.e. smoothed) illumination
 * correction and also, with appropriate flags, a fringe correction 
 * image.
 *
 * 06/27/2012: modified by V. Kindratenko
 *  - Added code to ignore hot/dead (BADPIX_BPM) pixels
 *  - Introduced -ignoredeadpixels to ignore hot/dead pixels
 *
 * 05/24/2013: modified by V. Kindratenko*
 *  - Added -fast option
 */

#include <time.h>
#include <string.h>
#include <getopt.h>

#include "imsupport.h"
#include "argutils.h"

#define  MAXPIXELS 500
static const char *svn_id = "$Id$";

void print_usage(char *program_name)
{
  printf("%s <input file> <options>\n",program_name);
  printf("  Smoothing Options\n");
  printf("    -minsize <pixels>\n");
  printf("    -maxsize <pixels>\n");
  printf("    -median (default)\n");
  printf("    -average\n");
  printf("    -ranseed <#>\n");
  printf("    -ignoredeadpixels\n");
  printf("  Output Options\n");
  printf("    -output_illum <file>\n");
  printf("    -output_fringe <file>\n");
  printf("  Testing Options\n");
  printf("    -scaleregion <Xmin:Xmax,Ymin,Ymax>\n");
  printf("    -illumcor_compare <template>\n");
  printf("    -fringecor_compare <template>\n");
  printf("    -verbose <0-3>\n");
  printf("    -help    (print usage and exit)\n");
  printf("    -version (print version and exit)\n");
  printf("  Performance Options\n");
  printf("    -fast\n");
}

static int flag_median  = YES;
static int flag_aver    = NO;
static int flag_deadpix = NO;
static int flag_fast    = NO;

int main(int argc,char *argv[])
{
  char	filter[100]="",comment[10000],longcomment[10000],event[10000],
    scaleregion[100],inname_temp[1000],command_line[1000];
  int	i,j,x,y,k,l,flag_fringe=NO,flag_illum=YES,
    minsize=4,maxsize=10,loc,dx,dy,xmin,xmax,ymin,ymax,
    count,localcount,flag_verbose=1,
    nvec,xlen,ylen,totpix,mkpath(),ccdnum=0,ncompare=0,
    flag_illumcor_compare=NO,flag_fringecor_compare=NO,
    scaleregionn[4]={500,1500,1500,2500},scalenum;
  static int status=0;
  long	ranseed=0;
  float	*vecsort=NULL,ran1(),subsample_fraction,*randnum=NULL,offset,rms,maxdev,
    mode,fwhm,scalefactor,*scalesort;
  double	value;
  desimage input,illum,fringe,illumcor_template,fringecor_template;	
  void	rd_desimage(),shell(),reportevt(),image_compare(); 
  time_t	tm;
  enum { OPT_MINSIZE =1 ,OPT_MAXSIZE,OPT_MEDIAN,OPT_AVERAGE,OPT_RANSEED,OPT_IGNOREDEADPIXELS,
	 OPT_OUTPUT_ILLUM,OPT_OUTPUT_FRINGE,OPT_SCALEREGION,OPT_ILLUMCOR_COMPARE,
	 OPT_FRINGECOR_COMPARE,OPT_VERBOSE,OPT_HELP,OPT_VERSION};
  

  if (argc<2) {
    print_usage(argv[0]);
    exit(0);
  }

  command_line[0] = '\0';
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
      if (flag_verbose<0 || flag_verbose>3) {
	sprintf(event,"Verbose level out of range %d . Reset to 2",
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
    static struct option illumcor_options[] =
      {
	{"minsize",          required_argument, 0,            OPT_MINSIZE},
	{"maxsize",          required_argument, 0,            OPT_MAXSIZE},
	{"ranseed",          required_argument, 0,            OPT_RANSEED},
	{"illumcor_compare", required_argument, 0,            OPT_ILLUMCOR_COMPARE},
	{"fringecor_compare",required_argument, 0,            OPT_FRINGECOR_COMPARE},
        {"output_illum",     required_argument, 0,            OPT_OUTPUT_ILLUM},
        {"output_fringe",    required_argument, 0,            OPT_OUTPUT_FRINGE},
	{"verbose",          required_argument, 0,            OPT_VERBOSE},
	{"scaleregion",      required_argument, 0,            OPT_SCALEREGION},
	{"version",          no_argument,       0,            OPT_VERSION},
	{"help",             no_argument,       0,            OPT_HELP},
	{"average",          no_argument,       &flag_aver,   YES},
	{"median",           no_argument,       &flag_median, YES},
	{"ignoredeadpixels", no_argument,       &flag_deadpix,YES},
	{"fast",             no_argument,       &flag_fast,   YES},
	{0,0,0,0}
      };


    int clopx = 0;

    clop = getopt_long_only(argc,argv,"",illumcor_options,&clopx);

    if(clop == -1)
      break;
    switch(clop){
    case 0:
      // For straightforward flags
      if(illumcor_options[clopx].flag != 0)
	break;
      printf("Option %s is set",illumcor_options[clopx].name);
      if(optarg)
	printf(" with %s",optarg);
      printf(".\n");
      break;
    case OPT_MINSIZE: // -minsize
      sscanf(optarg,"%d",&minsize);
      break;
    case OPT_MAXSIZE: // -maxsize
      sscanf(optarg,"%d",&maxsize);
      break;
    case OPT_RANSEED: // -ranseed
      sscanf(optarg,"%ld",&ranseed);
      break;
    case OPT_ILLUMCOR_COMPARE: // -illumcor_compare
      cloperr = 0;
      flag_illumcor_compare=YES;
      sprintf(illumcor_template.name,"%s",optarg);
      if (!strncmp(&(illumcor_template.name[strlen(illumcor_template.name)-5]),".fits",5)  
	  && !strncmp(&(illumcor_template.name[strlen(illumcor_template.name)-8]),".fits.gz",8))  {
	sprintf(event,"Template image must be FITS: %s",illumcor_template.name);
	reportevt(flag_verbose,STATUS,5,event);
	exit(1);
      }
      break;
    case OPT_FRINGECOR_COMPARE: // -fringecor_compare
      flag_fringecor_compare=YES;
      sprintf(fringecor_template.name,"%s",optarg);
      if (!strncmp(&(fringecor_template.name[strlen(fringecor_template.name)-5]),".fits",5)  
	  && !strncmp(&(fringecor_template.name[strlen(fringecor_template.name)-8]),".fits.gz",8))  {
	sprintf(event,"Template image must be FITS: %s",fringecor_template.name);
	reportevt(flag_verbose,STATUS,5,event);
	exit(1);
      }
      break;
    case OPT_OUTPUT_ILLUM: // -output_illum
      sprintf(illum.name,"!%s",optarg);
      flag_illum=YES;
      break;
    case OPT_OUTPUT_FRINGE: // -output_fringe
      sprintf(fringe.name,"!%s",optarg);
      flag_fringe=YES;
      break;      
    case OPT_SCALEREGION: // -scaleregion
      sprintf(scaleregion,"%s",optarg);
      decodesection(scaleregion,scaleregionn,flag_verbose);
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
    
  if(flag_median && flag_aver){
    reportevt(flag_verbose,STATUS,5,"Error: Cannot choose both MEDIAN and AVERAGE.");
    exit(1); 
  }
    
  /* ********************************** */
  /* ********************************** */
  /* *******  READ INPUT IMAGE  ******* */
  /* ********************************** */
  /* ********************************** */

  sprintf(input.name,"%s",inname_temp);
  rd_desimage(&input,READONLY,flag_verbose);

  /* check whether input image produced by DESMKSKY */
  headercheck(&input,filter,&ccdnum,"DESMKSKY",flag_verbose);
  sprintf(event,"Input image FILTER=%s and CCDNUM=%d: %s",filter,ccdnum,
	  input.name);
  reportevt(flag_verbose,STATUS,1,event);

  /* *********************************************************** */
  /* *********************************************************** */
  /* **************  SET UP TO CREATE ILLUMINATION IMAGE ******* */
  /* *********************************************************** */
  /* *********************************************************** */

  reportevt(flag_verbose,STATUS,1,"Creating illumination correction");
  illum.npixels=input.npixels;
  illum.nfound=input.nfound;
  for (i=0;i<illum.nfound;i++) illum.axes[i]=input.axes[i];
  illum.bitpix=FLOAT_IMG;
  illum.image=(float *)calloc(illum.npixels,sizeof(float));
  if (illum.image==NULL) {
    reportevt(flag_verbose,STATUS,5,"Calloc of illum.image failed");
    exit(0);
  }
  illum.mask=(short *)calloc(illum.npixels,sizeof(short));
  if (illum.mask==NULL) {
    reportevt(flag_verbose,STATUS,5,"Calloc of illum.mask failed");
    exit(0);
  }
  for (i=0;i<illum.npixels;i++) illum.mask[i]=input.mask[i];

  /* ************************************************************ */
  /* *********** set up for median sorting if needed ************ */
  /* ************************************************************ */
  if (flag_median) {
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
  /* smooth the input image to create the illumination correction */
  /* ************************************************************ */
  for (y=0;y<illum.axes[1];y++) {
    if (y%200==1) {printf(".");fflush(stdout);}
    dy=maxsize;
    if (y-dy<0)  dy=y;
    if (y+dy>=illum.axes[1])  dy=illum.axes[1]-y-1;
    if (dy<minsize) dy=minsize;
    ymin=y-dy;if (ymin<0) ymin=0;
    ymax=y+dy;if (ymax>illum.axes[1]) ymax=illum.axes[1];
    for (x=0;x<illum.axes[0];x++) {
      dx=maxsize;
      if (x-dx<0)  dx=x;
      if (x+dx>=illum.axes[0])  dx=illum.axes[0]-x-1;
      if (dx<minsize) dx=minsize;
      xmin=x-dx;if (xmin<0) xmin=0;
      xmax=x+dx;if (xmax>illum.axes[0]) xmax=illum.axes[0];
      loc=x+y*illum.axes[0];
      illum.image[loc]=0.0f;
      /* beginning of processing loop */
      if (flag_aver) {
	count=0;value=0.0;
	for (k=ymin;k<ymax;k++) for (l=xmin;l<xmax;l++) {
	  if ((flag_deadpix==YES) && (input.mask[l+k*input.axes[0]]&BADPIX_BPM)) continue;
	  value+=input.image[l+k*input.axes[0]];
	  count++;
	}
	if (count > 0)
	  illum.image[loc]=value/(float)count;
      }
      else if (flag_median) {
	count=0;
	ylen=ymax-ymin; 
	xlen=xmax-xmin; 
	totpix=ylen*xlen;
	/* use all the pixels */
	if (!ranseed || ylen*xlen<MAXPIXELS) 
	  for (k=ymin;k<ymax;k++) for (l=xmin;l<xmax;l++) {
	    if ((flag_deadpix==YES) && (input.mask[l+k*input.axes[0]]&BADPIX_BPM)) continue;
	    vecsort[count++]=input.image[l+k*input.axes[0]];
	  } 
	else { /* randomly choose the pixels */
	  localcount = 0;
	  while (count<MAXPIXELS && localcount<MAXPIXELS) {
	    k=ymin+(int)(totpix*randnum[localcount])/xlen;
	    if (k>=ymax) k=ymax-1;
	    l=xmin+(int)(totpix*randnum[localcount])%xlen;
	    if (l>=xmax) l=xmax-1;
	    localcount++;
	    if ((flag_deadpix==YES) && (input.mask[l+k*input.axes[0]]&BADPIX_BPM)) continue;
	    vecsort[count++]=input.image[l+k*input.axes[0]]; 
	  }
	}
	if (count > 0) {
          if (flag_fast)
            illum.image[loc] = quick_select(vecsort, count);
          else
          {
	    /* sort */
	    shell(count,vecsort-1);
	    /* odd or even number of pixels */
	    if (count%2) illum.image[loc]=vecsort[count/2];
	    else illum.image[loc]=0.5*(vecsort[count/2]+vecsort[count/2-1]);
	  }
	}
      }
      /* end of processing loop */
    }
  }
  printf("\n");

  /* set up to retrieve scale of images */
  scalenum=(scaleregionn[1]-scaleregionn[0]+1)*(scaleregionn[3]-
						scaleregionn[2]+1);
  scalesort=(float *)calloc(scalenum,sizeof(float));
  if (scalesort==NULL) {
    sprintf(event,"Calloc for scalesort failed");
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }

  if (flag_fast)
    retrievescale_fast(&illum,scaleregionn,scalesort,flag_verbose,
		&scalefactor,&mode,&fwhm);
  else
    retrievescale(&illum,scaleregionn,scalesort,flag_verbose,
		&scalefactor,&mode,&fwhm);

  /* ************************************************************ */
  /* ************************************************************ */
  /* ************** Test image against template ***************** */
  /* ************************************************************ */
  /* ************************************************************ */
  if (flag_illumcor_compare) {
    /*  Read template image */
    rd_desimage(&illumcor_template,READONLY,flag_verbose);
    /* check image for ccdnumber */
    headercheck(&illumcor_template,filter,&ccdnum,"DESMKICR",
		flag_verbose);
    /* first compare images */
    rms=maxdev=offset=0.0;
    ncompare=0;
    image_compare(&illum,&illumcor_template,&offset,&rms,&maxdev,
		  &ncompare,flag_verbose);
    /* issue STATUS events according to differences measured */

    /* ***************************************************** */
    /* ***************************************************** */
    /* ***************************************************** */
  }


  /* **************************************************** */
  /* *********  SET UP TO CREATE  FRINGE IMAGE ********** */
  /* **************************************************** */

  if (flag_fringe) {
    reportevt(flag_verbose,STATUS,1,"Creating fringe correction");
    fringe.npixels=input.npixels;
    fringe.nfound=input.nfound;
    for (i=0;i<fringe.nfound;i++) fringe.axes[i]=input.axes[i];
    fringe.bitpix=FLOAT_IMG;
    fringe.image=(float *)calloc(fringe.npixels,sizeof(float));
    if (fringe.image==NULL) {
      reportevt(flag_verbose,STATUS,5,"Calloc of fringe.image failed");
      exit(0);
    }
    fringe.mask=(short *)calloc(fringe.npixels,sizeof(short));
    if (fringe.mask==NULL) {
      reportevt(flag_verbose,STATUS,5,"Calloc of fringe.mask failed");
      exit(0);
    }
    for (i=0;i<fringe.npixels;i++) fringe.mask[i]=input.mask[i];
    for (i=0;i<fringe.npixels;i++) 
      if (!input.mask[i]) fringe.image[i]=input.image[i]-illum.image[i];
      else fringe.image[i]=0.0f;
    printf("\n");
  }

  if (flag_fast)
    retrievescale_fast(&fringe,scaleregionn,scalesort,flag_verbose,
		&scalefactor,&mode,&fwhm);
  else
    retrievescale(&fringe,scaleregionn,scalesort,flag_verbose,
		&scalefactor,&mode,&fwhm);

  /* ************************************************************ */
  /* ************** Test image against template ***************** */
  /* ************************************************************ */
  if (flag_fringecor_compare) {
    /*  Read template image */
    rd_desimage(&fringecor_template,READONLY,flag_verbose);
    /* check image for ccdnumber */
    headercheck(&fringecor_template,filter,&ccdnum,"DESMKFRG",
		flag_verbose);
    /* first compare images */
    rms=maxdev=offset=0.0;
    ncompare=0;
    image_compare(&fringe,&fringecor_template,&offset,&rms,&maxdev,
		  &ncompare,flag_verbose);
    /* issue STATUS events according to differences measured */

    /* ***************************************************** */
    /* ***************************************************** */
    /* ***************************************************** */
  }




  /* ******************************************** */
  /* ********** WRITE OUTPUT IMAGE(S) *********** */
  /* ******************************************** */
  if (flag_illum) {
    sprintf(event,"Writing results to %s",illum.name+1);
    reportevt(flag_verbose,STATUS,1,event);

    /* make sure path exists for new image */
    if (mkpath(illum.name,flag_verbose)) {
      sprintf(event,"Failed to create path to file: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    else {
      sprintf(event,"Created path to file: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,1,event);
    }

    /* create the file */
    if (fits_create_file(&illum.fptr,illum.name,&status)) {
      sprintf(event,"Creating file failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }

    /* create image extension */
    if (fits_create_img(illum.fptr,FLOAT_IMG,2,illum.axes,&status)) {
      sprintf(event,"Creating image failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
	  
    /* write the corrected image*/
    if (fits_write_img(illum.fptr,TFLOAT,1,illum.npixels,illum.image,
		       &status)) {
      sprintf(event,"Writing image failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
	  
    /* write basic information into the header */
    if (fits_update_key_str(illum.fptr,"OBSTYPE","illumcor",
			    "Observation type",&status)) {
      sprintf(event,"Writing keyword OBSTYPE=illumcor failed: %s",
	      illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);  
    }
    if (fits_update_key_str(illum.fptr,"FILTER",filter,
			    "Filter name(s)",&status)) {
      sprintf(event,"Writing keyword FILTER=%s failed: %s",
	      filter,illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (fits_update_key_lng(illum.fptr,"CCDNUM",ccdnum,
			    "CCD number",&status)) {
      sprintf(event,"Writing keyword CCDNUM=%d failed: %s",
	      ccdnum,illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* Write information into the header describing the processing */
    /* get system time */
    tm=time(NULL);
    sprintf(comment,"%s",asctime(localtime(&tm)));
    comment[strlen(comment)-1]=0;
    if (fits_write_key_str(illum.fptr,"DESMKICR",comment,
			   "Created illumination correction",&status)) {
      sprintf(event,"Writing processing history failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    sprintf(longcomment,"DESDM:");
    //    for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);
    sprintf(longcomment,"%s %s",longcomment,command_line);
    reportevt(flag_verbose,STATUS,1,longcomment);
    if (fits_write_comment(illum.fptr,longcomment,&status)) {
      sprintf(event,"Writing longcomment failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (fits_write_key_str(illum.fptr,"DES_EXT","IMAGE",
			   "Image extension", &status)) {
      sprintf(event,"Writing DES_EXT=IMAGE failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* close the corrected image */
    if (fits_close_file(illum.fptr,&status)) {
      sprintf(event,"File close failed: %s",illum.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }


  if (flag_fringe) {
    sprintf(event,"Writing results to %s",fringe.name+1);
    reportevt(flag_verbose,STATUS,1,event);
    /* make sure path exists for new image */
    if (mkpath(fringe.name,flag_verbose)) {
      sprintf(event,"Failed to create path to file: %s",fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    else {
      sprintf(event,"Created path to file: %s",fringe.name+1);
      reportevt(flag_verbose,STATUS,1,event);
    }
    /* create the file */
    if (fits_create_file(&fringe.fptr,fringe.name,&status)) {
      sprintf(event,"Creating file failed : %s",fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* create image extension */
    if (fits_create_img(fringe.fptr,FLOAT_IMG,2,fringe.axes,&status)) {
      sprintf(event,"Creating image failed : %s",fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* write the corrected image*/
    if (fits_write_img(fringe.fptr,TFLOAT,1,fringe.npixels,fringe.image,
		       &status)) {
      sprintf(event,"Writing image failed : %s",fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* write basic information into the header */
    if (fits_write_key_str(fringe.fptr,"OBSTYPE","fringecor",
			   "Observation type",&status)) {
      sprintf(event,"Writing Keyword OBSTYPE=fringecor failed : %s",
	      fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);  
    }
    if (fits_update_key_str(fringe.fptr,"FILTER",filter,
			    "Filter name(s)",&status)) {
      sprintf(event,"Writing Keyword FILTER=%s failed : %s",
	      filter,fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (fits_update_key_lng(fringe.fptr,"CCDNUM",ccdnum,
			    "CCD number",&status)) {
      sprintf(event,"Writing Keyword CCDNUM=%d failed : %s",
	      ccdnum,fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }

    /* Write information into the header describing the processing */
    /* get system time */
    tm=time(NULL);
    sprintf(comment,"%s",asctime(localtime(&tm)));
    comment[strlen(comment)-1]=0;
    if (fits_write_key_str(fringe.fptr,"DESMKFRG",comment,
			   "Created fringe correction",&status)) {
      sprintf(event,"Writing processing history failed : %s",
	      fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    sprintf(longcomment,"DESDM:");
    for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);
    reportevt(flag_verbose,STATUS,1,longcomment);
    if (fits_write_comment(fringe.fptr,longcomment,&status)) {
      sprintf(event,"Writing longcomment failed : %s",
	      fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (fits_write_key_str(fringe.fptr,"DES_EXT","IMAGE",
			   "Image extension",&status)) {
      sprintf(event,"Writing keyword DES_EXT=IMAGE failed : %s",
	      fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* close the fringe image */
    if (fits_close_file(fringe.fptr,&status)) {
      sprintf(event,"Closing file failed : %s",
	      fringe.name+1);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }

  printf("\n");

  return(0);
}

#undef  MAXPIXELS 
