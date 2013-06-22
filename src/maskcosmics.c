/*--main
Basic syntax: maskcosmics <input image> <options>
  Input data:
    <input image> = Input FITS image to be used to mask cosmics rays.
  Options:
    -crays 
    -crfract <cosmic ray flux fraction> (default 0.20)
    -crsig2 <cosmic variance> (default 2)
    -verbose <0-3>
    -help    (print usage and exit)
    -version (print version and exit)

Summary:
  This routine takes an input MEF image (image, bpm, weight map) and identifies
  (defined by doCosmics routine) every pixel touched by a cosmics ray and 
  that is crsig2 sigmas above sky level. An additional 1 pixel buffer
  around all pixels identified is also added. All cosmics rays found with
  their repective buffer are masked in the bad pixel mask (BPM) plane of 
  the input image as well as in the weight map plane (a value of zero
  is given to the masked cosmics rays written in the weight map image plane).
  Routine search for the FWHM on the image header and adjusts crfract and crsig2
  in the following manner:

  if FWHM < 3.3 pixels, it reports STATUS 3 warning and does not try to find
  cosmics rays.
  if 3.3 >= FWHM < 4.0 pixels, it defines crfract = 0.15 and crsig2 = 1.5 for
  optimal cosmics rays detection in this FWHM range.
  if FWHM >= 4.0 pixels, it defines crfract = 0.2 and crsig2 = 2.0 for optimal
  cosmics rays detections in this FWHM range.
  If no FWHM present in image header, it defaults to crfract = 0.1 and 
  crsig2 =1.0, but this values are not optimal for cosmics rays detection.


Detailed Description:

    -crays
     maskcosmics code will be executed, Must be used or code will exit.

    -crfract <crfract>
     Cosmics rays flux fraction:
     If this option is used, it will be passed to the code, but it might get
     overwritten depending of the exitence of the FWHM in the image header, or the
     value of the FWHM (see summary for explanation).

    -crsig2 <crsig2>
     Sigma above sky level:
     If this option is used, it will be passed to the code, but it might get
     overwritten depending of the exitence of the FWHM in the image header, or the
     value of the FWHM (see summary for explanation).

     -verbose (0-3)
     Sets the verbosity level for reporting actions being taken as the processing progresses.                                                                                                                                                                                                                  
*/

#include "imsupport.h"
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "argutils.h"


#define NXSIZE 2048
#define NYSIZE 4096
#define NSIZE NXSIZE*NYSIZE

/* variables used in function dosky() */
#define NXSKY 4
#define NYSKY 8
#define NSKYSAMP 128

#define BADPIX_CGROW 128

//Parameter for Qsort
#define INSERTION_SORT_BOUND 16 /* boundary point to use insertion sort */
/* Utility Routines */


/* Global Variables */
float skybrite, skysigma, estgain;

static const char *svn_id = "$Id$";

void print_usage(char *program_name)
{
  printf("%s <input image> <options>\n",program_name);
  printf("  Masking Options\n");
  printf("    -crays \n");
  printf("    -crfract <cosmic ray flux fraction> (default 0.2)\n");
  printf("    -crsig2 <cosmic variance> (default 2)\n");
  printf("  Output Options\n");
  printf("    -verbose <0-3>\n");
  printf("    -help    (print usage and exit)\n");
  printf("    -version (print version and exit)\n");
 }

static int flag_horiz    = NO;

int stack[2048*4096];
int naxis1, naxis2,naxis;

/* Parameters for doComic routine. This routines determines which pixels are
   potential cosmics rays and mask them */


double LOW_FACT = 5.0;
int MIN_SAMPLE = 100;
int MIN_LINE = 3;
double WEIGHT_MIN = 1.e-6;
double SATURATED = 1.e99;
double WGT_FACT = 0.5;
int BORDER=15;
double STOL = 20.0;
int dobadCol=0;
int domaskCol=0;


 /* This structure is used by doSky() */
 typedef struct
 {
   int npixX;
   int npixY;
   double weightmin;
   int minsample;
   int minline;
   float saturated;
   float underflow;
   int nxbin, nybin;
   int itype;
   int ordx;
   int ordy;
   double skysol[4];
   double gainsol[4];
   double rms;
   double sigma;
 } skypar;

skypar sky;

int main(int argc,char *argv[])
{
  /*int MakeMask(int argc,char *argv[])
    {*/
   char	filter[100]="",comment[1000],longcomment[10000],event[10000],
     scaleregion[100],keycard[100],outputname[1000],astrostdsfile[10000],
     imagetype[1000],exposstr[1000],command_line[1000], updated_command_line[1000];
   char outputnamewithtempl[2500],inname_temp[1000];
   int	i,j,x,y,k,l,flag_fringe=NO,flag_illum=YES,xp,yp,th,horiztrails=0,
     minsize=4,maxsize=10,loc,dx,dy,xmin,xmax,ymin,ymax,scalenum,bitpix,
     count,flag_verbose=1,
     nvec,xlen,ylen,totpix,mkpath(),ccdnum=0,ncompare=0,
     flag_nofwhm=0,flag_illumcor_compare=NO,flag_fringecor_compare=NO,
     scaleregionn[4]={500,1500,1500,2500},ctr,
                        keysexist,numcrays,interp,xpospix,ypospix,
                        pixrow,pixcol,xval,yval;
   static int status=0;
   long	ranseed=0,seed = -30;
   float	*vecsort,ran1(),subsample_fraction,*randnum,offset,rms,maxdev,exposure,
     fwhm,scalefactor,*scalesort,mean=0,krad,growrad,gasdev(),std,sigma,regionmode,
     goodpixx[10],goodpixy[10],**pixval,outy,erry,badrow,badcolumn,xmask,ymask,
     skybrite,skysigma,saturate;
   double	value;
   desimage input,cray,output,tempimage;	
   void	rd_desimage(),shell(),reportevt(),image_compare(),
     retrievescale(),printerror();
   void    polin2(),creategrid(),headervalue(),pixelhisto(),getfloatheader();
   /* process definitions from makeWeight */
   int doCosmic(float *image,short *bpm,float *weight, float *tempweight, double crfract, double crsig2);
   int doStat(float *image,short *bpm,float *weight);
   int doWeight(float *image,short *bpm,float *weight);
   void dogrow(int masksrch,int npix,int maskset,float *image,short *bpm,float *weight);
   void spread(int masksrch,double sig,int maskset,float *image,short *bpm,float *weight);
   /* end of process definitions */

   time_t	tm;
   float	interp_noise=1.0,interp_fwhm=0.0;
   int	flag_output=0,flag_crays=0,flag_stars=0,flag_srcgrowrad=0, flag_crfract = 0, flag_crsig2 = 0,
     num = 0,filemode=0,goodpix,badpix,hdutype;
   enum {OPT_CRAYS=1,OPT_CRFRACT,OPT_CRSIG2,OPT_FLAG_HORIZTRAILS,
         OPT_NOINTERPOLATE,OPT_OUTPUT,OPT_VERBOSE,OPT_HELP,OPT_VERSION};

   
   float crfract,crsig2;

   int error;
   int mask;

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

   for (i=1;i<argc;i++){
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
     static struct option mask_options[] =
       {
     {"output",           required_argument, 0,          OPT_OUTPUT},
     {"verbose",          required_argument, 0,              OPT_VERBOSE},
     {"crfract",          required_argument, 0,              OPT_CRFRACT},
     {"crsig2",           required_argument, 0,              OPT_CRSIG2},
     {"crays",            no_argument,       0,              OPT_CRAYS},
     {"version",          no_argument,       0,              OPT_VERSION},
     {"help",             no_argument,       0,              OPT_HELP},
     {"flag_horiztrails", no_argument,       &flag_horiz,    NO},
     {0,0,0,0}
       };

     int clopx = 0;
     clop = getopt_long_only(argc,argv,"",
                 mask_options,&clopx);
     if(clop == -1)
       break;
     switch(clop){
     case 0:
       // For straightforward flags
       if(mask_options[clopx].flag != 0)
     break;
       printf("Option %s is set",mask_options[clopx].name);
       if(optarg)
         printf(" with %s",optarg);
       printf(".\n");
       break;
     case OPT_CRFRACT: // -crfract
         flag_crfract = YES;
         sscanf(optarg,"%f",&crfract);
         //if (crfract = '\0') {
         //reportevt(flag_verbose,STATUS,5,"Option -crfract requires an argument.");
         //exit(1);
         //}
       break;
     case OPT_CRSIG2: // -crsig2
         flag_crsig2 = YES;
         sscanf(optarg, "%f", &crsig2);
         //if (crsig2 = '\0') {
         //reportevt(flag_verbose,STATUS,5,"Option -crsig2 requires an argument.");
         //exit(1);
         //}
       break;
     case OPT_CRAYS: // -crays
       /*sprintf(cray.name,"%s",optarg);*/
       flag_crays=YES;
       break;
       //case OPT_STARS: // -stars
       //sprintf(astrostdsfile,"%s",optarg);
       //flag_stars=YES;
       //break;
       //case OPT_OUTPUT: // -output
       //sprintf(outputname,"%s",optarg);
       //flag_output=YES;
       //break;
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
     reportevt(flag_verbose,STATUS,5,"Missing required input FITS image.");
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
  
  
  
  /* ********************************** */
  /* *******  READ INPUT IMAGE  ******* */
  /* ********************************** */
  
  sprintf(input.name,"%s",inname_temp);
 
  if (flag_output) { 
    rd_desimage(&input,READONLY,flag_verbose);
   }
  else {
    rd_desimage(&input,READWRITE,flag_verbose);
  }


  fits_file_mode(input.fptr,&filemode,&status);

  headervalue(&input,"FILTER",filter,filemode);
    
  sprintf(event,"Input image FILTER=%s\n",filter,input.name);

  reportevt(flag_verbose,STATUS,1,event);
  
  headervalue(&input,"OBSTYPE",imagetype,filemode);
  getfloatheader(&input,"EXPTIME",&exposure,filemode);
  /*	printf ("OBSTYPE=%s\n",imagetype); */ 
  if(!strcmp(imagetype,"remap") || (flag_horiz ==1))horiztrails = 1 ;

  if (input.axes[0] != NXSIZE || input.axes[1] != NYSIZE)
    {
      sprintf(event,"Error: Wrong size of image NXSIZE=%i, NYSIZE=%i, image %s\n", input.axes[0], input.axes[1], input.name);
      reportevt(flag_verbose,STATUS,5,event);      
    } else 
    {
      naxis1 = input.axes[0];
      naxis2 = input.axes[1];
      naxis=2;
      printf("naxis=%li naxis1=%li naxis2=%li \n",
             naxis,naxis1,naxis2);
    }

  /**********************************************************/
  /*************  Read Header from Nth Extension ************/
  /**********************************************************/

  getfloatheader(&input,"SKYSIGMA",&skysigma,filemode);
  getfloatheader(&input,"SKYBRITE",&skybrite,filemode);
  getfloatheader(&input,"SATURATE",&saturate,filemode);
  getfloatheader(&input,"FWHM",&fwhm,filemode);
  flag_nofwhm=0;
  if (fits_read_key_flt(input.fptr,"FWHM",&fwhm,comment,&status) ==KEY_NO_EXIST) {
     flag_nofwhm=1;
     sprintf(event,"Keyword FWHM not defined: %s",input);
     reportevt(flag_verbose,STATUS,3,event);
  }

  if (flag_verbose>2) {
    sprintf(event,"Input image SKYSIGMA=%5.2f, SKYBRITE=%5.2f, SATURATE=%5.2f and FWHM=%5.2f: %s",skysigma,skybrite,saturate,fwhm,input.name);
    reportevt(flag_verbose,STATUS,1,event);
  }

  if (skybrite>0.0 && skysigma>0.0) estgain = skybrite/(skysigma*skysigma);
  else estgain = 0.0;
  sprintf(event,"Estimated gain=%f skybrite=%f skysigma=%f saturate=%f \n",estgain,skybrite,skysigma,saturate);
  reportevt(3, STATUS, 1, event);

  if (SATURATED>saturate)
    {
      printf("Saturation level requested is greater than level in fits header.  SATURATE=%f keyword value.\n",saturate);
      SATURATED = saturate;
    }
 
  /**********************************************************/
  /* Set crfract and crsig2 based on the seeing (FWHM) unless explicitly set on the command line. */
  /**********************************************************/

  /*    
    #if seeing <= 3.0pix No CR detection done to image
    #if if 3.0 < FWHM < 4.0 pix, crfract = 0.15, crsig2 = 2
    #If seeing >= 4.0 pix crfract = 0.2, crsig2 = 2
  */

  if (!flag_crfract){
    crfract=0.1;
    if (flag_nofwhm){
        sprintf(event,"No seeing measurement present.  Defaulting to CRFRACT=%.2f\n",crfract);
        reportevt(flag_verbose,STATUS,3,event);
    }
  }else{
    if (fwhm > 3.3 && fwhm < 4.0) {
      crfract = 0.15;
    }else if (fwhm >= 4.0) {
      crfract = 0.2;
    }else{
      sprintf(event,"Image with FWHM=%.2f (<3.3 pix) is not suitable for automatic detection of CR\n",fwhm);
      reportevt(flag_verbose,STATUS,3,event);
      exit(0);
    }
  }
  sprintf(event," Value for CRFRACT set automaically based on seeing to CRFACT=%.2f\n",crfract);
  reportevt(flag_verbose,STATUS,1,event);
  

  if (!flag_crsig2){
    crsig2=1.0;
    if (flag_nofwhm){
      sprintf(event,"No seeing measurement present.  Defaulting to CRSIG2=%.2f\n",crsig2);
      reportevt(flag_verbose,STATUS,3,event);
    }
  }else{
    if (fwhm > 3.3 && fwhm < 4.0){
      crsig2 = 1.5;
    }else if (fwhm >= 4.0) {
      crsig2 = 2.0;
    }else{
      sprintf(event,"Image with FWHM=%.2f (<3.3 pix) is not suitable for automatic detection of CR\n",fwhm);
      reportevt(flag_verbose,STATUS,3,event);
      exit(0);
    }
  }
  sprintf(event," Value for CRSIG2 set automaically based on seeing to CRSIG2=%.2f\n",crsig2);
  reportevt(flag_verbose,STATUS,1,event);


  sprintf(event,"Values used according to FWHM from image are: crfract=%.2f and crsig2=%.2f \n", crfract, crsig2);
  reportevt(flag_verbose,STATUS,1,event);

 
  /* create an output image with same parameters as input image */
  output=input;
  tempimage=input;
  
  /* copy the correct output name into place */
  if (!flag_output) { 
    sprintf(output.name,"%s",input.name);
  
  /* create image, mask and varim arrays and copy from input */
  output.image=(float *)calloc(output.npixels,sizeof(float));
  if (output.image==NULL) reportevt(flag_verbose,STATUS,5,
                                    "Calloc of output.image failed");
  output.mask=(short *)calloc(output.npixels,sizeof(short));
  if (output.mask==NULL) reportevt(flag_verbose,STATUS,5,
		  "Calloc of output.mask failed");
  output.varim=(float *)calloc(output.npixels,sizeof(float));
  if (output.varim==NULL) reportevt(flag_verbose,STATUS,5,
                                    "Calloc of output.varim failed");
  tempimage.varim=(float *)calloc(output.npixels,sizeof(float));
  if (tempimage.varim==NULL) reportevt(flag_verbose,STATUS,5,
                                       "Calloc of tempimage.varim failed");
  for (i=0;i<output.npixels;i++) {
    output.image[i]=input.image[i];
    output.mask[i]=input.mask[i];
    output.varim[i]=input.varim[i];
    tempimage.varim[i]=input.varim[i];
  }
    
  
  
  
  /* Get Sky level */
  //  doSky(image,bpm,weight);
  sky.npixX = NXSIZE;
  sky.npixY = NYSIZE;
  sky.weightmin = WEIGHT_MIN;
  sky.minsample = MIN_SAMPLE;
  sky.minline = MIN_LINE;
  sky.saturated = SATURATED;
  sky.underflow = 0.0;

  // printf("sky values %d %d %f %d %f %f %f \n", sky.npixX,sky.npixY,sky.weightmin,sky.minsample,sky.minline,sky.underflow,sky.saturated);

  error = doSky(2,1,&sky,output.image,output.mask,output.varim);
  if (error)
    {
      sprintf(event,"Error in Computing sky level for image: %s", input.name);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }


  }
  /* *************************************** */
  /* *********  COSMIC RAY SECTION  ******** */
  /* *************************************** */
  /*******************************************/

    
  if (flag_crays) {
    if (flag_verbose) 
      {
        sprintf(event,"Masking cosmic rays in %s\n",input.name);
        reportevt(flag_verbose,STATUS,1,event);
      }
    
    /* *************************************** */
    /* ********  MASK COSMIC RAY PIXELS ****** */
    /* *************************************** */
    
    //Calculate weight map
    if (flag_verbose==3) reportevt(flag_verbose,STATUS,1,
                                   "Now Calculating Weight mask");

    doWeight(output.image,output.mask,tempimage.varim);

    //Spread saturated                                                                  
    spread(BADPIX_SATURATE,1,BADPIX_SATURATE,output.image,output.mask,output.varim);
  
    if (flag_verbose==3) reportevt(flag_verbose,STATUS,1,
                                   "Now masking cosmic ray pixels");

    doCosmic(output.image, output.mask, output.varim,tempimage.varim, crfract, crsig2);
    
    spread(BADPIX_CRAY, 1000.0, BADPIX_CRAY, output.image, output.mask, output.varim);

    dogrow(BADPIX_CRAY,1,BADPIX_CGROW,output.image,output.mask,output.varim);
    

    /*
    doStat(output.image,output.mask,output.varim);
    
    
    if (flag_verbose) {
      sprintf(event,"Number of cosmic ray pixels masked= %d",numcrays);
      reportevt(flag_verbose,QA,1,event);
    }
    */

    printf(" Finished cosmic ray masking \n");
  }
  

  /* ******************************************** */
  /* **********   WRITE OUTPUT IMAGE  *********** */
  /* ******************************************** */


  /* Start building command line used*/
  strcpy(updated_command_line,argv[0]);
  strcat(updated_command_line," ");
  strcat(updated_command_line,input.name);

  if (flag_crays) {
    strcat(updated_command_line, " -crays");
  }
  if (flag_crfract) {
    sprintf(event," -crfract %.2f",crfract);
    strcat(updated_command_line,event);
  }
  if (flag_crsig2) {
    sprintf(event," -crsig2 %.2f",crsig2);
    strcat(updated_command_line,event);
  }


  sprintf(event,"Writing results to %s",output.name);
  reportevt(flag_verbose,STATUS,1,event);

  /* make sure path exists for new image */
  if (mkpath(output.name,flag_verbose)) {
    sprintf(event,"Failed to create path to file: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0);
  }
  else {
    sprintf(event,"Created path to file: %s",output.name);
    reportevt(flag_verbose,STATUS,1,event);
  }

  output.fptr = input.fptr;
  /* write the corrected image*/
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.image,&status)) {
    sprintf(event,"Writing image failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
	  
  /* Write information into the header describing the processing */
  /* get system time */
  tm=time(NULL);
  sprintf(comment,"%s",asctime(localtime(&tm)));
  comment[strlen(comment)-1]=0;
  if (fits_write_key_str(output.fptr,"DESCRMSK",comment, "Masked the image",&status)) {
    sprintf(event,"Writing processing history failed: %s", output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  sprintf(longcomment,"DESDM:");
  //  for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);
  sprintf(longcomment,"%s %s",longcomment,updated_command_line);
  reportevt(flag_verbose,STATUS,1,longcomment);
  if (fits_write_history(output.fptr,longcomment,&status)) {
    sprintf(event,"Writing longcomment failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  fits_movrel_hdu(output.fptr, 1, NULL, &status);
  /*if (fits_write_img(output.fptr,TUSHORT,1,output.npixels,output.mask, &status)) {*/
  if (fits_write_img(output.fptr,TUSHORT,1,output.npixels,output.mask, &status)) {
    sprintf(event,"Writing image mask failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  if (fits_update_key_str(output.fptr,"DES_EXT","MASK", "Extension type",&status)) {
    reportevt(flag_verbose,STATUS,5,"Setting DES_EXT=MASK failed");
    printerror(status);
  }

  /* now store the variance image that has been created or updated */
  fits_movrel_hdu(output.fptr, 1, NULL, &status);
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.varim, &status)) {
    sprintf(event,"Writing variance image failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  if (fits_update_key_str(output.fptr,"DES_EXT","WEIGHT", "Extension type",&status)) {
    sprintf(event,"Writing DES_EXT=WEIGHT failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }


  /* close the corrected image */
  if (fits_close_file(output.fptr,&status)) {
    sprintf(event,"File close failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }


  printf("End of maskcosmics program");      
        
  printf("\n");
  return(0);
}



void dogrow(int masksrch,int npix,int maskset,float *image,short *bpm,float *weight)
{
  //A routine to mask pixels that are within npix of a pixel that has one or more bits
  //set according to masksrch (i.e., the and with masksrch is non-zero)  
  int n;
  int in;
  int ix, iy;
  int ixl, iyl, ixu, iyu;
  int kount[2]={0,0};
  int mask2;

  mask2 = masksrch | maskset;
  for (n=0;n<naxis1*naxis2;++n)
    {
      if ((bpm[n]&masksrch)==0) continue;
      ++kount[0];
      iy = n/naxis1;
      ix = n - naxis1*iy;
      ixl = ix - npix;
      if (ixl<0) ixl = 0;
      iyl = iy - npix;
      if (iy<0) iyl = 0;
      ixu = ix + npix;
      if (ixu>=naxis1) ixu = naxis1 - 1 ;
      iyu = iy + npix;
      if (iyu>=naxis2) iyu = naxis2 - 1;
      for (ix=ixl;ix<=ixu;++ix)
        {
          for (iy=iyl;iy<=iyu;++iy)
            {
              in = iy*naxis1 + ix;
              if (bpm[in]&mask2) continue;
              ++kount[1];
              /* this is a hack. It originally sets the bpm bit to "maskset"
                 but we don;t want to add a new bit to the badpixel mask
                 and want to mark neighbour cosmic rays as BADPIX_CRAY.
                 bpm[in] |= maskset;
                 In order to do this, I am hacking the above assigment to
                 masksrch instead of maskset (how originally was) */
              bpm[in] |= maskset;
              weight[in] = 0.0;
            }
        }
    }
  printf("Grow found %i pixels with mask=%#08x \n",kount[0],masksrch);
  // printf("Grow added %i pixels with mask=%#08x \n",kount[1],maskset);
  return;
}

int doStat(float *image,short *bpm,float *weight)
{
  //A routine to compute statistics for masked pixels and print the results to the log file                   
  int i, ip, ix, iy;
  int kount[32];
  int kbad=0;
  short bpmval;

  for (i=0;i<32;++i) kount[i] = 0;

  for (ix=0;ix<NXSIZE;++ix)
    {
      for (iy=0;iy<NYSIZE;++iy)
        {
          ip = iy*NXSIZE + ix;
          //Bad pixel mask                                                                                        
          bpmval = bpm[ip];
          if (bpmval) ++kbad;
          for (i=0;i<32;++i) if ((bpmval>>i) & 1) ++kount[i];
        }
    }

  printf("\n----Masked Pixel Summary---\n");
  printf("Total masked pixels=%i\n",kbad);
  printf("Bad pixels (bit 0)=%i\n",kount[0]);
  printf("Saturated pixels=%i\n",kount[1]);
  printf("Interpoled pixels=%i\n",kount[2]);
  printf("Underflow pixels=%i\n",kount[3]);
  printf("Cosmic ray pixels=%i\n",kount[4]);
  printf("Masked star pixels=%i\n",kount[5]);
  printf("Satellite trail pixels=%i\n",kount[6]);
  printf("Saturated spread pixels=%i\n",kount[7]);
  printf("Saturated grow pixels=%i\n",kount[8]);
  printf("Border pixels=%i\n",kount[9]);
  printf("NAN pixels=%i \n",kount[10]);
  printf("High pixels=%i \n",kount[11]);
  printf("Low pixels=%i \n",kount[12]);
  printf("Low weight pixels=%i \n",kount[13]);
  printf("makeWeight Cosmic ray=%i\n",kount[14]);
  printf("Cosmic Ray Spread=%i\n",kount[15]);
  printf("Cosmic Ray Grow=%i\n",kount[16]);
  printf("\n");

  return(0);
}


int doWeight(float *image,short *bpm,float *weight)
{
  //A routine to compute the pixel weight map and flag bad pixels. 
  //Bad pixels are one of the following 
  //1.  In the border region (BORDER key word)
  //2.  Not a number (NAN or INF)
  //3. Above a maximum (SATURATE key word)
  //4. Below a minimum (NSKYSIG standard dev's below sky level)
  //5. Extremely low weight (less than MINWEIGHT) 

  int i;
  int ip, ix, iy;
  int bpmval;
  float imval;
  double gain, skylvl, sigma;
  for (ix=0;ix<NXSIZE;++ix)
    {
      for (iy=0;iy<NYSIZE;++iy)
        {
          ip = iy*NXSIZE + ix;
          skylvl = sky.skysol[0] + sky.skysol[1]*ix + sky.skysol[2]*iy;
          if (ix<BORDER || iy<BORDER) bpm[ip] |= BADPIX_BPM;
            if (ix>=NXSIZE-BORDER || iy>=NYSIZE-BORDER) bpm[ip] |=  BADPIX_BPM;
          //Image value                                                        
          imval = image[ip];

          //Check for INF and NAN values
          if (0.0*imval==0.0)
            {
              //Set a floor on the minimum value (for error calculation)
              if (imval<skylvl) imval=skylvl;
              //Saturated pixels -> weight=0 (regardless of bpm)
              if (imval>SATURATED)
                {
                  bpm[ip] |= BADPIX_SATURATE;
                }
            }
          else
            {
              image[ip] = 0.0;
              imval = skylvl;
              bpm[ip] |= BADPIX_LOW;
            }

          //Good pixels have a minimum weight
          if (weight[ip]<WEIGHT_MIN)
            {
              if (bpm[ip]==0) bpm[ip] |= BADPIX_LOW;
            }

          //Get local gain value
          gain = sky.gainsol[0] + sky.gainsol[1]*ix + sky.gainsol[2]*iy;
          //Now look at bpm 
          bpmval = bpm[ip];
          if (bpmval==0)
            {
              //Normal case: pixel not masked
              weight[ip] = gain/imval;
              sigma = sqrt(1.0/weight[ip]);
              if (image[ip]<skylvl-LOW_FACT*sigma) bpm[ip] |= BADPIX_LOW;
            }
          else weight[ip] = 0.0;

        }
    }
  return(0);
}




int doCosmic(float *image,short *bpm,float *weight,float *tempweight, double crfract, double crsig2)
{
  int ix, iy, ip;
  double gx, gy, gtot, phi;
  double skylvl;
  int kount=0;
  double p0, px, py;
  //double crfract;
  //double crsig2;
  double sigx2, sigy2;
  double crfract2, ofract;
  int debug;

  //crsig2 = CR_SIG2;
  //crfract = CR_FRACT;

  crfract2 = crfract*crfract;
  ofract = 1.0 - crfract;
  //  int debug;                                                                                              
  //  FILE *fhist;                                                                                            
  // fhist = fopen("fhist.dat","w"); 
  //  printf("Start cosmic.\n");
                                                                              
  for (iy=16;iy<NYSIZE-15;++iy)
    {
      for (ix=16;ix<NXSIZE-15;++ix)
        {
          ip = iy*NXSIZE + ix;
          if (bpm[ip]) continue;
          if (bpm[ip+1]) continue;
          if (bpm[ip+NXSIZE]) continue;
          p0 = image[ip];
          px = image[ip+1];
          py = image[ip+NXSIZE];
          skylvl = sky.skysol[0] + ix*sky.skysol[1] + iy*sky.skysol[2];
          if (p0<skylvl) p0=skylvl;
          if (px<skylvl) px=skylvl;
          if (py<skylvl) py=skylvl;
          if (px>p0)
            {
              gx = crfract*px - p0 + ofract*skylvl;
              sigx2 = crfract2/tempweight[ip+1]+1.0/tempweight[ip];
            }
          else
            {
              gx = crfract*p0 - px + ofract*skylvl;
              sigx2 = crfract2/tempweight[ip]+1.0/tempweight[ip+1];
            }
          if (py>p0)
            {
              gy = crfract*py - p0 + ofract*skylvl;
              sigy2 = crfract2/tempweight[ip+NXSIZE]+1.0/tempweight[ip];
            }
          else
            {
              gy = crfract*p0 - py + ofract*skylvl;
              sigy2 = crfract2/tempweight[ip]+1.0/tempweight[ip+NXSIZE];
            }
          if (gx<0.0) gx=0.0;
          if (gy<0.0) gy=0.0;

          gtot = gx*gx/sigx2 + gy*gy/sigy2;
          /*                                                                                  
      debug = 0;                                                                                              
      if (ix>=362 && ix<=368 && iy>=2159 && iy<=2165) debug=1;                                                
      if (debug)                                                                                              
        {                                                                                                     
          printf("ix=%i iy=%i gx=%f gy=%f  \n",ix+1,iy+1,gx,gy);                                              
          printf("image %f %f %f \n",image[ip],image[ip+1],image[ip+NXSIZE]);                                 
          printf("gtot=%f sky=%f sigx2=%f sigy2=%f \n",gtot,skylvl,sigx2,sigy2);                              
          //printf("px=%f p0=%f crfract=%f ofract=%f \n",px,p0,crfract,ofract);                               
          //if (p0>px) px=p0;                                                                                 
          //if (p0>py) py=p0;                                                                                 
          //fprintf(fhist,"%f %f %f\n",sqrt(gtot),px,py);                                                     
        }                                                                                                     
          */

          if (gtot>crsig2)
            {
              ++kount;
              bpm[ip] |= BADPIX_CRAY;
              weight[ip] = 0;
            }
        }
    }
  printf("Found %i cosmic ray pixels.\n",kount);
  return(0);
}

void spread(int masksrch,double sig,int maskset,float *image,short *bpm,float *weight)
{

  //Takes the set of points that have a bit in masksrch set                                                   
  //Find all touching points more than sig = # sigma above background at set bpm bit(s) maskset               

  int ix, iy;
  int ixl, iyl, ixu, iyu;
  int n;
  int in;
  int nstack = 0;
  int mask2;
  double skylvl;
  double pixval, sigsq;
  int nmask=0;
  /*int naxis1,naxis2;*/
  char event[1000];
  int flag_verbose=1;
  int status=0;

  sigsq = sig*sig;
  mask2 = masksrch | maskset;


  for (n=0;n<NSIZE;++n)
    {
      if (!(bpm[n]&masksrch)) continue;
      stack[nstack] = n;
      ++nstack;
    }
  printf("Spread found %i masked pixels with mask=%#08x.\n",nstack,masksrch);


  for (n=0;n<nstack;++n)
    {
      iy = stack[n]/naxis1;
      ix = stack[n] - naxis1*iy;
      ixl = ix - 1;
      if (ixl<=0) ixl = 0;
      iyl = iy - 1;
      if (iyl<=0) iyl = 0;
      ixu = ix + 1;
      if (ixu>=naxis1) ixu = naxis1 - 1;
      iyu = iy + 1;
      if (iyu>=naxis2) iyu = naxis2 - 1;
      for (ix=ixl;ix<=ixu;++ix)
        {
          for (iy=iyl;iy<=iyu;++iy)
            {
              in = iy*naxis1 + ix;
              if (bpm[in]&mask2) continue;

              skylvl = sky.skysol[0] + sky.skysol[1]*ix + sky.skysol[2]*iy;

              pixval = image[in] - skylvl;             
              if (pixval<=0.0) continue;
              if (pixval*pixval*weight[in]<sigsq) continue;
              bpm[in] |= maskset;
              weight[in] = 0.0;
              ++nmask;
              stack[nstack] = in;
              ++nstack;
              //Stack overflow indicates a bug.  Should never happen                                    
              /*if (nstack>NSIZE) errExit(flog,"Stack overflow in spread. \n",20);*/
              if (nstack>NSIZE) {
                sprintf(event,"Stack overflow in spread \n");
                reportevt(flag_verbose,STATUS,5,event);
                printerror(status);
                }
            }
        }
  printf("Spread masked %i pixels above %f sigma with mask=%#08x \n",nmask,sig,maskset);
  return;
    }
}

void   creategrid(x1,y1,pixval,output,mean,std)  	    
     float x1,y1,**pixval,*mean,*std;
     desimage *output ;
{
  int n,num=0,pixnum,numval[10][10],xval,yval,i,j,pixbin,dx,dy ;
  int xpos,ypos;
  float *region;

  region = (float *)calloc(output->npixels,sizeof(float));
  /*x1,y1 are positions around which interpolations need to be done*/
  for (i=0;i<=9;i++) {
    for (j=0 ; j<=9 ; j++) {
      pixval[i+1][j+1]  = 0 ;
      numval[i][j] = 0 ;
    }
  }
  /* get the maximum and minimum values to be used for looping around the given point  */

  for (dy=-50;dy<50;dy++)
    for (dx=-50;dx<50;dx++)
      {
	/* only use pixel values which have cosmic ray mask set */
	/*      do bitwise Nanding of output.mask and BADPIX_CRAY */
	/* get only use pixel values whose X and Y-position is +/- 50 units around (x1,y1)  */
	xpos = dx + x1 ;
	ypos = dy + y1 ;
	xval = (xpos-x1+50)*10/100 ;
	yval = (ypos-y1+50)*10/100 ;
	i = (int) (ypos*output->axes[0] + xpos) ;
	if ((xval>9) || (yval>9))printf("Warning %d %d\n",xval,yval); 
	if (xval >9) xval = 9 ;
	if (yval >9) yval = 9 ;
	if (!(output->mask[i] & BADPIX_CRAY) & (output->image[i] >10) ) 
	  if ((i>=0) && (i < output->npixels)) { 
	    *mean += output->image[i] ;
	    region[num] = output->image[i] ;
	    num+= 1 ;
	    pixval[xval+1][yval+1] += output->image[i] ;
	    numval[xval][yval] += 1 ;
	  }		
      }
  *mean /= num ;
  *std = 0 ;
  for(i=0;i<num;i++)*std += pow((region[i]-*mean),2);
  *std = sqrt(*std)/(num-1) ;  
  for (i=0;i<=9;i++) {
    for (j=0 ; j<=9 ; j++) {
      if (numval[i][j] > 0)pixval[i+1][j+1]   /= numval[i][j] ;
    }
  }
  free(region);
}
 

void getfloatheader(image,inputval,value,mode)
     desimage *image;
     int mode;
     char *inputval ;
     float *value;
{
  char currentimage[1000] ;
  static int status=0;
  int hdutype ;

  sprintf(currentimage,"%s",image->name);
  if (!fits_open_image(&image->fptr, currentimage, 
		       mode, &status))
    if (fits_read_key_flt(image->fptr,inputval,value,NULL,
			  &status))
      {*value=0; }
  /* fits_close_file(image->fptr,&status) ; */
  return ;
}


void headervalue(image,inputval,value,mode)
     desimage *image;
     int mode;
     char *inputval,*value ;
{
  char currentimage[1000] ;
  static int status=0;
  int hdutype ;
  sprintf(currentimage,"%s",image->name);
  if (!fits_open_image(&image->fptr, currentimage, 
		       mode, &status))
    if (fits_read_key_str(image->fptr,inputval,value,NULL,
			  &status))
      {sprintf(value,"undefined"); }

  /* fits_close_file(image->fptr,&status) ; */
  return ;
}


/*main code */

/*int main(int argc,char *argv[])
{
  return(MakeMask(argc,argv));
}
*/




/* Routines (doSky and fitsky, skysub) are from SNe DIffimg pipeline */

int doSky(int verbose,int dogain,skypar *sky,float *image,short *bpm,float *weight)
{
  //A routine to calculate the sky background assuming a polynomial form
  //A linear form is assumed, but it could easily be extended to higher order, if necessary
  int gindex[NSKYSAMP*NSKYSAMP];
  double pixel[NSKYSAMP*NSKYSAMP];
  int ix, iy, i, j;
  int in, io;
  int nsample;
  double err;
  double skymap[NYSKY][NXSKY], skysig[NYSKY][NXSKY], gain[NYSKY][NXSKY];
  //  double skyz;
  double x, y;
  int n;
  int ie;
  double ave, sigtot, sigma1, sigma2;
  double lower, upper;
  int nl, nu;
  int ixref, iyref;
  int nxp, nyp, nxs, nys;
  int nyg;
  int ntotp;
  double rms;
  int npixX, npixY;
  double wgtmin;
  int minsamp;
  //  double Y[3];
  int error;
  float imval;
  float saturated;
  float underflow;


  npixX = sky->npixX;
  npixY = sky->npixY;
  wgtmin = sky->weightmin;
  minsamp = sky->minsample;
  saturated = sky->saturated;
  underflow = sky->underflow;

  // printf("inside sky saturated %f\n",saturated);

  nxp = npixX/NXSKY;
  nyp = npixY/NYSKY;
  nxs = nxp/NSKYSAMP;
  nys = nyp/NSKYSAMP;
  ntotp = 0;
  sigtot = 0.0;


  for (iy=0;iy<NYSKY;++iy)
    {
      iyref = iy*nyp;
      for (ix=0;ix<NXSKY;++ix)
	{
	  ixref = ix*nxp;
	  io = 0;
	  for (j=0;j<nyp;j+=nys)
	    {
	      for (i=0;i<nxp;i+=nxs)
		{
		  in = npixX*(iyref+j) + ixref+i;
		  //first check image for NAN or INF values
		  imval = image[in];
		  if (isnan(imval)) continue;
                  if (isinf(imval)) continue;
		  if (imval>saturated) continue;
		  if (imval<underflow) continue;
		  if (weight[in]<=wgtmin) continue;
		  if (bpm[in]) continue;
		  pixel[io] = image[in];
		  gindex[io] = io;
		  ++io;
		}
	    }
	  nsample = io;
	  if (nsample<minsamp) 
	    {
	      skymap[iy][ix] = 0.0;
	      skysig[iy][ix] = 0.0;
	      gain[iy][ix] = 1.0;
	      continue;
	    }
      
	  Qsort(gindex,pixel,nsample-1);

	  n = 0.5*nsample;
	  ie = gindex[n];
	  ave = pixel[ie];
	  //	  printf("Pixel average ie=%i %f\n",ie,pixel[ie]);
	  n = 0.1587*nsample;
	  ie = gindex[n];
	  sigma1 = ave - pixel[ie];
	  //printf("Lower ie=%i %f\n",ie,pixel[ie]);
	  //Refine estimate
	  lower = ave - 4.0*sigma1;
	  for (n=0;n<nsample;++n)
	    {
	      ie = gindex[n];
	      if (pixel[ie]>lower) break;
	    }
	  nl = n;
	  upper = ave + 4.0*sigma1;
	  for (n=nsample-1;n>0;--n)
	    {
	      ie = gindex[n];
	      if (pixel[ie]<upper) break;
	    }
	  nu = n;
	  //	  printf(" New lower limit=%i new upper=%i\n",nl,nu);
	  n = 0.5*(nl+nu);
	  ie = gindex[n];
	  ave = pixel[ie];
	  
      //printf("ie=%i %f\n",ie,pixel[ie]);
	  
      n = 0.1587*(nu-nl) + nl;
      ie = gindex[n];
	  sigma1 = ave - pixel[ie];
	  
      // printf("ie=%i %f\n",ie,pixel[ie]);
	  
      n = 0.8413*(nu-nl) + nl;
	  ie = gindex[n];
	  sigma2 = pixel[ie] - ave;
	  err = sigma1/sqrt(nu-nl);

	  // printf("err=%i %f %f %f\n",ie,err,sigma1,sigma1*sigma1);
	  
	  skymap[iy][ix] = ave;
	  skysig[iy][ix] = sigma1;
	  //This is estimate of the gain based on sky fluctuations.  It should be correct if:
          // 1.  The quantum efficiency is the same for all pixels
          // 2.  The fluctuations are dominated by photon statistics
          // For a more precise calculation, one needs the raw exposures
	  gain[iy][ix] = ave/(sigma1*sigma1);
             
      //if (verbose>=2) printf ("x=%i y=%i sky=%.2f +/- %.2f lower-sigma=%.2f upper-sigma=%.2f gain=%.3f\n",
	  //	   ix,iy,skymap[iy][ix],err,skysig[iy][ix],sigma2,gain[iy][ix]);
	  ++ntotp;
	  sigtot += (sigma1*sigma1);
           
	}
    }
  sky->sigma=sqrt(sigtot/ntotp);
  //Fit for sky level
  sky->itype=0;

 
  error = fitsky(sky,&skymap[0][0],&skysig[0][0]);
  //printf("fitsky\n");


  if (error) return(error);
  if (verbose>=2)
    {
      printf("Sky sigma=%.2f sky0=%.2f skyx=%.2e skyy=%.2e rms=%.2f\n",
			      sky->sigma,sky->skysol[0],sky->skysol[1]*npixX,sky->skysol[2]*npixY,sky->skysol[3]);
    }

  //Fit for effective gain
  if (dogain)
    {
      sky->itype=1;
      error = fitsky(sky,&gain[0][0],&skysig[0][0]);
      if (verbose>=2)
	{
	  printf("Sky gain0=%.2f gainx=%.2e gainy=%.2e rms=%.2f\n",
		  sky->gainsol[0],sky->gainsol[1]/npixX,sky->gainsol[2]/npixY,sky->gainsol[3]);
	}
    }

  return(error);
}

void skysub(skypar *sky,float *image)
{
  int npixX, npixY;
  int ix, iy, in;
  double x, y;

  printf("skysub\n");

  npixX = sky->npixX;
  npixY = sky->npixY;

  for (iy=0;iy<npixY;++iy)
    {
      y = iy;
      for (ix=0;ix<npixX;++ix)
	{
	  x = ix;
	  in = npixX*iy + ix;
	  image[in] = image[in] - sky->skysol[0] - x*sky->skysol[1] - y*sky->skysol[2];
	}
    }

  return;
}

int fitsky(skypar *sky,double *data,double *sigma)
{
  int i, j;
  int ix, iy;
  double x, y;
  double chisq, d, dev, rms;
  double desc;
  int kx, ky;
  int ndata;
  int minline;

  double H[3][3], Y[3];
  int indx[10];
  double det;
  int icon;


  minline = sky->minline;

  //Check for enough x values  
  kx = 0;
  for (ix=0;ix<NXSKY;++ix)
    {
      ky = 0;
      for (iy=0;iy<NYSKY;++iy) if (sigma[iy*NXSKY+ix]>0.0) ++ky;
      if (ky>0) ++ kx;
     }
  if (kx<minline)
    {
      printf("Found %i x values for sky fit. Required %i\n",kx,minline);
      return(1);
    }

  //Check for enough y values
  ky = 0;
  for (iy=0;iy<NYSKY;++iy)
    {
      kx = 0;
      for (ix=0;ix<NXSKY;++ix) if (sigma[iy*NXSKY+ix]>0.0) ++kx;
      if (kx>0) ++ ky;
     }
  if (ky<minline)
    {
      printf("Found %i y values for sky fit. Required %i\n",ky,minline);
      return(2);
    }

  for (i=0;i<3;++i)
    {
      Y[i] = 0.0;
      for (j=0;j<3;++j) H[i][j]=0.0;
    }

  for (iy=0;iy<NYSKY;++iy)
    {
      y = iy + 0.5;
      for (ix=0;ix<NXSKY;++ix)
	{
	  if (sigma[iy*NXSKY+ix]==0.0) continue;
	  x = ix + 0.5;
	  H[0][0] += 1.0;
	  H[0][1] += x;
	  H[0][2] += y;
	  H[1][1] += x*x;
	  H[1][2] += x*y;
	  H[2][2] += y*y;
	  d = data[iy*NXSKY+ix];
	  Y[0] += d;
	  Y[1] += x*d;
	  Y[2] += y*d;
	}
    }
  H[1][0] = H[0][1];
  H[2][0] = H[0][2];
  H[2][1] = H[1][2];
  //Solve for parameters 
  ludcmp(&H[0][0],3,3,indx,&desc,&icon);
  lubksb(&H[0][0],3,3,indx,Y);
  
  chisq = 0.0;
  ndata = 0;
  for (iy=0;iy<NYSKY;++iy)
    {
      y = iy + 0.5;
      for (ix=0;ix<NXSKY;++ix)
	{
	  if (sigma[iy*NXSKY+ix]==0.0) continue;
	  x = ix + 0.5;
	  dev = data[iy*NXSKY+ix]-Y[0]-Y[1]*x-Y[2]*y;
	  chisq += (dev*dev);
	  ++ndata;
	}
    }
  rms = sqrt(chisq/ndata);
  if (sky->itype==0)
    {
      sky->skysol[0] = Y[0];
      //Normalize slope from per bin to per pixel
      sky->skysol[1] = Y[1]*NXSKY/sky->npixX;
      sky->skysol[2] = Y[2]*NYSKY/sky->npixY;
      sky->skysol[3] = rms;
    }
  else
    {
      sky->gainsol[0] = Y[0];
      //Normalize slope from per bin to per pixel
      sky->gainsol[1] = Y[1]*NXSKY/sky->npixX;
      sky->gainsol[2] = Y[2]*NYSKY/sky->npixY;
      sky->gainsol[3] = rms;
    }
  return(0);
}



void ludcmp(double* a, const int n, const int ndim, int* indx,
            double* d, int* icon)
{
  /* Local variables */
  int a_dim1, a_offset, i1, i2, i3;
  int imax = 0;
  int  i, j, k;
  double aamax, dum, sum;
  double vv[1000] = { 0.0 };

  /* Function Body */
  *icon = 1;
  *d = 1.0;


  for (i = 0; i < n; ++i)
    {
      aamax = 0.0;

      for (j = 0; j < n; ++j)
        {
          if (fabs(a[i + j * ndim]) > aamax)
            {
              aamax = fabs(a[i + j * ndim]);
            }
        }

      if (aamax == 0.0)
        {
          printf("LU decomposition (ludcmp.c) : Singular matrix, but continue. %d %d %d %d %d\n",i,j,i1,i2,ndim);
          *icon = -10;
          return;

        }

      vv[i] = 1.0 / aamax;
    }
  i1 = n;
  for (j = 0; j < n; ++j)
    {
      if (j > 0)
        {
          for (i = 0; i < j; ++i)
            {
              sum = a[i + j * ndim];
              if (i > 0)
                {
                  i3 = i - 1;
                  for (k = 0; k < i; ++k)
                    {
                      sum -= a[i + k * ndim] * a[k + j * ndim];
                    }
                  a[i + j * ndim] = sum;
                }
            }
        }

      aamax = 0.0;

      for (i = j; i < n; ++i)
        {
          sum = a[i + j * ndim];
          if (j > 0)
            {

              for (k = 0; k < j; ++k)
                {
                  sum -= a[i + k * ndim] * a[k + j * ndim];
                }
              a[i + j * ndim] = sum;
            }

          dum = vv[i] * fabs(sum);
          if (dum >= aamax)
            {
              imax = i;
              aamax = dum;
            }
        }

      if (j != imax)
        {

          for (k = 0; k < n; ++k)
            {
              dum = a[imax + k * ndim];
              a[imax + k * ndim] = a[j + k * ndim];
              a[j + k * ndim] = dum;
            }
          *d = - *d;
          vv[imax] = vv[j];
        }

      indx[j] = imax;
      if (j != (n-1))
        {
          if (fabs(a[j + j * ndim]) <= 1.e-20)
            {
              if (a[j + j * ndim] > 0.0)
                {
                  a[j + j * ndim] = 1.0e-20;
                }
              else
                {
                  a[j + j * ndim] = -1.0e-20;
                }
            }

          dum = 1.0 / a[j + j * ndim];

          for (i = j + 1; i < n; ++i)
            {
              a[i + j * ndim] *= dum;
            }
        }
    }

  if (fabs(a[n-1 + (n-1) * ndim]) <= 1.0e-20)
    {
      if (a[n-1 + (n-1) * ndim] > 0.0)
        {
          a[n-1 + (n-1) * ndim] = 1.0e-20;
        }
      else
        {
          a[n-1 + (n-1) * ndim] = -1.0e-20;
        }
    }
  return;
} /* ludcmp */


void lubksb(const double* a, const int n, const int ndim,
            const int* indx, double* b)
{

  /* Local variables */
  int i1, i2;
  int i, j, ii, ll;
  double sum;


  /* Function Body */
  ii = -1;
  i1 = n;
  for (i = 0; i < i1; ++i)
    {
      ll = indx[i];
      sum = b[ll];
      b[ll] = b[i];
      if (ii != -1)
        {
          for (j = ii; j < i; ++j)
            {
              sum -= a[i + j * ndim] * b[j];
            }
        }
      else
        {
          if (sum != 0.0)
            {
              ii = i;
            }
        }
      b[i] = sum;
    }

  for (i = n-1; i >= 0; --i)
    {
      sum = b[i];
      if (i < n-1)
        {

          for (j = i + 1; j < n; ++j)
            {
              sum -= a[i + j * ndim] * b[j];
            }
        }
      b[i] = sum / a[i + i * ndim];
    }
  return;
} /* lubksb */


void Qsort(int gindex[],double value[], int last)
{
  //Quick sort routine
  //Sort according to value.  Last element in array is value[last]
  //On entry gindex[n] = n.  On exit value[gindex[n]] will be ordered for n=0,1,2,... 

  int stack_pointer = 0;
  int first_stack[32];
  int last_stack[32];
  int ifirst, ilast, imed, idown, iup;
  int first=0;
  for (;;)
    {
      if (last - first <= INSERTION_SORT_BOUND)
        {
          /* for small sort, use insertion sort */
          int indx;
          int prev_val = gindex[first];
          int cur_val;

          for (indx = first + 1; indx <= last; ++indx)
            {
              cur_val = gindex[indx];
              if (value[prev_val]>value[cur_val])
                {
                  /* out of order: array[indx-1] > array[indx] */
                  int indx2;
                  gindex[indx] = prev_val; /* move up the larger item first */

                  /* find the insertion point for the smaller item */
                  for (indx2 = indx - 1; indx2 > first; )
                    {
                      int temp_val = gindex[indx2 - 1];
                      if (value[temp_val]>value[cur_val])
                        {
                          gindex[indx2--] = temp_val;
                          /* still out of order, move up 1 slot to make room */
                        }
                      else
                        break;
                    }
                  gindex[indx2] = cur_val; /* insert the smaller item right here */
                }
              else
                {
                  /* in order, advance to next element */
                  prev_val = cur_val;
                }
            }
        }
      else
        {
          int pivot;

          /* try quick sort */
          {
            int temp;
            int med = (first + last) >> 1;
            /* Choose pivot from first, last, and median position. */
            /* Sort the three elements. */
            temp = gindex[first];
            ilast = gindex[last];
            if (value[temp]>value[ilast])
              {
                gindex[first] = gindex[last]; gindex[last] = temp;
              }
            temp = gindex[med];
            ifirst = gindex[first];
            if (value[ifirst]>value[temp])
              {
                gindex[med] = gindex[first]; gindex[first] = temp;
              }
            temp = gindex[last];
            imed = gindex[med];
            if (value[imed]>value[temp])
              {
                gindex[last] = gindex[med]; gindex[med] = temp;
              }
            pivot = gindex[med];
          }
          {
            int up;
            {
              int down;
              /* First and last element will be loop stopper. */
              /* Split array into two partitions. */
              down = first;
              up = last;
              for (;;)
                {
        do
          {
            ++down;
            idown = gindex[down];
          } while (value[pivot]>value[idown]);
        do
          {
            --up;
            iup = gindex[up];
          } while (value[iup]>value[pivot]);


        if (up > down)
          {
            int temp;
            /* interchange L[down] and L[up] */
            temp = gindex[down]; gindex[down]= gindex[up]; gindex[up] = temp;
          }
        else
          break;
                }
            }
            {
              int len1; /* length of first segment */
              int len2; /* length of second segment */
              len1 = up - first + 1;
              len2 = last - up;
              /* stack the partition that is larger */
              if (len1 >= len2)
                {
                  first_stack[stack_pointer] = first;
                  last_stack[stack_pointer++] = up;

                  first = up + 1;
                  /*  tail recursion elimination of                                                                         
                   *  Qsort(gindex,fun_ptr,up + 1,last)                                                                                         */
                }
              else
                {
                  first_stack[stack_pointer] = up + 1;
                  last_stack[stack_pointer++] = last;

                  last = up;
                  /* tail recursion elimination of                                                                                              * Qsort(gindex,fun_ptr,first,up)                                                                                             */
                }
            }
            continue;
          }
          /* end of quick sort */
        }
      if (stack_pointer > 0)
        {
          /* Sort segment from stack. */
          first = first_stack[--stack_pointer];
          last = last_stack[stack_pointer];
        }
      else
        break;
    } /* end for */
  return;
}
