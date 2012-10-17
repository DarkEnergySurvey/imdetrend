/*$Id$*/
/*
 mkmask

 	takes an input MEF image (image, bpm, weight map)
	along with an SExtractor OBJECT image marking the cosmic ray
	pixels and sets the cosmic ray bits in the bpm as well as interpolating
	over the cosmic ray pixels to replace them.

	Future plans:
    * Could also take an input star list (from USNO-B?) and apply bright
      star masks.
    * Could use a mapping for scattered light to mask light from stars that
      are off the field of view.
    * Could carry out a Hough transform to look for and remove satellite
      trails.

*/

#include "imsupport.h"
#include <getopt.h>

#include "argutils.h"

#define LINEAR_INTERP 1
#define MOFFAT_INTERP 2
#define GAUSS_INTERP 3

static const char *svn_id = "$Id$";

void print_usage(char *program_name)
{
  printf("%s <input image> <options>\n",program_name);
  printf("  Masking Options\n");
  printf("    -crays <cosmicrayimage>\n");
  printf("    -srcgrowrad <growrad>\n");
  printf("    -stars <astrostdsfile> \n");
  printf("    -flag_horiztrails \n"); 
  printf("    -nointerpolate \n");
  printf("  Output Options\n");
  printf("    -output <newimage>\n");
  printf("    -verbose <0-3>\n");
  printf("    -help    (print usage and exit)\n");
  printf("    -version (print version and exit)\n");
}

static int flag_nointerp = NO;
static int flag_horiz    = NO;

int MakeMask(int argc,char *argv[])
{
  char	filter[100]="",comment[1000],longcomment[10000],event[10000],
    scaleregion[100],keycard[100],outputname[1000],astrostdsfile[10000],
    imagetype[1000],exposstr[1000],command_line[1000];
  char outputnamewithtempl[2500],inname_temp[1000];
  int	i,j,x,y,k,l,flag_fringe=NO,flag_illum=YES,xp,yp,th,horiztrails=0,
    minsize=4,maxsize=10,loc,dx,dy,xmin,xmax,ymin,ymax,scalenum,
    count,flag_verbose=1,
    nvec,xlen,ylen,totpix,mkpath(),ccdnum=0,ncompare=0,
    flag_illumcor_compare=NO,flag_fringecor_compare=NO,
    scaleregionn[4]={500,1500,1500,2500},ctr,
					   keysexist,numcrays,interp,xpospix,ypospix,
					   pixrow,pixcol,xval,yval;
  short   *grow;
  static int status=0;
  long	ranseed=0,seed = -30;
  float	*vecsort,ran1(),subsample_fraction,*randnum,offset,rms,maxdev,exposure,
    fwhm,scalefactor,*scalesort,mean=0,krad,growrad,gasdev(),std,sigma,regionmode,
    goodpixx[10],goodpixy[10],**pixval,outy,erry,badrow,badcolumn,xmask,ymask;
  double	value;
  desimage input,cray,output;	
  void	rd_desimage(),shell(),reportevt(),image_compare(),
    retrievescale(),printerror();
  void    polin2(),creategrid(),mkstarmask(),headervalue(),pixelhisto(),getexposuretime();	
  time_t	tm;
  float	interp_noise=1.0,interp_fwhm=0.0;
  int	flag_output=0,flag_crays=0,flag_stars=0,flag_srcgrowrad=0,
    num = 0,filemode=0,goodpix,badpix;
  enum {OPT_CRAYS=1,OPT_SRCGROWRAD,OPT_STARS,OPT_FLAG_HORIZTRAILS,OPT_NOINTERPOLATE,
	OPT_OUTPUT,OPT_VERBOSE,OPT_HELP,OPT_VERSION};
  
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
	{"crays",            required_argument, 0,              OPT_CRAYS},
	{"srcgrowrad",       required_argument, 0,              OPT_SRCGROWRAD},
	{"stars",            required_argument, 0,              OPT_STARS},
        {"output",           required_argument, 0,              OPT_OUTPUT},
	{"verbose",          required_argument, 0,              OPT_VERBOSE},
	{"version",          no_argument,       0,              OPT_VERSION},
	{"help",             no_argument,       0,              OPT_HELP},
	{"nointerpolate",    no_argument,       &flag_nointerp, YES},
	{"flag_horiztrails", no_argument,       &flag_horiz,    YES},
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
    case OPT_CRAYS: // -crays
      sprintf(cray.name,"%s",optarg);
      flag_crays=YES;
      break;
    case OPT_SRCGROWRAD: // -srcgrowrad
      flag_srcgrowrad=YES;
      sscanf(optarg,"%f",&growrad);
      if (growrad>5) {
	sprintf(event,"Growrad value of %f is too large.  Resetting growrad to 5.\n",growrad);
	reportevt(flag_verbose,STATUS,3,event);
	growrad=5;
      }			
      break;
    case OPT_STARS: // -stars
      sprintf(astrostdsfile,"%s",optarg);
      flag_stars=YES;
      break;
    case OPT_OUTPUT: // -output
      sprintf(outputname,"%s",optarg);
      flag_output=YES;
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
  if (!flag_output || !strcmp(input.name,outputname)) { 
    rd_desimage(&input,READWRITE,flag_verbose);
  }
  else {
    rd_desimage(&input,READONLY,flag_verbose);
  }
  fits_file_mode(input.fptr,&filemode,&status);
  headervalue(&input,"FILTER",filter,filemode);
  sprintf(event,"Input image FILTER=%s and CCDNUM=%d: %s",filter,ccdnum, input.name);
  reportevt(flag_verbose,STATUS,1,event);
  
  headervalue(&input,"OBSTYPE",imagetype,filemode);
  getexposuretime(&input,"EXPTIME",&exposure,filemode);
  /*	printf ("OBSTYPE=%s\n",imagetype); */ 
  if(!strcmp(imagetype,"remap") || (flag_horiz ==1))horiztrails = 1 ;
  
  /* create an output image with same parameters as input image */
  output=input;
  
  /* copy the correct output name into place */
  if (!flag_output || !strcmp(input.name,outputname)) { 
    sprintf(output.name,"%s",input.name);
  }
  else {
    sprintf(output.name,"%s",outputname);
    
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
    for (i=0;i<output.npixels;i++) {
      output.image[i]=input.image[i];
      output.mask[i]=input.mask[i];
      output.varim[i]=input.varim[i];
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
    
    /* create the output file using input file as template */
    strcpy(outputnamewithtempl, "!");
    strcat(outputnamewithtempl, output.name);
    strcat(outputnamewithtempl, "(");
    strcat(outputnamewithtempl, input.name);
    strcat(outputnamewithtempl, ")");
    if (fits_create_file(&output.fptr,outputnamewithtempl,&status)) {
      sprintf(event,"File creation failed: %s",outputnamewithtempl);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
  }
  /* do star masking*/
  if (flag_stars){
    if (flag_verbose)reportevt(flag_verbose,STATUS,1,"Implementing bright stars masking"); 
    mkstarmask(&input,&output,filter,astrostdsfile,horiztrails,exposure,flag_verbose,flag_nointerp);
  }
  /*	printf("finished star masking\n"); */
  /* *************************************** */
  /* *********  COSMIC RAY SECTION  ******** */
  /* *************************************** */
  /*******************************************/
  
  if (flag_crays) {
    if (flag_verbose) 
      {
	sprintf(event,"Masking cosmic rays in %s using OBJECT image %s",input.name,cray.name);
	reportevt(flag_verbose,STATUS,1,event);
      }
    /*	    sprintf("Masking cosmic rays in %s using OBJECT image %s",input.name,
      cray.name)); */
    rd_desimage(&cray,READONLY,flag_verbose);
    
    /* *************************************** */
    /* *******  READ COSMIC RAY IMAGE  ******* */
    /* *************************************** */
    if (flag_verbose==3) reportevt(flag_verbose,STATUS,1,
				   "Reading cosmic ray OBJECT image");
    rd_desimage(&cray,READONLY,flag_verbose);
    
    /* make sure same size as input image */
    if (input.axes[0]!=cray.axes[0] || input.axes[1]!=cray.axes[1]) {
      sprintf(event,"Cosmic ray OBJECT image %ldX%ld different size than input image:  %ldX%ld )",
	      cray.axes[0],cray.axes[1],input.axes[0],input.axes[1]);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    
    
    /* *************************************** */
    /* ********  MASK COSMIC RAY PIXELS ****** */
    /* *************************************** */
    if (flag_verbose==3) reportevt(flag_verbose,STATUS,1,
				   "Now masking cosmic ray pixels");
    
    numcrays=0;
    for (i=0;i<cray.npixels;i++) {
      if ((cray.image[i]>0.0) && (output.mask[i] ==0)) {
	output.mask[i] |= BADPIX_CRAY;
	output.varim[i] = 0.0 ;
	
      }
    }
    /*get median and variance of entire section of image */
    scaleregionn[0] =  50; 
    scaleregionn[1] = output.axes[0] -50 ;
    scaleregionn[2] = 50; ;
    scaleregionn[3] = output.axes[1] -51  ;
    if (scaleregionn[0]<0)scaleregionn[0]= 0;
    if (scaleregionn[0]>output.axes[0])scaleregionn[0]= output.axes[0];
    if (scaleregionn[2]<0)scaleregionn[2]= 0;
    if (scaleregionn[3]>output.axes[1])scaleregionn[3]= output.axes[1];
    scalenum=(scaleregionn[1]-scaleregionn[0])*
      (scaleregionn[3]-scaleregionn[2]);
    scalesort=(float *)calloc(scalenum,sizeof(float));
    pixelhisto(&output,scaleregionn,scalesort,flag_verbose,&scalefactor,&regionmode,&sigma) ;
    if (sigma<7.0)sigma=7.0 ;
    free(scalesort);


    grow = (short *)calloc(output.npixels,sizeof(short));
    /****************************/ 
    /** grow radius section ****/	
    /***************************/
	
    if ((flag_srcgrowrad==YES) && (flag_crays==YES)) {
      if (flag_verbose==3) 
	reportevt(flag_verbose,STATUS,1,
		  "Beginning source grow radius filtering"); 
      printf("starting grow radius filtering %f \n",growrad);
      for (i=1;i<output.npixels;i++) 
	if (output.mask[i] & BADPIX_CRAY) {
	  xp = i%output.axes[0] ;
	  yp = i/output.axes[0] ;
	  /*loop over circle*/	
	  if((xp>=2) && (yp>=2) && (abs(xp-output.axes[0])>2) && abs(yp-output.axes[1])>2)
	    {
	      for (krad = 0;krad<growrad;krad+=0.2)
		for(th= 0;th<=360;th+=1)
		  {
		    xmask = krad*cos(th*3.1416/180.0);
		    ymask = krad*sin(th*3.1416/180.0);        
		    xpospix = (int) (xmask + xp) ;
		    ypospix = (int) (ymask + yp) ;
		    j = ypospix*output.axes[0] + xpospix ;
		    /*do not set cosmic ray mask at this stage, else it will loop over all pixels in the image*/
		    if (!(output.mask[j] & BADPIX_CRAY))grow[j] = 1 ; 
		  }		
	    } }
    }		

    if (flag_crays) {
      /* ****************************************** */
      /* *********  INTERPOLATION SECTION  ******** */
      /* ****************************************** */
      /* pixels used for interpolating */

      /*   do cosmic ray interpolation on input.image */
      /* probably want to introduce an interpolated image */
      ctr = 0 ;
      for (i=0;i<output.npixels;i++) {
	/*       mask all pixels which have a grow-radius flag set */	
	if (!(output.mask[i] & BADPIX_CRAY) && (grow[i] == 1))
	  {output.mask[i] |= BADPIX_CRAY ; output.varim[i]= 0.0 ;} 
	/*      do bitwise anding of output.mask and BADPIX_CRAY */
	if (output.mask[i] & BADPIX_CRAY) {
	  numcrays++;
	  /*get  X and Y position of the bad pixel badrow = x position; badcolumn = y position */
	  badrow =     i%cray.axes[0] ;		  
	  badcolumn =  i/cray.axes[0] ;
	  /*		printf("Bad row %f %f\n",badrow,badcolumn); */
	  /*    pixval = matrix(1,10,1,10) ; */
	  /* get the absiccaa and oordinate values around the desired points for interpolation */
	  /*
	    for (k=0 ; k<=9; k++)
	    {		  
	    goodpixx[k] = 10*(k-4.5) + badrow   ;
	    goodpixy[k] = 10*(k-4.5) + badcolumn ;
	    } */
	  /*create a 10 X 10 grid around (badrow,badcolumn)  to be used for interpolation*/

	  ctr ++ ;
	  /*call routine to get grid values for interpolation */
	  /*			if ((badrow>1) && (badcolumn>1)) printf("calling creategrid %f %f\n",badrow,badcolumn); */
	  /*creategrid(badrow,badcolumn,pixval,&output,&mean,&std)  ; 	     */
	  /* Using numerical recipies routine polin2 to do the interpolation */
	  /*	polin2(goodpixx-1,goodpixy-1,pixval,10,10,badrow,
	    badcolumn,&outy,&erry); */ 
	  /*	output.image[i] = outy ; */
			
	  /*	scaleregionn[0] =  badrow  - 30; 
	    scaleregionn[1] = badrow   + 30 ;
	    scaleregionn[2] =  badcolumn - 30 ;
	    scaleregionn[3] =  badcolumn + 30   ;
	    if (scaleregionn[0]<0)scaleregionn[0]= 0;
	    if (scaleregionn[0]>output.axes[0])scaleregionn[0]= output.axes[0];
	    if (scaleregionn[2]<0)scaleregionn[2]= 0;
	    if (scaleregionn[3]>output.axes[1])scaleregionn[3]= output.axes[1];
	    scalenum=(scaleregionn[1]-scaleregionn[0])*
	    (scaleregionn[3]-scaleregionn[2]);
	    scalesort=(float *)calloc(scalenum,sizeof(float));
	    pixelhisto(&output,scaleregionn,scalesort,0,&scalefactor,&regionmode,&sigma) ;
	    free(scalesort);
	    if(sigma<10)sigma=10.0; */
	  /*	output.image[i] = outy + sigma*gasdev(&seed);  */
	  output.image[i] = scalefactor +  sigma*gasdev(&seed);  
	  /*  free_matrix(pixval,1,10,1,10) ; */ 
				
	}
      }
    }


    if (flag_verbose) {
      sprintf(event,"Number of cosmic ray pixels masked= %d",numcrays);
      reportevt(flag_verbose,QA,1,event);
    }

    printf(" Finished cosmic ray masking \n");
  }
	
  /* ******************************************** */
  /* **********   WRITE OUTPUT IMAGE  *********** */
  /* ******************************************** */
  filemode=0;
  /* write the corrected image*/
  if (fits_write_img(output.fptr,TFLOAT,1,output.npixels,output.image,
		     &status)) {
    sprintf(event,"Writing image failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
	  
  /* Write information into the header describing the processing */
  /* get system time */
  tm=time(NULL);
  sprintf(comment,"%s",asctime(localtime(&tm)));
  comment[strlen(comment)-1]=0;
  if (fits_write_key_str(output.fptr,"DESMKMSK",comment, "Masked the image",&status)) {
    sprintf(event,"Writing processing history failed: %s", output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }
  sprintf(longcomment,"DESDM:");
  //  for (i=0;i<argc;i++) sprintf(longcomment,"%s %s",longcomment,argv[i]);
  sprintf(longcomment,"%s %s",longcomment,command_line);
  reportevt(flag_verbose,STATUS,1,longcomment);
  if (fits_write_comment(output.fptr,longcomment,&status)) {
    sprintf(event,"Writing longcomment failed: %s",output.name);
    reportevt(flag_verbose,STATUS,5,event);
    printerror(status);
  }

  fits_movrel_hdu(output.fptr, 1, NULL, &status);
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
  printf("end of program");      
        
  printf("\n");
  return(0);
}



#undef LINEAR_INTERP
#undef MOFFAT_INTERP
#undef GAUSS_INTERP



 
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
 

void getexposuretime(image,inputval,value,mode)
     desimage *image;
     int mode;
     char *inputval ;
     float *value;
{
  char currentimage[1000] ;
  static int status=0;
  int hdutype ;
  float val;
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

int main(int argc,char *argv[])
{
  return(MakeMask(argc,argv));
}

