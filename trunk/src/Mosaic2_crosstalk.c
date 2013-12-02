/*$Id$*/
/* Mosaic2_convert <> <> ....
   
   This program will transform Mosaic2 images from the Blanco 4m into 
   8 FITS files, each containing the imaging data from a single CCD. 
   Initial WCS parameters are simply copied directly from the original 
   image header Mosaic2 images with tnx information deleted and no attempt
   to load in an initial PV distortion model made.
   
   The programs alters the OBSTYPE and FILTER header keywords to conform to 
   DESDM standards, and leaves a history in the header.
  
*/
#include <string.h>
#include <time.h>
#include <math.h>
#include <getopt.h>

#include "imsupport.h"
#include "argutils.h"

#define INCLUDECOMMENTS 0
#define DONTINCLUDECOMMENTS 1

static const char *svn_id = "$Id$";


static int flag_sat = 0;
static int flag_overscan = 0;

void print_usage(char *program_name){
  printf("%s <infile.fits> <root/outfile> <options>\n",program_name);
  printf("  -crosstalk <crosstalk matrix- file>\n");
  printf("  -fixim <gain/rdnoise/saturate file>\n");
  printf("  -photflag <0 or 1>\n");
  printf("  -satmask\n");
  printf("  -overscansample <-1 for MEDIAN, 0 for MEAN, 1 for MEAN w/MINMAX>, \n");
  printf("  -overscanfunction < -50 < N < -1 for cubic SPLINE, 0 for LINE_BY_LINE, 1 < N < 50 for legendre polynomial>\n");
  printf("  -overscanorder <1-6>, default=1 <order for legendre polynomial>\n");
  printf("  -overscantrim <ncols>, default=0 <trim ncols at both edges of overscan, default = 0>\n");
  printf("  -overscan\n");
  printf("  -verbose <0-3>\n");
  printf("  -help (print usage and exit)\n");
  printf("  -version (print version and exit)\n");
}

int Mosaic2XTalk(int argc,char *argv[])
{
  char	filename[500],nfilename[500],nfilename_temp[500],event[500],newimagename[500],
    comment[1000],filetype[100],longcomment[10000],newname[500],
    filter[1000],obstype[1000],tr[300],*keycard,command_line[1000],
    ctype1[1000],ctype2[1000],trash[200],tag1[50],tag2[50],
    /* *********************************************** */
    /* **** define Detector sections for each CCD **** */
    /* *********************************************** */
    detsec[8][30]={"\'[1:2048,1:4096]\'","\'[2049:4096,1:4096]\'",
		   "\'[4097:6144,1:4096]\'","\'[6145:8192,1:4096]\'",
		   "\'[1:2048,4097:8192]\'","\'[2049:4096,4097:8192\']",
		   "\'[4097:6144,4097:8192]\'","\'[6145:8192,4097:8192]\'"},
    /*  ***********************************************  */
    /*  ***  List the Headers Keywords NOT to be ******  */
    /*  ***      copied from 0th Header          ******  */
    /*  ***********************************************  */
    *zeroheader,**exclkey,*nthheader,*nthp1header,*striparchiveroot(),
      exclist[100][10]={
      "SIMPLE","BITPIX","NAXIS","NAXIS1",
      "NAXIS2","EXTEND","NEXTEND","CHECKSUM",
      "DATASUM","CHECKVER","GCOUNT","PCOUNT",
      "BZERO","BSCALE","INHERIT","XTENSION",
      "TFIELDS","TFORM1","ZIMAGE","ZTILE1",
      "ZTILE2","ZCMPTYPE","ZNAME1","ZVAL1",
      "ZNAME2","ZVAL2","EXTNAME","ZTENSION",
      "ZBITPIX","ZNAXIS","ZNAXIS1","ZNAXIS2",
      "ZPCOUNT","ZGCOUNT","DCREATED","TTYPE1",
      "ZHECKSUM","RADECSYS","RADECEQ","RA    ",
      "DEC   ","TIMESYS","WAT0_001","WAT1_001",
      "WAT1_002", "WAT1_003", "WAT1_004", "WAT1_005",
      "WAT2_001", "WAT2_002", "WAT2_003", "WAT2_004",
      "WAT2_005", "ATM1_1","ATM2_2","ATV1",
      "ATV2","LTM1_1","LTM2_2","LTV1",
      "LTV2","DTM1_1","DTM2_2","DTV1",
      "DTV2","EXTVER","IMAGEID","PREFLASH",
      "WCSASTRM","BPM","CCDSEC"
    };
    int	nexc=71,numnthkeys,numwrittenkeys=0,keyflag,keynum,nkeys,
      totkeynum,numzerokeys,numnthp1keys,naxis,
      status=0,anynull,nfound,i,
      hdunum,x,y,hdutype,chdu,j,k,len,
      flag,wcsdim=2,chip_i,count,frameidlen,
      ccdnum,locin,locout,locin2,num_sata=0,num_satb=0,
      flag_crosstalk=0,flag_verbose=2,ext1,ext2,flag_fixim=0,bitpix,
      flag_phot=0,mkpath(),pixels,npixels,
      fpixel,oldpixel=0,getheader_flt(), getheader_str(), flag_osorder = 0;
    long	axes[2],naxes[2];
    float	gain_a,gain_b,rdnoise_a,rdnoise_b,saturate_a,saturate_b,
      significance,uncertainty,value,xoff[63],yoff[63],
      nullval,*outdata,**indata,
      /*  *************************************************** */
      /*  *** array for storing crosstalk coefficients ****** */
      /*  *************************************************** */
      xtalk[16][16],fix_gain[16],fix_rdnoise[16],fix_saturate[16];
    desimage input_image;
    desimage output_image;
    overscan_config osconfig;
    short *outmaskdata;
    float satA = 0.0;
    float satB = 0.0;
    double	xtalkcorrection=0.0;
    
    
    /*  ******************************************************** */
    /*  ******* other variables needed for the processing ****** */
    /*  ******************************************************** */
    void	printerror(),reportevt(),fixheader_flt(),
      fixheader_str(),init_desimage(),decodesection();
    fitsfile *fptr,*nfptr;
    FILE	*inp,*pip, *foffset;
    time_t	tm;
    enum { OPT_CROSSTALK=1,OPT_FIXIM,OPT_PHOTFLAG,OPT_SATMASK,OPT_OVERSCANSAMPLE,OPT_OVERSCANFUNCTION,
	   OPT_OVERSCANORDER,OPT_OVERSCANTRIM,OPT_OVERSCAN,OPT_VERBOSE,OPT_HELP,OPT_VERSION};
    
    osconfig.sample   = 0;
    osconfig.function = 0;
    osconfig.order    = 1;
    osconfig.trim     = 0;
    osconfig.debug    = 0;


    if(build_command_line(argc,argv,command_line,1000) <= 0){
      reportevt(2,STATUS,1,"Failed to record full command line.");
    }

    /* RAG: Added to print version of code to standard output (for logs) */
    sprintf(event,"%s",svn_id);
    reportevt(2,STATUS,1,event);
    reportevt(2,STATUS,1,command_line);

    if (argc<2) {
      print_usage(argv[0]);
      exit(0);
    }

    /* ****************************************************** */
    /* ****************  Process Command Line *************** */
    /* ****************************************************** */



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
    flag_sat = 0;
    flag_overscan = 0;
    while(1){
      int curind = optind;
      static struct option xtalk_options[] =
	{
	  {"crosstalk",        required_argument, 0,             OPT_CROSSTALK},
	  {"photflag",         required_argument, 0,             OPT_PHOTFLAG},
	  {"overscansample",   required_argument, 0,             OPT_OVERSCANSAMPLE},
	  {"overscanfunction", required_argument, 0,             OPT_OVERSCANFUNCTION},
	  {"overscanorder",    required_argument, 0,             OPT_OVERSCANORDER},
	  {"overscantrim",     required_argument, 0,             OPT_OVERSCANTRIM},
	  {"fixim",            required_argument, 0,             OPT_FIXIM},
	  {"verbose",          required_argument, 0,             OPT_VERBOSE},
	  {"version",          no_argument,       0,             OPT_VERSION},
	  {"help",             no_argument,       0,             OPT_HELP},
	  {"satmask",          no_argument,       &flag_sat,       1},
	  {"overscan",         no_argument,       &flag_overscan,  1},
	  {0,0,0,0}
	};
      int clopx = 0;

      clop = getopt_long_only(argc,argv,"",xtalk_options,&clopx);

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
	for (j=0;j<16;j++) 
	  for (k=0;k<16;k++) 
	    xtalk[j][k]=0.0;
	if(optarg){
	  /* read through crosstalk matrix file */
	  if (flag_verbose==3) {
	    sprintf(event,"Reading file %s",optarg);
	    reportevt(flag_verbose,STATUS,1,event);
	  }
	  inp=fopen(optarg,"r");
	  if (inp==NULL) {
	    sprintf(event,"Crosstalk file %s not found",optarg);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(0);
	  }
	  while (fgets(trash,200,inp)!=NULL) {
	    if (!strncmp(trash,"im",2)) {
	      /* first replace parentheses with spaces */
	      for (j=0;j<strlen(trash);j++) 
		if (!strncmp(&(trash[j]),"(",1) || 
		    !strncmp(&(trash[j]),")",1)) trash[j]=32; 
	      sscanf(trash,"%s %s %f %f %f",tag1,tag2,
		     &value,&uncertainty,&significance);
	      sscanf(&(tag1[2]),"%d",&ext1);
	      sscanf(&(tag2[2]),"%d",&ext2);
	      xtalk[ext1-1][ext2-1]=value;
	    }
	  }
	  if (fclose(inp)) {
	    sprintf(event,"File close failed: %s",optarg);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(0);
	  }
	  if (flag_verbose==3) {
	    for (j=0;j<16;j++) {
	      printf("  j=%2d  ",j);
	      for (k=0;k<16;k++) {
		printf("  %7.5f",xtalk[j][k]);
	      }
	      printf("\n");
	    }
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
      case OPT_FIXIM: // -fixim
	cloperr = 0;
	flag_fixim=1;
	j=0;
	if(optarg){
	  inp=fopen(optarg,"r");
	  if (inp==NULL) {
	    sprintf(event,"Fixim file %s not found",optarg);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	  while (fgets(trash,200,inp)!=NULL) {
	    if (strncmp(trash,"#",1)) {		
	      sscanf(trash,"%s %f %f %f",tr,fix_gain+j,
		     fix_rdnoise+j,fix_saturate+j);
	      j++;
	    }
	  }
	  if (fclose(inp)) {
	    sprintf(event,"File close failed: %s",optarg);
	    reportevt(flag_verbose,STATUS,5,event);
	    exit(1);
	  }
	}
	else cloperr = 1;
	if(cloperr){
	  reportevt(flag_verbose,STATUS,5,
		    "Option -fixim requires an argument.");
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
      sprintf(nfilename_temp,"%s",argv[optind]);
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
      sprintf(event,"File open failed: %s",filename);
      reportevt(flag_verbose,STATUS,3,event);
      /* Flavor=F (fpacked) */
      if (fits_open_file(&fptr,newname,READONLY,&status)) {
	sprintf(event,"File open failed: %s",newname);
	reportevt(flag_verbose,STATUS,3,event);
	status=0;
	sprintf(newname,"%s.gz",filename); 
	/* Flavor=G (gzipped) */
	if (fits_open_file(&fptr,newname,READONLY,&status)) {
	  sprintf(event,"File open failed: %s",newname);
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
    status=0;
    if (fits_get_num_hdus(fptr,&hdunum,&status)) {
      sprintf(event,"Reading HDUNUM failed: %s",filename);
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (flag_verbose) {
      sprintf(event,"%s has %d HDUs",filename,hdunum);
      reportevt(flag_verbose,STATUS,1,event);
    }
    if (hdunum!=17) {
      sprintf(event,"%d HDUs not standard Mosaic2 format",hdunum);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
	
    /* allocate vector to hold 16 images */
    indata=(float **)calloc(hdunum-1,sizeof(float *));
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
    status=0;
    if (fits_hdr2str(fptr,1,exclkey,nexc,&zeroheader,&numzerokeys,&status)) {
      sprintf(event,"Could not read kewords");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /* cut off last bogus "END" line */
    zeroheader[strlen(zeroheader)-80]=0;
    if (flag_verbose) {
      printf("  Read %d header keywords from 0th header\n",numzerokeys);
      if (flag_verbose==3) {
	printf("  *************************************************\n");
	printf("  %s",zeroheader);
	printf("*************************************************\n");
      }
    }

    /* read the OBSTYPE and use it to determine <filetype> */
    status=0;
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
      else if (!strcmp(obstype,"dome flat") || !strcmp(obstype,"DOME FLAT")
	       || !strcmp(obstype,"flat") || !strcmp(obstype,"FLAT"))
	sprintf(filetype,"raw_dflat");
      else if (!strcmp(obstype,"sky flat") || !strcmp(obstype,"SKY FLAT"))
	sprintf(filetype,"raw_tflat");	   
      else if (!strcmp(obstype,"dark") || !strcmp(obstype,"DARK"))
	sprintf(filetype,"raw_dark");
      else if (!strcmp(obstype,"focus") || !strcmp(obstype,"FOCUS"))
	sprintf(filetype,"raw_focus");
      else {
	sprintf(event,"OBSTYPE=%s not recognized",obstype);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
    }

    /* read the FILTER and convert to accepted type */
    status=0;
    if (fits_read_key_str(fptr,"FILTER",filter,comment,&status)) {
      sprintf(event,"Header Keyword FILTER not found");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    /*   ***********************************************  */
    /*   ***** Currently only set up for 11 filters **** */
    /*   ***********************************************  */
    if      (!strncmp(filter,"g",1)) sprintf(filter,"\'g\'");
    else if (!strncmp(filter,"r",1)) sprintf(filter,"\'r\'");
    else if (!strncmp(filter,"i",1)) sprintf(filter,"\'i\'");
    else if (!strncmp(filter,"z",1)) sprintf(filter,"\'z\'");
    else if (!strncmp(filter,"R",1)) sprintf(filter,"\'R\'");
    else if (!strncmp(filter,"B",1)) sprintf(filter,"\'B\'");
    else if(!strncmp(filter,"VR",2)) sprintf(filter,"\'VR\'");
    else if (!strncmp(filter,"V",1)) sprintf(filter,"\'V\'");
    else if (!strncmp(filter,"I",1)) sprintf(filter,"\'I\'");
    else if (!strncmp(filter,"Y",1)) sprintf(filter,"\'Y\'");
    else if (!strncmp(filter,"u",1)) sprintf(filter,"\'u\'");
    else if (!strncmp(filter,"U",1)) sprintf(filter,"\'U\'");
    else if (!strncmp(filter,"815",3)) sprintf(filter,"\'815\'");
    else if (!strncmp(filter,"823",3)) sprintf(filter,"\'823\'");
    else if(!strncmp(filter,"ha ",3)) sprintf(filter,"\'ha\'");
    else if(!strncmp(filter,"halp",4)) sprintf(filter,"\'ha8\'");
    else if(!strncmp(filter,"ha8",3)) sprintf(filter,"\'ha8\'");
    else if(!strncmp(filter,"H-alpha+8nm",11))sprintf(filter,"\'ha8\'");
    else if(!strncmp(filter,"SII",3)) sprintf(filter,"\'s2\'");
    else if(!strncmp(filter,"o3",2)) sprintf(filter,"\'o3\'");
    else if(!strncmp(filter,"o3off",5)) sprintf(filter,"\'Ooff\'");
    else if(!strncmp(filter,"C",1))sprintf(filter,"\'C\'");
    else if(!strncmp(filter,"D51",3))sprintf(filter,"\'D51\'");
    else if(!strncmp(filter,"M",1))sprintf(filter,"\'M\'");
    else if(!strncmp(filter,"DDS",3))sprintf(filter,"\'D51\'");
    else if(!strncmp(filter,"DDO",3))sprintf(filter,"\'D51\'");
    else {
      sprintf(event,"Unrecognized FILTER type %s",filter);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    /* apply fix to zeroheader */
    fixheader_str(zeroheader,"FILTER","FILTER",filter,flag_verbose);
    fixheader_str(zeroheader,"OBSTYPE","OBSTYPE",filetype,flag_verbose);
    /**********************************************************/
    /**********  Read All Image Data into Memory **************/
    /**********    Enables Crosstalk Correction  **************/
    /**********************************************************/

    if (flag_verbose) printf("  Reading image data from %s\n  ",
			     filename);

    for (i=0;i<hdunum-1;i++) {

      /* Move to the correct HDU */
      status=0;
      if (fits_movabs_hdu(fptr,i+2,&hdutype,&status)) {
	sprintf(event,"Move to HDU=%d failed: %s",i+2,filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      status=0;
      if (fits_get_img_param(fptr,2,&bitpix,&naxis,axes,&status)) {
	sprintf(event,"Image params not found in %s",filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /* double check image size */
      pixels  = axes[0]*axes[1];
      fpixel   = 1; nullval  = 0.0;
      if (i) {
	if (pixels!=oldpixel) {
	  sprintf(event,"Image extensions have different sizes:  %d  vs %d",
		  pixels,oldpixel);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
      }
      else oldpixel=pixels;

      indata[i]=(float *)calloc(pixels,sizeof(float));
      if (!indata[i]) {
	sprintf(event,"Error: Could not calloc indata[%d] %d",i,pixels);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }

      /* read the CHDU image  */
      status=0;
      if (fits_read_img(fptr,TFLOAT,fpixel,pixels,&nullval,
			indata[i],&anynull,&status)) {
	sprintf(event,"Reading image extension %d failed: %s",i,
		filename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      if (flag_verbose) {printf(".");fflush(stdout);}

    }
    if (flag_verbose) {printf("\n");fflush(stdout);}


    /**********************************************************/
    /**********     Cycle through extensions        ***********/
    /**********    preparing and writing data       ***********/
    /**********************************************************/

    naxes[0]=2148;naxes[1]=4096;
    npixels=naxes[0]*naxes[1];
    outdata=(float *)calloc(npixels,sizeof(float));
    if (!outdata) {
      sprintf(event,"Error: Could not calloc outdata %d",npixels);
      reportevt(flag_verbose,STATUS,5,event);
      exit(0);
    }
    /* Allocate and initialize the mask data */
    outmaskdata = NULL;
    if(flag_sat){
      outmaskdata=(short *)calloc(npixels,sizeof(short));
      if (!outmaskdata) {
	sprintf(event,"Error: Could not calloc outmaskdata %d",npixels);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
    }
    /* jump to 2nd HDU  */
    status=0;
    if (fits_movabs_hdu(fptr,2,&hdutype,&status)) {
      sprintf(event,"Move to HDU=2 failed");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);  
    }

    /* now cycle through the HDU's parsing properly and writing */
    /* useful FITS file */
    for (i=2;i<=hdunum;i+=2) {
      /* remember that for i=2 the image number is 0 */
      k=i-2;
      
      /* 
	 Use DESIMAGE data structure to store information about
	 the input image.  This sets input_image image buffers
	 for all of the processing done in this utility. 
      */
      init_desimage(&input_image);
      init_desimage(&output_image);
      input_image.image = outdata;
      input_image.mask  = outmaskdata;
      input_image.npixels = npixels;
      input_image.axes[0] = naxes[0];
      input_image.axes[1] = naxes[1];
      fits_get_hdu_num(fptr,&chdu);
      if (flag_verbose) printf("  Currently located at HDU %d of %d\n",
			       chdu,hdunum);
      if (chdu!=i) {
	sprintf(event,"Not located at correct HDU (%d instead of %d)",
		chdu,i);
	reportevt(flag_verbose,STATUS,5,event);
	exit(0);
      }
      
      
      /* ****************************************************** */
      /* ***********  Read Header from Nth Extension ********** */
      /* ****************************************************** */
      status=0;
      if (fits_hdr2str(fptr,INCLUDECOMMENTS,exclkey,nexc,&nthheader,
		       &numnthkeys,&status)) {
	sprintf(event,"Header %d not read",i-1);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* cut off last bogus "END" line */
      nthheader[strlen(nthheader)-80]=0;
      if (flag_verbose) {
	printf("  Read %d header keywords\n",numnthkeys);
	if (flag_verbose==3) {
	  printf("  *************************************************\n");
	  printf("  %s",nthheader);
	  printf("*************************************************\n");
	}
      }
      status = 0;
      /* Get a bunch of useful info from the AMP(A)-specific header */
      if(getheader_flt(nthheader,"SATURATE",&input_image.saturateA,flag_verbose)){
	sprintf(event,"Header Keyword for saturate not found for amp A.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_flt(nthheader,"GAIN",&input_image.gainA,flag_verbose)){
	sprintf(event,"Header Keyword for gain not found for amp A.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_flt(nthheader,"RDNOISE",&input_image.rdnoiseA,flag_verbose)){
	sprintf(event,"Header Keyword for rdnoise not found for amp A.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      sprintf(input_image.datasec,"'[51:2098,1:4096]'");
      sprintf(input_image.ampseca,"'[1:1024,1:4096]'");
      sprintf(input_image.biasseca,"'[1:50,1:4096]'");
      sprintf(input_image.trimsec,"'[51:2098,1:4096]'");
      decodesection((input_image.biasseca), (input_image.biassecan),
		    flag_verbose);
      decodesection((input_image.datasec), (input_image.datasecn),flag_verbose);
      decodesection((input_image.ampseca),  (input_image.ampsecan), 
		    flag_verbose);
      decodesection((input_image.trimsec), (input_image.trimsecn),flag_verbose);
      
      
      /* ****************************************************** */
      /* **************  Move to Nth+1 Extension ************** */
      /* ****************************************************** */
      
      status=0;
      
      if (fits_movrel_hdu(fptr,1,&hdutype,&status)) {
	sprintf(event,"Move to next HDU failed");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      
      /* ****************************************************** */
      /* **********  Read Header from Nth+1 Extension ********* */
      /* ****************************************************** */
      
      status=0;
      if (fits_hdr2str(fptr,INCLUDECOMMENTS,exclkey,nexc,&nthp1header,
		       &numnthp1keys,&status)) {
	sprintf(event,"Header %d not read",i);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      /* cut off last bogus "END" line */
      nthp1header[strlen(nthp1header)-80]=0;
      if (flag_verbose) {
	printf("  Read %d header keywords\n",numnthp1keys);
	if (flag_verbose==3) {
	  printf("  *************************************************\n");
	  printf("  %s",nthp1header);
	  printf("*************************************************\n");
	}
      }
      status = 0;
      /* Get a bunch of useful info from the AMP(B)-specific header */
      if(getheader_flt(nthheader,"SATURATE",&input_image.saturateB,flag_verbose)){
	sprintf(event,"Header Keyword for saturate not found for amp B.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_flt(nthheader,"GAIN",&input_image.gainB,flag_verbose)){
	sprintf(event,"Header Keyword for gain not found for amp B.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      if(getheader_flt(nthheader,"RDNOISE",&input_image.rdnoiseB,
		       flag_verbose)){
	sprintf(event,"Header Keyword for rdnoise not found for amp B.");
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);	  
      }
      sprintf(input_image.ampsecb,"'[1025:2048,1:4096]'");
      sprintf(input_image.biassecb,"'[2099:2148,1:4096]'");
      decodesection((input_image.biassecb), (input_image.biassecbn),
		    flag_verbose);
      decodesection((input_image.ampsecb),  (input_image.ampsecbn), 
		    flag_verbose);


      //      gain_a=gain_b=rdnoise_a=rdnoise_b=saturate_a=saturate_b=0.0;
      /* alter these values if fixim flag is set */
      if (flag_fixim) {
	input_image.gainA    = fix_gain[i-2];
	input_image.gainB    = fix_gain[i-1];
	input_image.rdnoiseA = fix_rdnoise[i-2];
	input_image.rdnoiseB = fix_rdnoise[i-1];
	input_image.saturateA = fix_saturate[i-2];
	input_image.saturateB = fix_saturate[i-1];
      }
      if(flag_verbose){
	sprintf(event,"GainA: %.2f RDNoiseA %.2f SATURATA: %.0f  GainB: %.2f RDNoiseB %.2f SATURATB: %.0f",
		input_image.gainA,input_image.rdnoiseA,input_image.saturateA,
		input_image.gainB,input_image.rdnoiseB,input_image.saturateB); 
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"AMPSECA: %s BIASSECA: %s AMPSECB: %s BIASSECB: %s",
		input_image.ampseca,input_image.biasseca,input_image.ampsecb,
		input_image.biassecb); 
	reportevt(flag_verbose,STATUS,1,event);
	sprintf(event,"DATASEC: %s TRIMSEC: %s",
		input_image.datasec,input_image.trimsec); 
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* ************************************************************* */
      /* ************** copy amp images into ccd images ************** */
      /* ************** apply the crosstalk correction *************** */
      /* ************************************************************* */
      
      /* Reset the mask if enabled */
      if(flag_sat){
	for(y = 0;y < npixels;y++)
	  outmaskdata[y] = 0;
      }
      num_sata = 0;
      num_satb = 0;
      /* copy the amp images into the new ccd image */
      for (y=0;y<naxes[1];y++) for (x=0;x<naxes[0];x++) {
	  locout=y*naxes[0]+x;
	  /* copy overscan region for first amplifier */
	  if (x<50) {
	    locin=y*axes[0]+x+1062;
	    locin2=y*axes[0]+(axes[0]-(x+1062)-1);
	    outdata[locout]=indata[k][locin];
	    /* detect saturated pixels if enabled */
	    if(flag_sat && (outdata[locout] > input_image.saturateA)){
	      outmaskdata[locout] = BADPIX_SATURATE;
	      num_sata++;
	    }
	    if (flag_crosstalk)   {/* apply crosstalk correction */
	      for (j=0;j<16;j+=2) outdata[locout]-=
				    (xtalk[k][j]*indata[j][locin]+xtalk[k][j+1]*
				     indata[j+1][locin2]);
	    }
	  }
	  /* copy data from first amplifier */
	  else if (x<1074) {
	    locin=y*axes[0]+(x-50)+24;
	    locin2=y*axes[0]+(axes[0]-((x-50)+24)-1);
	    outdata[locout]=indata[k][locin];
	    if(flag_sat && (outdata[locout] > input_image.saturateA)){
	      outmaskdata[locout] = BADPIX_SATURATE;
	      num_sata++;
	    }
	    if (flag_crosstalk)   /* apply crosstalk correction */
	      for (j=0;j<16;j+=2) outdata[locout]-=
				    (xtalk[k][j]*indata[j][locin]+xtalk[k][j+1]*
				     indata[j+1][locin2]);
	  }
	  /* copy data from second amplifier */
	  else if (x<2098) {
	    locin=y*axes[0]+(x-1074)+64;
	    locin2=y*axes[0]+(axes[0]-((x-1074)+64)-1);
	    outdata[locout]=indata[k+1][locin];
	    if(flag_sat && (outdata[locout] > input_image.saturateB)){
	      outmaskdata[locout] = BADPIX_SATURATE;
	      num_satb++;
	    }
	    if (flag_crosstalk)   /* apply crosstalk correction */
	      for (j=0;j<16;j+=2) outdata[locout]-=
				    (xtalk[k+1][j]*indata[j][locin2]+xtalk[k+1][j+1]*
				     indata[j+1][locin]);
	  }
	  /* copy overscan region for second amplifier */
	  else {
	    locin=y*axes[0]+(x-2098);
	    locin2=y*axes[0]+(axes[0]-(x-2098)-1);
	    outdata[locout]=indata[k+1][locin];
	    if(flag_sat && (outdata[locout] > input_image.saturateB)){
	      outmaskdata[locout] = BADPIX_SATURATE;
	      num_satb++;
	    }
	    if (flag_crosstalk)   /* apply crosstalk correction */
	      for (j=0;j<16;j+=2) outdata[locout]-=
				    (xtalk[k+1][j]*indata[j][locin2]+xtalk[k+1][j+1]*
				     indata[j+1][locin]);
	  }
	}
      
      if(flag_verbose && flag_sat){
	sprintf(event,"image=%s,CCD[%d] NSATPIX=%d",
		filename,i-1,num_sata+num_satb);
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* PERFORM OVERSCAN if not disabled */
      if(flag_overscan){
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
	
	/* update dataset within the output image */
	sprintf(output_image.datasec,"[%d:%ld,%d:%ld]",1,output_image.axes[0],
		1,output_image.axes[1]);
	if(flag_verbose){
	  sprintf(event,"OVERSCAN: New DataSec = %s",output_image.datasec);
	  reportevt(flag_verbose,STATUS,1,event);	  
	}
      }
      else { 
	// If overscan was not used, then just straight copy the input image
	output_image.image      = input_image.image;
	output_image.mask       = input_image.mask;
	output_image.axes[0]    = input_image.axes[0];
	output_image.axes[1]    = input_image.axes[1];
	output_image.gainA      = input_image.gainA;
	output_image.gainB      = input_image.gainB;
	output_image.rdnoiseA   = input_image.rdnoiseA;
	output_image.rdnoiseB   = input_image.rdnoiseB;
	output_image.saturateA  = input_image.saturateA;
	output_image.saturateB  = input_image.saturateB;
	output_image.npixels    = input_image.npixels;  
	strncpy(output_image.datasec,input_image.datasec,100);
      }
      strncpy(output_image.trimsec,input_image.trimsec,100);
      strncpy(output_image.biasseca,input_image.biasseca,100);
      strncpy(output_image.biassecb,input_image.biassecb,100);
      strncpy(output_image.ampseca,input_image.ampseca,100);
      strncpy(output_image.ampsecb,input_image.ampsecb,100);
	

      /* ****************************************************** */
      /* *******       Alter Selected Keywords           ****** */
      /* ****************************************************** */
      fixheader_flt(nthheader,"GAIN","GAINA",output_image.gainA,flag_verbose);
      fixheader_flt(nthp1header,"GAIN","GAINB",output_image.gainB,flag_verbose);
      fixheader_flt(nthheader,"RDNOISE","RDNOISEA",output_image.rdnoiseA,
		    flag_verbose);
      fixheader_flt(nthp1header,"RDNOISE","RDNOISEB",output_image.rdnoiseB,
		    flag_verbose);
      fixheader_flt(nthheader,"SATURATE","SATURATA",output_image.saturateA,
		    flag_verbose);
      fixheader_flt(nthp1header,"SATURATE","SATURATB",output_image.saturateB,
		    flag_verbose);
      fixheader_str(nthheader,"AMPNAME","AMPNAMEA","",flag_verbose);
      fixheader_str(nthp1header,"AMPNAME","AMPNAMEB","",flag_verbose);
      fixheader_str(nthheader,"AMPSEC","AMPSECA",output_image.ampseca,
		    flag_verbose);
      fixheader_str(nthp1header,"AMPSEC","AMPSECB",output_image.ampsecb,
		    flag_verbose);
      fixheader_str(nthheader,"BIASSEC","BIASSECA",output_image.biasseca,
		    flag_verbose);
      fixheader_str(nthp1header,"BIASSEC","BIASSECB",output_image.biassecb,
		    flag_verbose);
      fixheader_str(nthheader,"DATASEC","DATASEC",output_image.datasec,
		    flag_verbose);
      fixheader_str(nthheader,"TRIMSEC","TRIMSEC",output_image.trimsec,
		    flag_verbose);
      fixheader_str(nthheader,"DETSEC","DETSEC",detsec[i/2-1],
		    flag_verbose);
      fixheader_str(nthheader,"FILTER","FILTER",filter,flag_verbose);
      fixheader_str(nthp1header,"FILTER","FILTER",filter,flag_verbose);
      if (!strcmp(filetype,"raw_obj")) { /* fix the WCS params */
	fixheader_str(nthheader,"CTYPE1","CTYPE1","\'RA---TAN\'",
		      flag_verbose);
	fixheader_str(nthheader,"CTYPE2","CTYPE2","\'DEC--TAN\'",
		      flag_verbose);
      }
      fixheader_str(nthheader,"OBSTYPE","OBSTYPE",filetype,flag_verbose);
      fixheader_str(nthp1header,"OBSTYPE","OBSTYPE",filetype,flag_verbose);
      if (flag_verbose==3) {
	printf("  *************************************************\n");
	printf("  %s",nthheader);
	printf("*************************************************\n");
	printf("  *************************************************\n");
	printf("  %s",nthp1header);
	printf("*************************************************\n");
      }
      
      /* **************************************************************** */
      /* **************** Write New Single CCD Image ******************** */
      /* ************ Store Processing History in Header **************** */
      /* **************************************************************** */
      
      /* make sure path exists first time through */
      if (i==2) {
	sprintf(nfilename,"%s_%02d.fits",nfilename_temp,i/2);
	if (mkpath(nfilename,flag_verbose)) {
	  sprintf(event,"Failed to create path to file: %s",nfilename);
	  reportevt(flag_verbose,STATUS,5,event);
	  exit(0);
	}
	else {
	  sprintf(event,"Created path to file: %s",nfilename);
	  reportevt(flag_verbose,STATUS,1,event);
	}
      }
      sprintf(nfilename,"!%s_%02d.fits",nfilename_temp,i/2);
      status=0;
      if (fits_create_file(&output_image.fptr,nfilename,&status)) {
	sprintf(event,"File creation of %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) {
	sprintf(event,"Opened image for CCD %d : %s ",i/2,nfilename+1);
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* create image extension  */
      status=0;
      if (fits_create_img(output_image.fptr,FLOAT_IMG,2,output_image.axes,&status)) {
	sprintf(event,"Creating image failed: %s",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /* first reserve room within the header for the known # of keywords */
      status=0;
      if (fits_set_hdrsize(output_image.fptr,numzerokeys+numnthkeys+15,&status)) {
	sprintf(event,"Reserving header space for %d keywords failed in %s",
		numzerokeys+numnthkeys+10,nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }


      /* copy zero header information into the new header */
      /* note that last keyword is null keyword and truncates the header! */
      /* only copy the header params that live within zerokeyword */
      status=0;
      if (fits_get_hdrpos(output_image.fptr,&totkeynum,&keynum,&status)) {
	sprintf(event,"Reading header position in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      numwrittenkeys=0;
      for (j=0;j<numzerokeys-1;j++) {
	/* copy the record into the header */
	status=0;
	keycard = zeroheader+j*80;
	if( strncmp(keycard,"COMMENT   FITS (Flexible Image",30) &&
	    strncmp(keycard,"COMMENT   and Astrophysics', v",30))
	  {

	    if (fits_insert_record(output_image.fptr,totkeynum+numwrittenkeys,
				   keycard,&status)) {
	      sprintf(event,"Zero Header insert in %s failed",nfilename);
	      reportevt(flag_verbose,STATUS,5,event);
	      printerror(status);
	    }
	    else
	      numwrittenkeys++;
	  }
      }
      if (flag_verbose) {
	sprintf(event,"  Copied %d keywords from 0th header\n",numwrittenkeys);
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* copy the header information from the Nth HDU of MEF */
      status=0;
      if (fits_get_hdrpos(output_image.fptr,&totkeynum,&keynum,&status)) {
	sprintf(event,"Reading header position in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose==3) {
	sprintf(event,"Currently at header position %d with %d keywords",
		keynum,totkeynum);
	reportevt(flag_verbose,STATUS,1,event);
      }
      numwrittenkeys=0;
      for (j=0;j<numnthkeys-1;j++) {
	/* don't write keyword if it is duplicate */
	keyflag=0;
	for (k=0;k<numzerokeys;k++) {
	  if (!strncmp(zeroheader+k*80,nthheader+j*80,8)) {
	    keyflag=1;
	    break;
	  }
	}
	if (!keyflag) {
	  /* copy the record into the header */
	  status=0;
	  keycard = nthheader+j*80;
	  if( strncmp(keycard,"COMMENT   FITS (Flexible Image",30) &&
	      strncmp(keycard,"COMMENT   and Astrophysics', v",30)){
	    if (fits_insert_record(output_image.fptr,totkeynum+numwrittenkeys,
				 keycard,&status)) {
	      sprintf(event,"Insert of %dth header keyword in %s failed:\n%s",
		      numwrittenkeys,nfilename,nthheader+j*80);
	      reportevt(flag_verbose,STATUS,5,event);
	      printerror(status);
	    }
	    else
	      numwrittenkeys++;
	  }
	}
      }
      if (flag_verbose){
	sprintf(event,"  Copied %d keywords from Nth header\n",
		numwrittenkeys);
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* copy the header information from the Nth+1 HDU of MEF */
      status=0;
      if (fits_get_hdrpos(output_image.fptr,&totkeynum,&keynum,&status)) {
	sprintf(event,"Reading header position in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose==3) 
	{
	  sprintf(event,"Currently at header position %d with %d keywords",
		  keynum,totkeynum);
	  reportevt(flag_verbose,STATUS,1,event);
	}
      numwrittenkeys=0;
      for (j=0;j<numnthp1keys-1;j++) 
	{
	  /* don't write keyword if it is duplicate */
	  keyflag=0;
	  for (k=0;k<numzerokeys;k++) 
	    {
	      if (!strncmp(zeroheader+k*80,nthp1header+j*80,8)) 
		{
		  keyflag=1;
		  break;
		}
	    }
	  for (k=0;k<numnthkeys;k++) 
	    {
	      if (!strncmp(nthheader+k*80,nthp1header+j*80,8)) 
		{
		  keyflag=1;
		  break;
		}
	    }
	  if (!keyflag) 
	    {
	      /* copy the record into the header */
	      status=0;
	      keycard = nthp1header+j*80;
	      if( strncmp(keycard,"COMMENT   FITS (Flexible Image",30) &&
		  strncmp(keycard,"COMMENT   and Astrophysics', v",30)){
		if (fits_insert_record(output_image.fptr,
				       totkeynum+numwrittenkeys,
				       keycard,&status)) {
		  sprintf(event,
			  "Insert of %dth header keyword in %s failed:\n%s",
			  numwrittenkeys,nfilename,keycard);
		  reportevt(flag_verbose,STATUS,5,event);
		  printerror(status);
		}
		else
		  numwrittenkeys++;
	      }
	    }
	}
      if (flag_verbose) {
	sprintf(event,"  Copied %d keywords from Nth+1 header\n",
		numwrittenkeys);
	reportevt(flag_verbose,STATUS,1,event);
      }
      
      /* write the new image */
      status=0;
      if (fits_write_img(output_image.fptr,TFLOAT,1,output_image.npixels,
			 output_image.image,&status))
	printerror(status);
      if (flag_verbose) printf("  Wrote image data to %s\n",nfilename+1);
      
      
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

      /* add PHOTFLAG */
      status=0;
      if (fits_update_key_lng(output_image.fptr,"PHOTFLAG",flag_phot,
			      "Night Photometric (1) or not (0)",&status)) {
	sprintf(event,"Keyword PHOTFLAG insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) printf("PHOTFLAG = %d\n",flag_phot);

      /* update OBSTYPE keyword */
      status=0;
      if (fits_update_key_str(output_image.fptr,"OBSTYPE",filetype,
			      "Type of observation",&status)) {
	sprintf(event,"Keyword OBSTYPE update in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) printf("  OBSTYPE = %s\n",filetype);

      /* add CDDNUM keyword */
      ccdnum=i/2;
      status=0;
      if (fits_insert_key_lng(output_image.fptr,"CCDNUM",ccdnum,"CCD number (1-8)",
			      &status)) {
	sprintf(event,"Keyword CCDNUM insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      /* ************************************************** */
      /* ***** leave processing history in the header ***** */
      /* ************************************************** */
      tm=time(NULL);
      sprintf(comment,"%s",asctime(localtime(&tm)));
      comment[strlen(comment)-1]=0;
      status=0;
      if (fits_write_key_str(output_image.fptr,"DESMSxtk",comment,
			     "Mosaic2 image conversion and crosstalk correction",&status)) {
	sprintf(event,"Keyword DESMSxtk insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if(flag_overscan){
	if (fits_write_key_str(output_image.fptr,"DESOSCN",comment, 
			       "overscan corrected",&status)) {
	  sprintf(event, 
		  "DESDM: Setting DESOSCN=overscan corrected failed: %s", 
		  nfilename);
	  reportevt(flag_verbose, STATUS, 5, event);
	  printerror(status);
	}
	sprintf(longcomment, "DESDM: (%d) ", 1);
	sprintf(longcomment, "%s %s ", longcomment, 
		striparchiveroot(filename));
	if (flag_verbose) {
	  sprintf(event, "%s", longcomment);
	  reportevt(flag_verbose, STATUS, 1, event);
	}
	if (fits_write_history(output_image.fptr, longcomment, &status)) {
	  sprintf(event, 
		  "DESDM: Writing overscan comment in output image failed: %s",
		  nfilename);
	  reportevt(flag_verbose, STATUS, 5, event);
	  printerror(status);
	}
      }
      if(flag_sat){ 
	if (fits_write_key_lng(output_image.fptr,"NSATPIX",(long)(num_sata+num_satb),
			       "Number of saturated pixels",&status)) {
	  sprintf(event,"Setting number of saturated pixels failed: %ld",(long)(num_sata+num_satb));
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
      sprintf(longcomment,"DESDM:");
      //      for (j=0;j<argc;j++) sprintf(longcomment,"%s %s",longcomment,argv[j]);
      sprintf(longcomment,"%s %s",longcomment,command_line);
      status=0;
      if (fits_write_history(output_image.fptr,longcomment,&status)) {
	sprintf(event,"DESDM: longcomment insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose==3) printf("  => %s\n",longcomment);
      status=0;
      if (fits_write_key_str(output_image.fptr,"DES_EXT","IMAGE",
			     "Image extension",&status)) {
	sprintf(event,"Keyword DES_EXT insert in %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }

      if(flag_sat){
	/* Create extension for the mask */
	if (fits_create_img(output_image.fptr,USHORT_IMG,2,
			    output_image.axes,&status)) {
	  sprintf(event,"Creating image mask failed.");
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
	/* Write the mask data */	  
	if (fits_write_img(output_image.fptr,TUSHORT,1,output_image.npixels,
			   output_image.mask,
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
	/* May be needed - maybe not ... duped from imcorrect mask dump */
	if (fits_update_key_lng(output_image.fptr,"EXTVER",2,
				"Extension version",&status)) {
	  reportevt(flag_verbose,STATUS,5,"Setting EXTVER=2 failed");
	  printerror(status); 
	}
      }
 
      /* close the new image */
      status=0;
      if (fits_close_file(output_image.fptr,&status)) {
	sprintf(event,"Closing %s failed",nfilename);
	reportevt(flag_verbose,STATUS,5,event);
	printerror(status);
      }
      if (flag_verbose) {
	sprintf(event,"Closed image %ld X %ld OBSTYPE = %s PHOTFLAG = %d : %s",
		output_image.axes[0],output_image.axes[1],filetype,
		flag_phot,nfilename+1);
	reportevt(flag_verbose,STATUS,1,event);
      }

      /* Move to the next HDU */
	  
      if (i<hdunum-1) {
	status=0;
	if (fits_movrel_hdu(fptr,1,&hdutype,&status)) {
	  sprintf(event,"Move to next HDU failed");
	  reportevt(flag_verbose,STATUS,5,event);
	  printerror(status);
	}
      }
    }
	

    /* ********************************************************* */
    /* ************** Close Input Image and Exit *************** */
    /* ********************************************************* */

    status=0;
    if (fits_close_file(fptr,&status)) {
      sprintf(event,"Closing Exposure failed");
      reportevt(flag_verbose,STATUS,5,event);
      printerror(status);
    }
    if (flag_verbose) {
      sprintf(event,"Closed image %s",filename);
      reportevt(flag_verbose,STATUS,1,event);
    }
    return(0);
}

void fixheader_str(char *hdr,char *oldkey,char *newkey,
		   char *value,int flag_verbose)
{
	
  int	len,lines,flag=0,keylen,i,nlinelen,j,cflag=0,nkeylen;
  char	event[500],nline[100],comment[80];
  void	reportevt();

  len=strlen(hdr);
  keylen=strlen(oldkey);
  nkeylen=strlen(newkey);if (nkeylen>8) nkeylen=8;
  lines=len/80;

  for (i=0;i<lines;i++) {
    if (!strncmp(oldkey,hdr+i*80,keylen)) {
      flag=1;
      if (value=="") strncpy(hdr+i*80,newkey,nkeylen);
      else {
	sprintf(nline,"%-8s= %-20s",newkey,value);
	nlinelen=strlen(nline);  if (nlinelen>79) nlinelen=79;
	strncpy(hdr+i*80,nline,nlinelen);
      }
      break;
    }
  }
  if (!flag && flag_verbose) {
    sprintf(event,"Keyword %s not found and updated to %s",oldkey,newkey);
    reportevt(flag_verbose,STATUS,5,event);
    exit(0); 
  }
}

void fixheader_flt(hdr,oldkey,newkey,value,flag_verbose)
  char *hdr,*oldkey,*newkey;
  float value;
  int flag_verbose;
{
  int	len,lines,flag=0,keylen,i,nlinelen,j,cflag=0;
  char	event[500],nline[100],comment[80];
  void	reportevt();
  len=strlen(hdr);
  keylen=strlen(oldkey);
  lines=len/80;
  for (i=0;i<lines;i++) {
    if (!strncmp(oldkey,hdr+i*80,keylen)) {
      flag=1;
      /* just change the keyword */
      if (fabs(value)<1.0e-6)strncpy(hdr+i*80,newkey,strlen(newkey));
      else { /* change keyword and value while preserving comment */
	//	if (!cflag) sprintf(comment,"");
	sprintf(nline,"%-8s= '%18.2f'",newkey,value);
	nlinelen=strlen(nline)-1;  if (nlinelen>79) nlinelen=79;
	strncpy(hdr+i*80,nline,strlen(nline)-1);
      }
      break;
    }
  }

  if (!flag && flag_verbose) {
    sprintf(event,"Keyword %s not found and updated to %s",oldkey,newkey);
    reportevt(flag_verbose,STATUS,5,event); 
    /*printf("******\n%s\n*******\n",hdr); */
    exit(0);  
  }

  else if (flag_verbose==3) {
    sprintf(event,"Keyword %s found and updated to %s",oldkey,newkey);
    reportevt(flag_verbose,STATUS,1,event);
  }

}

int getheader_flt(char *hdr,char *key,float *value,int flag_verbose)
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
      *value = atof(answer);
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
  else if (flag_verbose==3) {
    sprintf(event,"Keyword %s found to be %.2f",key,*value);
    reportevt(flag_verbose,STATUS,1,event);
  }
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
  else if (flag_verbose==3) {
    sprintf(event,"Keyword %s found to be %s",key,value);
    reportevt(flag_verbose,STATUS,1,event);
  }
  return(0);
}

int main(int argc,char *argv[])
{
  return(Mosaic2XTalk(argc,argv));
}
