///
/// \file MakeBleedMask.cc
/// \brief Command line interface for bleedtrail masking
/// \author Mike Campbell (mtcampbe@illinois.edu)
///
/// Performs bleedtrail processing on input images, marks
/// BADPIX_TRAIL for each bleedtrail pixel
/// 
/*--MakeBleedMask

Basic syntax:
mkbleedmask [-aghimz] [-b <factor> -f <value> -l <value> -n <npix> 
            -o <filename> -r <factor> -s <n> -t <npixels> 
            -v <level> -w <value> -x <filename> ] <infile> <outfile> 

Input data
  <infile>   = FITS image (hopefully with astrometric soln)
  <outfile>  = FITS image output 
Options:
  Background Determination:
   -g,--global_only:       Do not use local statistics.
   -s,--bgiters <n>:       Number of iterations for background estimation. 
                           (default=1)
   -b,--bgreject <factor>: Use specified <factor> as scalefactor for
                           background rejection. (default=5.0)
   -r,--growbox <factor>:  Factor by which to grow blob boxes for determining
                           a background level. (defualt=2)

  Filter Description Options:
   -t,--trail_length <pix>: Use specified <npixels> for trail structuring 
                            element in Y. (default=20)
   -n,--numtrail <pix>:     Number of contiguous pixels required for trail
                             detection. (default=trail_length/2)
   -f,--scalefactor <value>:Use specified <value> as scalefactor on background 
                            to detect bleedtrails. (default=10.0)

  Bleed Trail Options:
  -i,--interpolate:         Interpolate over bleedtrails.

  Saturated Object Options:
  -l,--starlevel <value>:   Use specified <value> as scalefactor on background 
                            to detect stars. (default=5.0)
  -w,--growrad <value>:     Factor by which to grow detected star radii. 
                            (default=1.0)
  -a,--interpolate-stars:   Do radial interpolation over saturated stars for 
                            which bleedtrails were detected (default=No)
  -m,--starmask:            Create a mask for the detected bright objects.
                            (default=No)
  -z,--zerostarweights:     Set weights to zero for masked stars. (default=No)

  Miscellaneous Output Options:
  -o,--saturated_objects <filename>: Filename for saturated star table.(None)
  -x,--trailboxes <filename>: Filename for trail box table.(None)

  -h,--help:                Prints long version of help.
  -v,--verb <level>:        Verbosity level [0-3], default = 1.

Summary:
  Detects bleedtrails in Y for incoming FITS image.
  Sets BADPIX_TRAIL for detected bleedtrails,
  and BADPIX_INTERP if interpolation is enabled,
  and BADPIX_STAR if star masking is enabled.
  0's out weight pixels for detected trails.
  Creates object table with WCS positions and radii.

Detailed Description:
  This code searches for and masks/interpolates over bleed trails by using a 
  linear filter that searches for values above the local background.  In 
  addition saturated objects are identified and masked.  The code is highly 
  tunable with a number of options which are described below.

  Background estimation options:
    The determination of the background level and RMS is made locally unless 
    the -g option is specified.  The size of the local region considered to
    estimate the background is set by the -r option.  When determining the 
    background and RMS levels the program uses a sigma-clipping algorithm 
    which rejects values that are more than N*sigma from the median (where 
    N is set by the -b option).  The maximum number of iterations that may be 
    used is set by the -s option.  

  Filter Description Options:
    The geometry of the filter used to detect bleed trails is a single pixel
    wide kernel with length set by the -t option.  The -n option determines
    how many contiguous pixels within the filter must be above the background
    by a factor of M times the RMS (where the -f option sets the value of M).
     
  Bleed Trail Options:
    The default behavior is to flag the pixels within the bleed trail but 
    if the -i option is given then pixels identified as part of the trail
    are also interpolated (using values of pixels adjacent to the trail).

  Saturated Object Options:
    In addition to flagging/interpolating bleed trails, saturated objects can
    be masked (by setting the -m option).  Saturated objects are identified 
    and then a radial profile is used to estimate the size of the object by
    searching for the point where the surrounding pixels are on average more
    than L times the RMS above the background (where L is set by the -l option).
    If the -w option is given then masked radius is further multiplied by
    a constant factor to increase the size of the region masked.  In addition
    to masking, when the -z option is used the weight image is altered to a 
    value of 0.0 over the masked areas.  Finally, if the -a option is used
    then the region around the star that is also part of the bleed trail is
    interpolated using the radial profile (Note: this last option can remove
    sources close to saturated objects and therefore is not recommended for 

  Miscellaneous Output Options:
    Two miscellaneous FITS tables can be output.  The -x options causes a FITS
    table to be written which describes the region (corners) affected by each 
    bleed trail.  Similarly, the -o options causes a similar output for 
    saturated object (giving the center and radius of the area masked around
    each object).  If a WSC solution is present and SCAMPFLG=0 then the values
    in these tables are given in decimal degrees and arcseconds.  Otherwise,
    pixel units are used in the tables.

Known Issues:
  Better handling of "edge bleed" (bleed trails that saturate to the chip edge and then 
  further affect the portion of the image near the amplifier.

*/

#include "ComLine.hh"
#include "FitsTools.hh"
#include "BleedTrails.hh"
//#include "Profiler.hh"

extern "C" {
#undef   OFF_T
#include "fitscat_mine.h"
#include "fitswcs.h"

static const char *svn_id = "$Id$";

}

///
/// @brief ComLineObject for mkbleedmask app
///
/// ComLine is a convenience object to help
/// deal with command-line arguments and
/// options.
class MakeBleedMaskComLine : public ComLineObject
{
public:
  MakeBleedMaskComLine()
    : ComLineObject()
  { Initialize(); };
  MakeBleedMaskComLine(const char *args[])
    : ComLineObject(args)
  { Initialize(); };
  void Initialize(){
    AddOption('a',"interpolate-stars");
    AddOption('b',"bgreject",2,"factor");
    AddOption('d',"edgesize",2,"npix");
    AddOption('e',"version");
    AddOption('f',"scalefactor",2,"value");
    AddOption('g',"global_only");
    AddOption('h',"help");
    AddOption('i',"interpolate");
    AddOption('j',"trailreject",2,"npix");
    AddOption('l',"starlevel",2,"value");
    AddOption('m',"starmask");
    AddOption('n',"numtrail",2,"npix");
    AddOption('o',"saturated_objects",2,"filename");
    AddOption('r',"growbox",2,"factor");
    AddOption('s',"bgiters",2,"n");
    AddOption('t',"trail_length",2,"npixels");
    AddOption('v',"verbose",2,"level");
    AddOption('w',"growrad",2,"value");
    AddOption('x',"trailboxes",2,"filename");
    AddOption('y',"holelevel",2,"value");
    AddOption('z',"zerostarweights");
    AddOption('D',"Debug");
    AddArgument("infile",1);
    AddArgument("outfile",1);
    AddHelp("interpolate-stars","Do radial interpolation over saturated stars for which bleedtrails were detected.");
    AddHelp("global_only","Do not use local statistics.");
    AddHelp("help","Prints this long version of help.");
    AddHelp("interpolate","Interpolate over bleedtrails.");
    AddHelp("trailreject","Reject bleedtrails of size <npix> or smaller. (1)");
    AddHelp("holelevel","Number of sigma below sky for hole detection. (3.0)");
    AddHelp("edgesize","Size of edges used to detect edgebleed. (15)");
    AddHelp("starmask","Create a mask for the detected bright objects. (No)");
    AddHelp("bgreject","Use specified <factor> as scalefactor for background rejection. (5.0)");
    AddHelp("scalefactor","Use specified <value> as scalefactor on background to detect bleedtrails. (1.0)");
    AddHelp("starlevel","Use specified <value> as scalefactor on background to detect stars. (5.0)");
    AddHelp("numtrail","Number of contiguous pixels required for trail detection. (trail_length/2)");
    AddHelp("saturated_objects","Filename for saturated star table.(None)");
    AddHelp("growbox","Factor by which to grow blob boxes for determining background level. (2)");
    AddHelp("growrad","Factor by which to grow detected star radii. (1.0)");
    AddHelp("bgiters","Number of iterations for background estimation. (1)");
    AddHelp("trail_length","Use specified <npixels> for trail structuring element in Y. (20)");
    AddHelp("verbose","Verbosity level [0-3], default = 1.");
    AddHelp("version","Print version and exit");
    AddHelp("trailboxes","Filename for trail box table.(None)");
    AddHelp("zerostarweights","Set weights to zero for masked stars.");
    AddHelp("Debug","Turn on debug/development mode.");
    //    AddHelp("starmask","Set BADPIX_STAR for bright stars from USNO-B.");
    std::ostringstream Ostr;
    Ostr << "Input FITS file.\n";
    AddArgHelp("infile",Ostr.str());
    Ostr.str("");
    Ostr << "Output FITS file.\n";
    AddArgHelp("outfile",Ostr.str());
    Ostr.str("");
    Ostr << "Detects bleedtrails in Y for incoming FITS image.\n"
	 << "Sets BADPIX_TRAIL for detected bleedtrails,\n"
	 << "and BADPIX_INTERP if interpolation is enabled,\n"
	 << "and BADPIX_STAR if star masking is enabled.\n";
    _description.assign(Ostr.str());
  };
};

///
/// @brief Attempts to determine actual star radius given a guess. 
///
/// @param mask                The bitmask plane
/// @param image               The image data
/// @param cx                  The center in image coordinates X
/// @param cy                  The center in image coordinates Y
/// @param star_r              Initial guess for stellar radius (result returned here, too)
/// @param Nx                  Number of pixels in X for image
/// @param Ny                  Number of pixels in Y for image
/// @param rejection_mask      Pixels with these bits set are not reliable
/// @param star_scalefactor    Number of sigma above which an object pixel is assumed
/// @param stats               Reference to a Morph::StatType holding image statistics
///
/// This routine makes a box that is 10 times larger than the initial guess and then
/// calculates the radial histogram.  Stepping out from the center, the first time the
/// bin median drops below that which is considered object levels, the radius is recorded
/// and returned in star_r.
///
void ModifyStarR(Morph::MaskDataType *mask,Morph::ImageDataType *image,
		 double cx,double cy,double &star_r,Morph::IndexType Nx,
		 Morph::IndexType Ny,Morph::MaskDataType rejection_mask,
		 double star_scalefactor,Morph::StatType &stats);


///
/// @brief Leave mkbleedmask processing history in header
///
/// @param Inimage Reference to a FitsTools::FitsImage object with the main (processed) image.
/// @param argv Arguments from the command line.
/// @param bleedpix Integer number of bleedtrail pixels detected.
/// @param do_interp Indicates whether bleedtrail interpolation was enabled.
/// @param do_star_interp Indicates whether star/object interpolation was enabled.
/// @param do_starmask Indicates whether star masking was enabled.
/// @param bleedmask_status Integer indicating whether there was some error in processing. 
///
int UpdateImageHeader(FitsTools::FitsImage &Inimage,const char *argv[],int nbleedpix,
		      bool do_interp,bool do_star_interp,bool do_starmask,int bleedmask_status);

///
/// @brief Writes the specified FITS image into the specified filename.
///
/// @param Inimage Reference to a FitsTools::FitsImage object with the image to write.
/// @param ofilename Reference to a string containing the filename into which to write the image.
/// @param flag_verbose Integer flag indicating verbosity level.
///
int WriteOutputImage(FitsTools::FitsImage &Inimage,const std::string &ofilename,int flag_verbose);

///
/// @brief Checks to see whether mkbleedmask has had critical errors
///
/// @param status The integer status of the processing
///
bool CheckStatus(int &status)
{
  return(status > 0);
}


///
/// @brief This is the main driver for mkbleedmask.
///
/// @param argv arguments from the command line.
///
int MakeBleedMask(const char *argv[])
{
  //  Profiler::ProfilerObj profiler;
  // ------- Setup and Command Line Parsing --------
  
  // Initialization of command line object
  MakeBleedMaskComLine comline;
  int comline_status = comline.ProcessCommandLine((const char **)argv);

  std::string program_name(comline.ProgramName());
  //  profiler.Init(program_name,0);
  //  profiler.FunctionEntry("Init");
  // Just print help on stdout and quit if user
  // did -h or --help. 
  if(comline.Narg()==1){
    std::cout << comline.ShortUsage() << std::endl;
    return(1);
  }
  if(!comline.GetOption("help").empty()){
    std::cout << comline.LongUsage() << std::endl;
    return(0);
  }
  if(!comline.GetOption("version").empty()){
    std::cout << svn_id << std::endl;
    return(0);
  }

  if(comline_status){
    std::cerr << comline.ErrorReport() << std::endl
	      << std::endl
	      << comline.ShortUsage()  << std::endl 
	      << std::endl;
    return(1);
  }
  
  // Find out what command line options were set
  bool do_interp         = !comline.GetOption("interpolate").empty();
  bool debug             = !comline.GetOption("Debug").empty();
  bool do_star_interp    = !comline.GetOption("interpolate-stars").empty();
  std::string stw        =  comline.GetOption("trail_length");
  std::string eds        =  comline.GetOption("edgesize");
  std::string hols       =  comline.GetOption("holelevel");
  std::string trj        =  comline.GetOption("trailreject");
  std::string sfac       =  comline.GetOption("scalefactor");
  std::string slfac      =  comline.GetOption("starlevel");
  std::string sgrb       =  comline.GetOption("growbox");
  std::string sgrr       =  comline.GetOption("growrad");
  std::string snum       =  comline.GetOption("numtrail");
  std::string sbgr       =  comline.GetOption("bgreject");
  std::string sverb      =  comline.GetOption("verbose");
  std::string sbgit      =  comline.GetOption("bgiters");
  bool do_starmask       = !comline.GetOption("starmask").empty();
  bool do_zero           = !comline.GetOption("zerostarweights").empty();
  bool global_stats_only = !comline.GetOption("global_only").empty();

  std::string SaturatedObjectFileName =  comline.GetOption("saturated_objects");
  std::string TrailBoxesFileName      =  comline.GetOption("trailboxes");

  int bleedmask_status = 0;

  // Parse verbosity level
  int flag_verbose = 1;
  if(!sverb.empty()){
    std::istringstream Istr(sverb);
    Istr >> flag_verbose;
    if(flag_verbose < 0 || flag_verbose > 3)
      flag_verbose = 1;
  }
  LX::ReportMessage(flag_verbose,STATUS,1,svn_id);
  
  // - Set up IO objects -
  //
  // Default to reading in from stdin
  std::ostringstream Out;
  // If user specified an input file, then switch to using
  // the specified input file.
  bool unknown_flags = false;
  std::vector<std::string> infile_names(comline.GetArgs());
  if(infile_names.size() != 2){
    Out << comline.ShortUsage() << std::endl
	<< program_name << ": Requires input file name and output file name." << std::endl;
    LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
    return(1);
  }

  // Parse normalization factor
  long trail_width = 0;
  if(!stw.empty()){
    std::istringstream Istr(stw);
    Istr >> trail_width;
    if(trail_width <= 0){
      Out << program_name << ": Invalid size for trail structuring element, " << trail_width 
	  << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }

  double scalefactor = 1.0;
  if(!sfac.empty()){
    std::istringstream Istr(sfac);
    Istr >> scalefactor;
    if(scalefactor < 0){
      Out << program_name << ": Invalid scalefactor for bleedtrail detection, " 
	  << scalefactor << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }
  
  double star_scalefactor = 5.0;
  if(!slfac.empty()){
    std::istringstream Istr(slfac);
    Istr >> star_scalefactor;
    if(star_scalefactor < 0){
      Out << program_name << ": Invalid scalefactor for star detection, " 
	  << star_scalefactor << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }
  
  Morph::ImageDataType ground_rejection_factor = 5.0;
  if(!sbgr.empty()){
    std::istringstream Istr(sbgr);
    Istr >> ground_rejection_factor;
  }
  
  Morph::ImageDataType box_growth_factor = 2.0;
  if(!sgrb.empty()){
    std::istringstream Istr(sgrb);
    Istr >> box_growth_factor;
    if(box_growth_factor < 0){
      Out << program_name << ": Invalid box growth factor, " << box_growth_factor << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }
  
  Morph::ImageDataType radius_growth_factor = 1.0;
  if(!sgrr.empty()){
    std::istringstream Istr(sgrr);
    Istr >> radius_growth_factor;
    if(radius_growth_factor < 0){
      Out << program_name << ": Invalid star radius growth factor, " << radius_growth_factor << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }
  Morph::ImageDataType rgf2 = radius_growth_factor*radius_growth_factor;

  int nbgit = 1;
  if(!sbgit.empty()){
    std::istringstream Istr(sbgit);
    Istr >> nbgit;
    if(nbgit < 0){
      Out << program_name << ": Invalid number of background iterations, " << nbgit << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }

  // This option is needed to ensure that the bleed trail size rejection is functioning properly
  int ntrail_reject = 1;
  if(!trj.empty()){
    std::istringstream Istr(trj);
    Istr >> ntrail_reject;
    if(ntrail_reject < 0){
      Out << program_name << ": Invalid trail rejection threshold, " 
	  << ntrail_reject << ", resetting to 1." << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
      ntrail_reject = 1;
    }
  }

  // This option is needed to ensure that the edgebleed detection is functioning properly
  int npix_edge = 15;
  if(!eds.empty()){
    std::istringstream Istr(eds);
    Istr >> npix_edge;
    if(npix_edge < 0){
      Out << program_name << ": Invalid edge size for edgebleed detection, " 
	  << npix_edge << ", resetting to 15." << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
      npix_edge = 15;
    }
  }

  // This option is needed to ensure that the edgebleed detection is functioning properly
  double hole_detection_level  = 3.0;
  if(!hols.empty()){
    std::istringstream Istr(hols);
    Istr >> hole_detection_level;
  }


  // Set up a list of keywords we want
  // to exclude from our headers
  char **exclusion_list;
  exclusion_list = new char * [43];
  exclusion_list[0] = (char *)"SIMPLE";
  exclusion_list[1] = (char *)"BITPIX";
  exclusion_list[2] = (char *)"NAXIS";
  exclusion_list[3] = (char *)"NAXIS1";
  exclusion_list[4] = (char *)"NAXIS2";
  exclusion_list[5] = (char *)"EXTEND";
  exclusion_list[6] = (char *)"NEXTEND";
  exclusion_list[7] = (char *)"CHECKSUM";
  exclusion_list[8] = (char *)"DATASUM";
  exclusion_list[9] = (char *)"CHECKVER";
  exclusion_list[10] =(char *)"GCOUNT";
  exclusion_list[11] =(char *)"PCOUNT";
  exclusion_list[12] =(char *)"BZERO";
  exclusion_list[13] =(char *)"BSCALE";
  exclusion_list[14] =(char *)"INHERIT";
  exclusion_list[15] =(char *)"XTENSION";
  exclusion_list[16] =(char *)"TFIELDS";
  exclusion_list[17] =(char *)"TFORM1";
  exclusion_list[18] =(char *)"ZIMAGE";
  exclusion_list[19] =(char *)"ZTILE1";
  exclusion_list[20] =(char *)"ZTILE2";
  exclusion_list[21] =(char *)"ZCMPTYPE";
  exclusion_list[22] =(char *)"ZNAME1";
  exclusion_list[23] =(char *)"ZVAL1";
  exclusion_list[24] =(char *)"ZNAME2";
  exclusion_list[25] =(char *)"ZVAL2";
  exclusion_list[26] =(char *)"EXTNAME";
  exclusion_list[27] =(char *)"ZTENSION";
  exclusion_list[28] =(char *)"ZBITPIX";
  exclusion_list[29] =(char *)"ZNAXIS";
  exclusion_list[30] =(char *)"ZNAXIS1";
  exclusion_list[31] =(char *)"ZNAXIS2";
  exclusion_list[32] =(char *)"ZPCOUNT";
  exclusion_list[33] =(char *)"ZGCOUNT";
  exclusion_list[34] =(char *)"DCREATED";
  exclusion_list[35] =(char *)"TTYPE1";
  exclusion_list[36] =(char *)"ZHECKSUM";
  exclusion_list[37] =(char *)"TTYPE2";
  exclusion_list[38] =(char *)"TTYPE3";
  exclusion_list[39] =(char *)"ZSIMPLE";
  exclusion_list[40] =(char *)"ZEXTEND";
  exclusion_list[41] =(char *)"TFORM2";
  exclusion_list[42] =(char *)"TFORM3";
  
  std::vector<std::string>::iterator inamei = infile_names.begin();

  //  profiler.FunctionExit("Init");
  //  profiler.FunctionEntry("FitsRead");
  // Step through the file lists and read in the FITS images

  // Declare and initialize FitsImage object which
  // will handle most of the FITS-specific interactions
  // under-the-covers.
  FitsTools::FitsImage Inimage;
  //  Inimage.SetProfiler(profiler);
  Inimage.SetOutStream(Out);
  Inimage.SetExclusions(exclusion_list,43);

  // Set input and output FITS filenames
  std::string filename(*inamei++);
  std::string ofilename(*inamei);
  if(flag_verbose){
    Out.str("");
    Out << "Searching for " << filename;
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  
  int status = 0;
  // If found, this will read the FITS file including
  // all data units, headers and images.
  if(Inimage.Read(filename,flag_verbose)){
    LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
    Out.str("");
    Out << "Failed to read " << filename << ".";
    LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
    return(1);
  }
  //  Inimage.Copy(ofilename,flag_verbose);
  // Close the input FITS file, we're done with it.
  Inimage.Close();
  // Let stdout know what file was actually read.
  // (it can be different due to the various valid extensions)
  if(flag_verbose){
    Out.str("");
    Out << "Found " << Inimage.DES()->name << ".";
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  // How many HDU's were read?  Report it to stdout if verbose
  int number_of_hdus = Inimage.Headers().size();
  std::vector<std::string>::iterator hi = Inimage.Headers().begin();
  if(flag_verbose){
    Out.str("");
    Out << "Read " << number_of_hdus << " data units from " 
	<< filename;
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  // If *really* verbose, dump the headers to stdout
  if(flag_verbose>4){
    while(hi != Inimage.Headers().end())
      {
	Out.str("");
	Out << "Size for image: " << Inimage.DES()->axes[0] << "X" 
	    << Inimage.DES()->axes[1] << ", number of pixels = " 
	    << Inimage.DES()->npixels << std::endl
	    << "Header " << hi-Inimage.Headers().begin()+1 << std::endl
	    << *hi << std::endl;
	LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
	hi++;
      }
  }

  //  profiler.FunctionExit("FitsRead");
  
  //  profiler.FunctionEntry("Setup");
  
  Morph::IndexType Nx = Inimage.DES()->axes[0];
  Morph::IndexType Ny = Inimage.DES()->axes[1];
  Morph::IndexType npix = Nx*Ny;
  Morph::IndexType Npix = npix;

  std::vector<Morph::MaskDataType> original_mask(npix);
  for(int omi = 0;omi < npix;omi++)
    original_mask[omi] = Inimage.DES()->mask[omi];
  
  // Flag bad values in the image (e.g. INF, NAN) - shouldn't have to do this,
  // but bad values in input images really mess up the statistics and cause
  // bleed detection to perform poorly.
  Morph::IndexType nbadpix = Morph::MaskBadPixelData(Inimage.DES()->image,
						     Inimage.DES()->mask,Nx,Ny,BADPIX_BPM);
  if(nbadpix){
    // Set weights to zero for all of BPM 
    if(Inimage.DES()->varim){
      for(Morph::IndexType nni = 0;nni < npix;nni++){
	if(Inimage.DES()->mask[nni]&BADPIX_BPM)
	  Inimage.DES()->varim[nni] = 0;
      }
    }
    Out.str("");
    Out << "Warning: Flagged " << nbadpix << " bad-valued pixels and set weights to zero."; 
    LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
    Out.str("");
  }

  // Set up trail_width - which is the length of the structuring element
  // in the scan direction. 
  // Check for FWHM keyword if width unset
  float fwhm = 0;
  if(trail_width == 0){
    fwhm = FitsTools::GetHeaderValue<float>(Inimage.ImageHeader(),"FWHM");
    if(fwhm > 0){
      fwhm *= 4;
      if(trail_width == 0)
	trail_width = static_cast<long>(fwhm);
    }
  }
  // Default to 20
  if(trail_width == 0)
    trail_width = 20;
  // Force symmetry
  if(trail_width%2)
    trail_width++;
  
  // The trail_arm is the number of pixels in each of +/- scan direction
  Morph::IndexType trail_arm = trail_width/2;
  if(!snum.empty()){
    std::istringstream Istr(snum);
    Istr >> trail_arm;
    if(trail_arm < 2){
      Out << program_name << ": Invalid number of pixels for trail detection, must be > 2. (" 
	  << trail_arm << ")"  << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }    
  }

  // This actually forms the structuring element for the trail
  // search, trail_structure.
  Morph::StructuringElementType trail_structuring_element;
  unsigned int trail_pix = trail_width/2;
  std::vector<long> trail_structure(trail_width,0);
  std::vector<long>::iterator trailit = trail_structure.begin();
  for(unsigned int i = trail_pix;i > 0;i--){
    *trailit++ = -static_cast<long>(Nx*i);
    Morph::ImageOffsetType filter_offset;
    filter_offset.first = 0;
    filter_offset.second = -i;
    trail_structuring_element.push_back(filter_offset);
  }
  //  Morph::ImageOffsetType fo;
  //  fo.first = fo.second = 0;
  //  trail_structuring_element.push_back(fo);
  for(unsigned int i = 1;i <= trail_pix;i++){
    *trailit++ = i*Nx;
    Morph::ImageOffsetType filter_offset;
    filter_offset.first = 0;
    filter_offset.second = i;
    trail_structuring_element.push_back(filter_offset);
  }

  
  int ccdnum = 0;
  ccdnum = FitsTools::GetHeaderValue<int>(Inimage.ImageHeader(),"CCDNUM");
  if(ccdnum < 0)
    ccdnum = 0;
  
  // Pixels with these bits set are ignored in bleedtrail detection
  short trail_rejection_mask  = BADPIX_CRAY | BADPIX_BPM;
  // Accept pixels with this bit set no matter what
  short trail_exception_mask  = 0;
  // Set this bit for detected trail pixels
  short trail_mask = BADPIX_TRAIL;
  // Reject pixels with these bits set when determining background levels
  short ground_rejection_mask = trail_rejection_mask | trail_mask | BADPIX_SATURATE | BADPIX_LOW;

 
  int bleed_status = 0; 
  int total_trail_pixels = 0;
  int method = 1;
  

  // Copy the incoming mask (or create one if it doesn't already exist)
  std::vector<Morph::MaskDataType> temp_mask(npix,0);
  if(!Inimage.DES()->mask){
    method = 0;
    Inimage.DES()->mask = &temp_mask[0];
    if(flag_verbose){
      Out.str("");
      Out << program_name << "::Warning: Incoming image had no mask, using temporary.";
      LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
    }
  } else {
    std::vector<Morph::MaskDataType>::iterator  tmi = temp_mask.begin();
    while(tmi != temp_mask.end()){
      Morph::IndexType index = tmi - temp_mask.begin();
      *tmi++ = Inimage.DES()->mask[index];
    }
  }
  
  // Dilate the BADPIX_SATURATE mask so that pixels that neighbor saturated
  // pixels are also marked as "saturated". At this point, the temp mask is
  // just a marker that tells the rest of the code where to look for bleed
  // trails, and this is a way to specify that pixels neighboring saturated
  // pixels are untrusted.
  //
  //  profiler.FunctionEntry("Dilation");
  std::vector<Morph::IndexType> blob_image(npix,0);
  std::vector<std::vector<Morph::IndexType> > saturated_blobs;
  if(debug){
    Morph::GetBlobsWithRejection(&temp_mask[0],Nx,Ny,BADPIX_SATURATE,BADPIX_BPM,
				 0,blob_image,saturated_blobs);
    std::vector<Morph::BlobType>::iterator bbbi   = saturated_blobs.begin();
    std::cout << "SATURATED BLOBS BEFORE DILATION:" << std::endl;
    while(bbbi != saturated_blobs.end()){
      int blobno = bbbi - saturated_blobs.begin() + 1;
      Morph::BlobType &blob = *bbbi++;
      int nblobpix = blob.size();
      std::cout << "Blob(" << blobno << ") size = " << nblobpix << std::endl;
    }
  } 
  // The structuring element used for dilation of the saturated pixel mask
  std::vector<long> structuring_element(8,0);
  structuring_element[0] = -(Nx+1);
  structuring_element[1] = -Nx;
  structuring_element[2] = -Nx + 1;
  structuring_element[3] = -1;
  structuring_element[4] = 1;
  structuring_element[5] = Nx-1;
  structuring_element[6] = Nx;
  structuring_element[7] = Nx+1;
  // Dilate twice.  This *should* speed up the bleedtrail detection because
  // now instead of using the bounding box, just the "saturated" mask can be 
  // used.
  Morph::DilateMask(&temp_mask[0],Nx,Ny,structuring_element,BADPIX_SATURATE);
  Morph::DilateMask(&temp_mask[0],Nx,Ny,structuring_element,BADPIX_SATURATE);

  // profiler.FunctionExit("Dilation");
  // Get the "blobs" of BADPIX_SATURATE pixels. The number of blobs is the 
  // initial number of saturated objects before bleedtrail rejection.
  //
  // - blob_image is an integer image wherein the value indicates which
  //   blob to which the pixel belongs.[1:nblobs] 0 = not part of a blob
  // - blobs contains the pixel indices for each blob [0:npix-1]
  //
  //  std::vector<Morph::IndexType> blob_image(npix,0);
  //  std::vector<std::vector<Morph::IndexType> > blobs;
  blob_image.resize(npix,0);
  saturated_blobs.resize(0);
  // profiler.FunctionEntry("GetBlobs");
  // This call populates the above two data structures
  Morph::GetBlobsWithRejection(&temp_mask[0],Nx,Ny,BADPIX_SATURATE,BADPIX_BPM,
			       0,blob_image,saturated_blobs);

  // profiler.FunctionExit("GetBlobs");

  if(debug){
    std::vector<Morph::BlobType>::iterator bbbi   = saturated_blobs.begin();
    std::cout << "AFTER DILATION:" << std::endl;
    while(bbbi != saturated_blobs.end()){
      int blobno = bbbi - saturated_blobs.begin() + 1;
      Morph::BlobType &blob = *bbbi++;
      int nblobpix = blob.size();
      std::cout << "Blob(" << blobno << ") size = " << nblobpix << std::endl;
    } 
  }
  if(flag_verbose){
    Out.str("");
    Out << "Found " << saturated_blobs.size() << " saturated blobs initially."; 
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
    Out.str("");  
  }
  
    
  // profiler.FunctionEntry("Statistics");
  // These data structures hold the bounding boxes for
  // each blob, and the image statistics within an extended
  // version of the bounding box.
  std::vector<Morph::StatType> box_stats;
  std::vector<Morph::BoxType>  blob_boxes;
  std::vector<Morph::BoxType>  bg_boxes;
  std::vector<Morph::BoxType>  candidate_trail_boxes;
  
  // Data structures for global image statistics
  Morph::StatType image_stats;
  Morph::BoxType image_box(4,0);
  image_box[0] = 0;
  image_box[1] = Nx-1;
  image_box[2] = 0;
  image_box[3] = Ny-1;

  Morph::IndexType image_npix_box;
  Morph::IndexType minpix = static_cast<Morph::IndexType>(.3*static_cast<double>(npix));
  Morph::IndexType niter = 0;
  
  // profiler.FunctionEntry("GetSkyB");
  // This determines the sky by rejecting bright pixels until the mean stops changing.
  if(Morph::GetSky(Inimage.DES()->image,&temp_mask[0],Nx,Ny,
  		   minpix,nbgit,ground_rejection_factor,1e-3,
  		   BADPIX_SATURATE | BADPIX_CRAY | BADPIX_BPM | BADPIX_STAR | 
  		   BADPIX_TRAIL,BADPIX_INTERP,image_stats,image_npix_box,niter,&Out)){
    std::ostringstream GSO;
    GSO << "Global image statistics failed:" << std::endl
	<< Out.str() << std::endl;
    LX::ReportMessage(flag_verbose,STATUS,4,GSO.str());
    bleedmask_status = 1;
    UpdateImageHeader(Inimage,argv,0,do_interp,do_star_interp,
		      do_starmask,bleedmask_status);
    return(WriteOutputImage(Inimage,ofilename,flag_verbose));
  } else if(flag_verbose) {
    Out << "Global image statistics converged after NITER=" << niter 
	<< (niter==1 ? " iteration." : " iterations.") << std::endl; 
    LX::ReportMessage(flag_verbose,QA,1,Out.str());
    Out.str("");
    Out << Util::stripdirs(filename) << "   BKGD=" << image_stats[Image::IMMEAN] << ", "
	<< "BKGD_SIGMA=" << image_stats[Image::IMSIGMA] << std::endl;
    LX::ReportMessage(flag_verbose,QA,1,Out.str());
  }
  // profiler.FunctionExit("GetSky");


  // ******** Edge Bleed Detection ************
  //
  // This vector will store the bounding boxes
  // for the edgebleeds (if any).
  std::vector<Morph::BoxType> edgebleed_boxes;
  
  // These y values will indicate the extent of
  // of the edgebleed damage (if any).
  std::vector<Morph::BoxType> edgebleed_detection_boxes; 
  Morph::IndexType y00 = 0;
  Morph::IndexType y01 = Ny;
  Morph::IndexType y10 = 0;
  Morph::IndexType y11 = Ny;
  Morph::IndexType ampx = Nx/2;
  Morph::IndexType nrows_detection = 20;

  // Using arbitrary number (200) here to be the max size
  // of an edgebleed bounding box's extent in Y.  
  Morph::IndexType ytop = Ny - (npix_edge+nrows_detection);
  Morph::IndexType ybottom = (npix_edge+nrows_detection);
  
  // For indicating whether edgebleed was 
  // actually detected on a given edge.
  bool bl_bleed = false;
  bool tl_bleed = false;
  bool br_bleed = false;
  bool tr_bleed = false;

  // -- These configuration values are used for  --
  // -- edgebleed detection and processing
  //  Morph::ImageDataType hole_detection_level = 3.0;
  short BADPIX_HOLE = BADPIX_TRAIL;
  int NPIX_EDGEBLEED = 500; // how large should the "low" section be to get caught
  int NPIXMAX_EDGEBLEED = 200000; // larger than this? probably not an edgebleed
  int XSIZE_EDGEBLEED = 100; // must be extended at least this far in X to be caught
  // --------------------------

  // define the levels at which the threshold the image
  // for detecting low pixels
  double hole_level = image_stats[Image::IMMEAN] - 5.0*image_stats[Image::IMSIGMA];
  double high_level = std::numeric_limits<Morph::ImageDataType>::max();

  Morph::BlobsType SatBlobs;
  std::vector<Morph::IndexType> SatBlobImage(npix,0);
  Morph::GetBlobs(&temp_mask[0],Nx,Ny,BADPIX_SATURATE,SatBlobImage,SatBlobs);
  SatBlobImage.resize(0);
  Morph::BlobsType::iterator satblobit = SatBlobs.begin();
  bool suspect_edge_bleed = false;
  bool suspect_edge_bleed_bottom = false;
  bool suspect_edge_bleed_top = false;
  bool read_register_on_bottom = ((ccdnum > 31) || (ccdnum == 0));
  bool read_register_on_top    = ((ccdnum < 32) || (ccdnum == 0));
  Out << "Found read registers on " 
      << (read_register_on_top ? (read_register_on_bottom ? 
				  "top *and* bottom (i.e. ccdnum not resolved)." : "top.") : "bottom.");
  LX::ReportMessage(flag_verbose,STATUS,1,Out.str());  

  while((satblobit != SatBlobs.end()) && !suspect_edge_bleed){
    Morph::BlobType &satblob(*satblobit++);
    Morph::BoxType box;
    std::sort(satblob.begin(),satblob.end());
    Morph::GetBlobBoundingBox(satblob,Nx,Ny,box);
    if((box[2] <= 20) && read_register_on_bottom)
      suspect_edge_bleed_bottom = true;
    if((box[3] >= (Ny - 20)) && read_register_on_top)
      suspect_edge_bleed_top = true;
  }
  suspect_edge_bleed = (suspect_edge_bleed_bottom || suspect_edge_bleed_top); 
  if(suspect_edge_bleed){
    // Give a quick status message to indicate potential edgebleed.
    if(flag_verbose){
      Out.str("");
      Out << "Checking for possible edgebleed on "
	  << ((suspect_edge_bleed_bottom && suspect_edge_bleed_top) ? "both edges " : 
	      (suspect_edge_bleed_bottom ? "bottom edge " : "top edge "))
	  << "because saturated blob(s) collide(s) with read registers.";
      LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
    }
    // Look through the rows by the read registers and find out
    // if there are some rows with large numbers of negative pixels.
    std::vector<Morph::IndexType> ndip(2,0);
    Morph::IndexType bottom_y = npix_edge;
    Morph::IndexType top_y = Ny - npix_edge - 1;
    Morph::IndexType bottom_row = (read_register_on_bottom ? 
				   (read_register_on_top ? top_y : bottom_y) : (top_y-14));
    Morph::IndexType top_row    = (read_register_on_top ? 
				   (read_register_on_bottom ? bottom_y : top_y) : (bottom_y+14));
    for(Morph::IndexType rowin = bottom_row;rowin <= top_row;rowin++){
      Morph::IndexType index_base = rowin*Nx;
      for(Morph::IndexType pixin = 0;pixin < ampx;pixin++){
	if(Inimage.DES()->image[index_base+pixin] < hole_level) ndip[0]++;
	if(Inimage.DES()->image[index_base+ampx+pixin] < hole_level) ndip[1]++;
	index_base++;
      }
    }
    bool detected_dip_left  = (ndip[0] > 500);
    bool detected_dip_right = (ndip[1] > 500);
    bool strong_dip_left    = (ndip[0] > 1000);
    bool strong_dip_right   = (ndip[1] > 1000);
    if(detected_dip_left || detected_dip_right){
      Out.str("");
      Out << "Number of Low Pixels near read registers: (" << ndip[0] << "," 
	  << ndip[1] << ")";
      LX::ReportMessage(flag_verbose,STATUS,1,Out.str());  
    } 
    if(!(strong_dip_left || strong_dip_right)){
      hole_detection_level *= 5.0;
      Out.str("");
      Out << "No edgebleed dips detected near read registers, resetting hole detection level to "
	  << hole_detection_level << ".";
      LX::ReportMessage(flag_verbose,STATUS,1,Out.str());  
    }
    hole_level = image_stats[Image::IMMEAN] - hole_detection_level*image_stats[Image::IMSIGMA];
    // This sets BADPIX_LOW (in a temporary mask) for the low pixels
    std::vector<Morph::MaskDataType> hole_mask(npix,0);
    std::vector<Morph::IndexType> hole_image(npix,0);
    std::vector<std::vector<Morph::IndexType> > hole_blobs;
    Morph::ThresholdImage(Inimage.DES()->image,&hole_mask[0],Nx,Ny,
			  hole_level,high_level,BADPIX_LOW,BADPIX_SATURATE);  

    // Erosion followed by Dilation will removed straggly connections 
    // but preserve the effective extent of the bounding box for the large
    // part of the "blob" - hopefully help to clean up these vastly overflagged
    // edgebleeds.
    Morph::ErodeMask(&hole_mask[0],Nx,Ny,structuring_element,BADPIX_LOW);
    Morph::DilateMask(&hole_mask[0],Nx,Ny,structuring_element,BADPIX_LOW);

    // Update the temporary mask so that it knows about the image's BPM
    for(int jj = 0;jj < npix;jj++)
      if(Inimage.DES()->mask[jj] & BADPIX_BPM)
	hole_mask[jj] |= BADPIX_BPM;
    // Get the blobs of BADPIX_LOW while respecing the BPM
    Morph::GetBlobsWithRejection(&hole_mask[0],Nx,Ny,BADPIX_LOW,BADPIX_BPM|BADPIX_SATURATE,
				 0,hole_image,hole_blobs);
  
    // Loop through the blobs and if a large (i.e. larger than NPIX_EDGEBLEED)
    // blob is found, then check whether it collides with the image boundary.
    // If so, it is a positive edgebleed detection;
    // If not, then warn about strange low section in image
    std::vector<Morph::BlobType>::iterator hbi   = hole_blobs.begin();
    Morph::IndexType edgebleed_y_lower = 0;
    Morph::IndexType edgebleed_y_upper = 0;
    while(hbi != hole_blobs.end()){
      int blobno = hbi - hole_blobs.begin() + 1;
      Morph::BlobType &blob = *hbi++;
      int nblobpix = blob.size();
      std::vector<Morph::IndexType> box;
      // Sort the blob before bounding box determination
      std::sort(blob.begin(),blob.end());
      Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
      int xsize = box[1] - box[0];
      if((((nblobpix > NPIX_EDGEBLEED) && (xsize > XSIZE_EDGEBLEED)) || (detected_dip_left || detected_dip_right))
	 && (nblobpix < NPIXMAX_EDGEBLEED)){ // blob is of proper size
	if(box[0] < ampx){ // "left" side of chip
	  if((box[2] <= ybottom) && (box[3] < ytop) && read_register_on_bottom){ // box collides with bottom image boundary (and not the top)
	    bl_bleed = true;
	    if(box[3] > y00)
	      y00 = box[3];
	  } else if((box[3] >= ytop) && (box[2] > ybottom) && read_register_on_top){ // box collides with top image boundary (and not the bottom)
	    tl_bleed = true;
	    if(box[2] < y01)
	      y01 = box[2];
	  } else { // box didn't collide with the image boundary; or it collided with both (it's weird)
	    // 	    Out.str("");
	    // 	    Out << "Warning: Detected anomalous low section in image near ("
	    // 		<< (box[0]+box[1])/2.0 << "," << (box[3]+box[2])/2.0 << ")";
	    //	    LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	  }
	} 
	if(box[1] >= ampx) { // "right" side of chip
	  if((box[2]  <= ybottom) && (box[3] < ytop) && read_register_on_bottom){ // box collides with bottom image boundary (and not the top)
	    br_bleed = true;
	    if(box[3] > y10)
	      y10 = box[3];
	  } else if((box[3] >= ytop) && (box[2] > ybottom) && read_register_on_top){ // box collides with top image boundary
	    tr_bleed = true;
	    if(box[2] < y11)
	      y11 = box[2];
	  } else{ // box didn't collide with the image boundary; it's weird
	    // 	    Out.str("");
	    // 	    Out << "Warning: Detected anomalous low section in image near ("
	    // 		<< (box[0]+box[1])/2.0 << "," << (box[3]+box[2])/2.0 << ")";
	    // 	    LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	  }
	}
      }
    }
    // ******** End Edge Bleed Detection ************
  
    // If possible edge bleed is detected, then recalculate image statistics
    // with edge bleed masked out - and then make sure to use global
    // image statistics and *not* local statistics for the rest of the
    // processing
    if(bl_bleed || tl_bleed || br_bleed || tr_bleed){
      // Dilate twice again if edge bleed.
      //    Morph::DilateMask(&temp_mask[0],Nx,Ny,structuring_element,BADPIX_SATURATE);
      //    Morph::DilateMask(&temp_mask[0],Nx,Ny,structuring_element,BADPIX_SATURATE);
      global_stats_only = true;
      std::vector<Morph::MaskDataType> temp_mask2(temp_mask.begin(),temp_mask.end());
      if(bl_bleed){
	Out.str("");
	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (0:" 
	    << ampx-1 << ",0:" << y00 << ")";
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	std::vector<Morph::IndexType> box(4,0);
	box[0] = 0;
	box[1] = ampx-1;
	box[2] = 0;
	box[3] = y00;
	edgebleed_boxes.push_back(box);
	for(Morph::IndexType by = box[2];by <= box[3];by++){
	  Morph::IndexType iminy = by*Nx;
	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	    Morph::IndexType imin = iminy + bx;
	    temp_mask2[imin] |= BADPIX_LOW;
	  }
	} 
      }
      if(tl_bleed){
	Out.str("");
	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (0:" << ampx-1 
	    << "," << y01 << ":" << Ny-1 << ").";
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	std::vector<Morph::IndexType> box(4,0);
	box[0] = 0;
	box[1] = ampx-1;
	box[2] = y01;
	box[3] = Ny-1;
	edgebleed_boxes.push_back(box);      
	for(Morph::IndexType by = box[2];by <= box[3];by++){
	  Morph::IndexType iminy = by*Nx;
	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	    Morph::IndexType imin = iminy + bx;
	    temp_mask2[imin] |= BADPIX_LOW;
	  }
	} 
      }
      if(br_bleed){
	Out.str("");
	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (" << ampx << ":" 
	    << Nx-1 << ",0:" << y10 << ").";
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());    
	std::vector<Morph::IndexType> box(4,0);
	box[0] = ampx;
	box[1] = Nx-1;
	box[2] = 0;
	box[3] = y10;
	edgebleed_boxes.push_back(box);
	for(Morph::IndexType by = box[2];by <= box[3];by++){
	  Morph::IndexType iminy = by*Nx;
	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	    Morph::IndexType imin = iminy + bx;
	    temp_mask2[imin] |= BADPIX_LOW;
	  }
	} 
      }
      if(tr_bleed){
	Out.str("");
	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (" << ampx << ":" 
	    << Nx-1 << "," << y11 << ":" << Ny - 1 << ")."; 
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());    
	std::vector<Morph::IndexType> box(4,0);
	box[0] = ampx;
	box[1] = Nx-1;
	box[2] = y11;
	box[3] = Ny-1;
	edgebleed_boxes.push_back(box);      
	for(Morph::IndexType by = box[2];by <= box[3];by++){
	  Morph::IndexType iminy = by*Nx;
	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	    Morph::IndexType imin = iminy + bx;
	    temp_mask2[imin] |= BADPIX_LOW;
	  }
	} 
      }
      Out.str("");
      if(Morph::GetSky(Inimage.DES()->image,&temp_mask2[0],Nx,Ny,
		       minpix,nbgit,ground_rejection_factor,1e-3,
		       BADPIX_SATURATE | BADPIX_CRAY | BADPIX_BPM | BADPIX_STAR | 
		       BADPIX_TRAIL | BADPIX_LOW,BADPIX_INTERP,image_stats,image_npix_box,niter,&Out)){
	std::ostringstream GSO;
	GSO << "Global image statistics failed:" << std::endl
	    << Out.str() << std::endl;
	LX::ReportMessage(flag_verbose,STATUS,4,GSO.str());
	bleedmask_status = 1;
	UpdateImageHeader(Inimage,argv,0,do_interp,do_star_interp,
			  do_starmask,bleedmask_status);
	return(WriteOutputImage(Inimage,ofilename,flag_verbose));
      } else if(flag_verbose){
	Out << "Global image statistics (recalculated due to strongly suspected edgebleed) converged (again) after NITER=" << niter 
	    << (niter==1 ? " iteration." : " iterations.") << std::endl; 
	LX::ReportMessage(flag_verbose,QA,1,Out.str());
	Out.str("");
	Out << Util::stripdirs(filename) << "   BKGD=" << image_stats[Image::IMMEAN] << ", "
	    << "BKGD_SIGMA=" << image_stats[Image::IMSIGMA] << std::endl;
	LX::ReportMessage(flag_verbose,QA,1,Out.str());
      }
//       bl_bleed = false;
//       br_bleed = false;
//       tl_bleed = false;
//       tr_bleed = false;
    }
  } // suspect_edge_bleed

  // Experimental - push saturated,interpolated pixels back up to detectable levels. This undos 
  // imcorrect's interpolation over 1-pixel gaps.
  double trail_level = image_stats[Image::IMMEAN] + 2.0*scalefactor*image_stats[Image::IMSIGMA];
  //   //  double trail_level = image_stats[Image::IMMAX];
  for(int pixi = 0;pixi < Nx*Ny;pixi++){
    if((Inimage.DES()->mask[pixi]&BADPIX_SATURATE)&&(Inimage.DES()->mask[pixi]&BADPIX_INTERP)){
      Inimage.DES()->mask[pixi] ^= BADPIX_INTERP;
      temp_mask[pixi]           ^= BADPIX_INTERP;
      Inimage.DES()->image[pixi] = trail_level;
    }
  }
  
  
  
  // Copy the input image into temp_image
  std::vector<Morph::ImageDataType> temp_image(npix,0);
  std::vector<Morph::ImageDataType>::iterator tii = temp_image.begin();
  while(tii != temp_image.end()){
    Morph::IndexType index = tii - temp_image.begin();
    *tii++ = Inimage.DES()->image[index];
  }
  
  // profiler.FunctionEntry("LocalStats");
  // Loop through each blob, get image stats inside extended 
  // blob bounding box.
  //
  std::vector<Morph::BlobType>::iterator bbbi   = saturated_blobs.begin();
  while(bbbi != saturated_blobs.end()){
    int blobno = bbbi - saturated_blobs.begin() + 1;
    Morph::BlobType &blob = *bbbi++;
    int nblobpix = blob.size();
    std::vector<Morph::IndexType> box;
    std::sort(blob.begin(),blob.end());
    Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
    blob_boxes.push_back(box);
    int bleed_over = (trail_width/2) - trail_arm + 1;
    box[2] = box[2] - static_cast<Morph::IndexType>(box_growth_factor*bleed_over);
    box[3] = box[3] + static_cast<Morph::IndexType>(box_growth_factor*bleed_over);
    if(box[2] < 0)   box[2] = 0;
    if(box[3] >= Ny) box[3] = Ny - 1;
    candidate_trail_boxes.push_back(box);
    
    Morph::IndexType xsize = box[1] - box[0] + 1;
    box[0] = box[0] - static_cast<Morph::IndexType>(box_growth_factor*xsize);
    box[1] = box[1] + static_cast<Morph::IndexType>(box_growth_factor*xsize);
    if(box[0] < 0)  box[0] = 0;
    if(box[1] >= Nx) box[1] = Nx - 1;
    bg_boxes.push_back(box);
    if(debug){
      std::cout << "BackGround Box(" << blobno << "): [" << box[0] << ":" << box[1] << ","
		<< box[2] << ":" << box[3] << "]" << std::endl;
    }
    Morph::IndexType npix_box;
    Morph::StatType stats;
    bool use_image_stats = false;
    if(global_stats_only){
      use_image_stats = true;
    } else {
      Morph::BoxStats(Inimage.DES()->image,&temp_mask[0],Nx,Ny,box,
		      BADPIX_SATURATE | BADPIX_CRAY |            // reject pixels with
		      BADPIX_BPM | BADPIX_STAR | BADPIX_TRAIL,   // these bits set
		      BADPIX_INTERP,                             // accept these
		      stats,npix_box);
      if(npix_box < 300)
	use_image_stats = true;
      else{
	if(debug){
	  std::cout << "Box had " << npix_box << " valid image pixels for statistics." << std::endl;
	  std::cout << "Initial Stats: (" << stats[Image::IMMIN] << "," << stats[Image::IMMAX]
		    << "," << stats[Image::IMMEAN] << "," << stats[Image::IMSIGMA] << std::endl;
	}
	if(flag_verbose >= 3){
	  Out.str("");
	  Out << "Getting local statistics for blob " << blobno;
	  LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
	}
	bool stats_converged = false;
	Morph::IndexType bbb = 0;
	double last_mean = stats[Image::IMMEAN];
	while((bbb < nbgit) && !stats_converged && !use_image_stats){
	  bbb++;
	  Morph::ImageDataType local_ground_rejection_level = stats[Image::IMMEAN] + 
	    ground_rejection_factor*stats[Image::IMSIGMA];
	  // Get stats again, but reject high pixels
	  Morph::BoxStats(Inimage.DES()->image,&temp_mask[0],Nx,Ny,box,
			  BADPIX_SATURATE | BADPIX_CRAY |                  // reject pixels with
			  BADPIX_BPM | BADPIX_STAR | BADPIX_TRAIL,         // these bits set
			  BADPIX_INTERP,                                   // accept these
			  0,local_ground_rejection_level,                  // rejection levels
			  stats,npix_box);
	  if(npix_box < 300){
	    use_image_stats = true;
	    if(flag_verbose >= 3){
	      Out.str("");
	      Out << "Reverting to global image statistics with " << npix_box 
		  << " valid pixels in local box."; 
	      LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	    }
	  }
	  double residual = std::abs(stats[Image::IMMEAN] - last_mean);
	  if(residual < 1e-3){
	    stats_converged = true;
	    if(flag_verbose >= 3){
	      Out.str("");
	      Out << "Local image statistics converged after " << bbb << " iteration"
		  << (bbb==1 ? "." : "s."); 
	      LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
	    }
	  }
	  last_mean = stats[Image::IMMEAN];
	}
      }
    }
    if(use_image_stats){
      stats = image_stats;
      // If the user wanted local statistics, warn that they could not be used
      // for a given blob.
      if(flag_verbose && !global_stats_only){
	Out.str("");
	Out << "Global image statistics used for blob " << blobno << ".";
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
      }
    }
    // Record statistics for this box in persistent data structure
    box_stats.push_back(stats);
  }
  // profiler.FunctionExit("LocalStats");
  // profiler.FunctionExit("Statistics");

  if(flag_verbose){
    Out.str("");
    Out << program_name  
	<< " using " << trail_width 
	<< " in the trail direction, and " << std::endl
	<< " searching for trails of " << trail_arm 
	<< " contiguous bright ( >" << scalefactor << " sigma) pixels, " 
	<< std::endl
	<< " and detecting stars at " << star_scalefactor << " sigma."; 
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  


  // profiler.FunctionExit("Setup");
    


  // profiler.FunctionEntry("BleedTrails");
  BleedTrailDetectionFilter bleed_trail_detection_filter;
  bleed_trail_detection_filter.SetStructuringElement(trail_structuring_element);
  bleed_trail_detection_filter.SetInputImageData(Inimage.DES()->image);
  bleed_trail_detection_filter.SetInputMaskData(&temp_mask[0]);
  bleed_trail_detection_filter.SetOutputImageData(&temp_image[0]);
  bleed_trail_detection_filter.SetOutputMaskData(Inimage.DES()->mask);
  bleed_trail_detection_filter.SetNumberOfPixels(npix);
  bleed_trail_detection_filter.SetNumberOfPixelsInX(Nx);
  bleed_trail_detection_filter.SetTrailDetectionThreshold(trail_arm);
  bleed_trail_detection_filter.SetTrailMask(trail_mask);
  bleed_trail_detection_filter.SetStarMask(BADPIX_STAR);
  bleed_trail_detection_filter.SetInterpolatedMask(BADPIX_INTERP);
  bleed_trail_detection_filter.SetProcessingRejectionMask(trail_rejection_mask);
  bleed_trail_detection_filter.SetDataRejectionMask(trail_rejection_mask);
  std::vector<Morph::ImageDataType> trail_levels;
  std::vector<Morph::ImageDataType> star_levels;
  std::vector<Morph::StatType>::iterator statit = box_stats.begin();
  while(statit != box_stats.end()){
    Morph::StatType &bstat(*statit++);
    trail_levels.push_back(bstat[Image::IMMEAN]+scalefactor*bstat[Image::IMSIGMA]);
    star_levels.push_back(bstat[Image::IMMEAN]+star_scalefactor*bstat[Image::IMSIGMA]);
  }
  bleed_trail_detection_filter.SetStarNoise(image_stats[Image::IMSIGMA]/5.0);
  bleed_trail_detection_filter.SetBackgroundNoise(image_stats[Image::IMSIGMA]/2.0);
  bleed_status = bleed_trail_detection_filter.DetectBleedsOnBlobs(saturated_blobs,trail_levels);
  bleed_trail_detection_filter.SetDataRejectionMask(trail_rejection_mask|trail_mask);
  std::vector<Morph::MaskDataType>::iterator  tmi = temp_mask.begin();
  while(tmi != temp_mask.end()){
    Morph::IndexType index = tmi - temp_mask.begin();
    *tmi++ = Inimage.DES()->mask[index];
  }
  //  std::vector<Morph::MaskDataType> tempmask2(temp_mask.begin(),temp_mask.end());
  //  bleed_trail_detection_filter.SetOutputMaskData(&tempmask2[0]);
  int star_status = bleed_trail_detection_filter.DetectStarsOnBlobs(saturated_blobs,star_levels);
  tmi = temp_mask.begin();
  while(tmi != temp_mask.end()){
    Morph::IndexType index = tmi - temp_mask.begin();
    *tmi++ = Inimage.DES()->mask[index];
  }
  // This will close up small 1-pixel gaps in the trail.
  Morph::BlobsType trail_blobs;
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,trail_mask,blob_image,trail_blobs);
  Morph::BlobsType::iterator tbi = trail_blobs.begin();
  while(tbi != trail_blobs.end()){
    Morph::BoxType box;
    Morph::BlobType &trailblob(*tbi++);
    std::sort(trailblob.begin(),trailblob.end());
    Morph::GetBlobBoundingBox(trailblob,Nx,Ny,box);
    for(int tj = box[2];tj <= box[3];tj++){
      for(int ti = box[0];ti <= box[1];ti++){
	if(((ti < (Nx-1)) && (ti > 0))){
	  int pixel_index    = tj*Nx+ti;
	  if(!(Inimage.DES()->mask[pixel_index]&BADPIX_TRAIL)){
	    int left_neighbor  = pixel_index - 1;
	    int right_neighbor = pixel_index + 1;
	    if(Inimage.DES()->mask[left_neighbor]&BADPIX_TRAIL &&
	       Inimage.DES()->mask[right_neighbor]&BADPIX_TRAIL){
	      Inimage.DES()->mask[pixel_index] |= BADPIX_TRAIL;
	      bleed_status++;
	    }
	  }
	}
      }
    }
  }
  trail_blobs.resize(0);
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,trail_mask,blob_image,trail_blobs);
  int interp_status = bleed_trail_detection_filter.LinearInterpolateOverTrailsInBlobs(trail_blobs);
  


  //   bleed_status = DetectBleedTrailsInBlobs(Inimage.DES()->image,
  //   					  &temp_image[0],
  //   					  &temp_mask[0],
  //   					  Inimage.DES()->mask,
  //   					  Nx,Ny,blobs,box_stats,
  //   					  trail_structure,
  //   					  trail_rejection_mask,trail_exception_mask,
  //   					  ground_rejection_mask,trail_exception_mask,
  //   					  trail_mask,
  //   					  BADPIX_INTERP,
  //   					  BADPIX_STAR,
  //   					  trail_arm,scalefactor,
  //   					  star_scalefactor,
  //   					  do_interp,std::cout);
  
  // profiler.FunctionExit("BleedTrails");
  total_trail_pixels = bleed_status;

  if(bleed_status < 0)
    { 
      LX::ReportMessage(flag_verbose,STATUS,4,Out.str());
      Out.str("");
      Out << "DetectBleedTrails failed on pixel " << bleed_status << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,4,Out.str());
      // Return non-error exit code to avoid disturbing the pipeline
      return(0);
    }

  if(flag_verbose){
    Out.str("");
    Out << "DetectBleedTrails flagged " << total_trail_pixels << " pixel" 
	<< (bleed_status != 1 ? "s." : ".") << std::endl;
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  

  //  // TEMPORARY - EXAMINE THE INTERMEDIATE MASK
  //  for(Morph::IndexType p = 0;p < npix;p++)
  //    Inimage.DES()->mask[p] = temp_mask[p];

  // Now the trail pixels are marked, and some trail pixels are marked
  // with BADPIX_STAR as well.  If enabled, do "star interpolation" 
  // (i.e. radial interpolation) for the pixels marked with BADPIX_STAR, 
  // instead of the run-of-the-mill linear interpolation in the row 
  // for the normal bleed trail pixels.
  
  // Step 1 - Detect Star Blobs
  // Get the "blobs" of BADPIX_STAR pixels. The number of blobs is the 
  // new number of saturated objects after bleedtrail rejection.
  //
  // - blob_image is an integer image wherein the value indicates the
  //   blob to which the pixel belongs.
  // - blobs contains the pixel indices for each blob
  //
  blob_image.resize(npix,0);
  saturated_blobs.resize(0);
  // This call populates the above two data structures
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,BADPIX_STAR,blob_image,saturated_blobs);

  if(flag_verbose){
    Out.str("");
    Out << "Found " << saturated_blobs.size() << " saturated stars."; 
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  
  // profiler.FunctionEntry("RejectSmall");
  // First pass rejects any trail of size <ntrail_reject> or smaller
  std::vector<Morph::BoxType>  trail_boxes;
  if(ntrail_reject > 0){
    blob_image.resize(npix,0);
    trail_blobs.resize(0);
    Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,trail_mask,blob_image,trail_blobs);
    bbbi   = trail_blobs.begin();
    while(bbbi != trail_blobs.end()){
      int blobno = bbbi - trail_blobs.begin() + 1;
      Morph::BlobType &blob = *bbbi++;
      int nblobpix = blob.size();
      if(nblobpix <= ntrail_reject){
	// Need to correct image mask to remove trail
	Morph::BlobType::iterator blobit = blob.begin();
	while(blobit != blob.end()){
	  Morph::IndexType pixind = *blobit++;
	  (Inimage.DES()->mask)[pixind] &= ~(trail_mask|BADPIX_INTERP);
	  temp_image[pixind] = Inimage.DES()->image[pixind];
	}
	continue;
      }
    }
  }
  // profiler.FunctionExit("RejectSmall");
  
  // profiler.FunctionEntry("EdgeBleed");
  // Second pass blobs the valid trails and checks
  // for potential edgebleed
  //
  // First get the trail blobs
  bool bl_bleed_trail = false;
  bool br_bleed_trail = false;
  bool tl_bleed_trail = false;
  bool tr_bleed_trail = false;

  blob_image.resize(npix,0);
  trail_blobs.resize(0);					
  bool edge_bleed_detected = false;
  std::vector<Morph::IndexType> edge_bleed_blob_indices;
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,trail_mask,blob_image,trail_blobs);
  // Loop through the trail blobs and 
  bbbi   = trail_blobs.begin();
  while(bbbi != trail_blobs.end()){
    int blobno = bbbi - trail_blobs.begin();
    Morph::BlobType &blob = *bbbi++;
    // Step 1 - Get blob bounding box
    std::vector<Morph::IndexType> box;
    std::sort(blob.begin(),blob.end());
    Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
    //  adds the box to our list of blob boxes
    trail_boxes.push_back(box);
    // if the box is large enough and comes near the edge, edgebleed is possible
    Morph::IndexType yside = box[3] - box[2];
    double xpos = static_cast<double>(box[1] + box[0])/2.0;
    bool by_side = ((xpos < (10+npix_edge)) || (std::abs(xpos-ampx) < 10) || ((Nx-xpos) < (10+npix_edge)));
    if(yside > 2000 || by_side){
      if(((box[2] <= npix_edge) && read_register_on_bottom) ||
	 ((box[3] >= (Ny - npix_edge-1)) && read_register_on_top)){
	edge_bleed_detected = true;
	if((box[2] <= npix_edge)){
	  if(box[0] < ampx){
	    bl_bleed_trail = true;
	    if(bl_bleed) // both possible and detected earlier
	      edge_bleed_blob_indices.push_back(blobno);
	  }
	  if(box[1] >= ampx){
	    br_bleed_trail = true;
	    if(br_bleed) // both possible and detected earlier
	      edge_bleed_blob_indices.push_back(blobno);
	  }
	}
	if (box[3] >= (Ny - npix_edge - 1)){
	  if(box[0] < ampx){
	    tl_bleed_trail = true;
	    if(tl_bleed) // both possible and detected earlier
	      edge_bleed_blob_indices.push_back(blobno);
	  }	    
	  if(box[1] >= ampx){
	    tr_bleed_trail = true;
	    if(tr_bleed) // both possible and detected earlier
	      edge_bleed_blob_indices.push_back(blobno);
	  }
	}
      }
    }
  }
  std::vector<Morph::IndexType>::iterator ebbii = edge_bleed_blob_indices.begin();
  while(ebbii != edge_bleed_blob_indices.end()){
    Morph::BlobType &edgebleedblob(trail_blobs[*ebbii++]);
    Morph::BlobType::iterator ebbi = edgebleedblob.begin();
    while(ebbi != edgebleedblob.end()){
      Morph::IndexType blob_pixel_index = *ebbi++;
      Morph::IndexType blob_pixel_x = blob_pixel_index%Nx;
      Morph::IndexType blob_pixel_y = blob_pixel_index/Nx;
      blob_pixel_x = Nx - blob_pixel_x - 1;
      blob_pixel_index = blob_pixel_y*Nx + blob_pixel_x;
      Inimage.DES()->mask[blob_pixel_index] |= BADPIX_TRAIL;
    }
  }
  blob_image.resize(0);
  trail_blobs.resize(0);
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,trail_mask,blob_image,trail_blobs);


  // Since edge bleeds have already been detected, this code just marks them
  // in the output mask.
  if(edge_bleed_detected && ((bl_bleed && bl_bleed_trail) || (br_bleed && br_bleed_trail) ||
			     (tl_bleed && tl_bleed_trail) || (tr_bleed && tr_bleed_trail))) {
    Morph::MaskDataType pixel_rejection_mask = BADPIX_SATURATE | BADPIX_CRAY | BADPIX_BPM | BADPIX_TRAIL | BADPIX_STAR;
    Out.str("");
    Out << "Warning: Edgebleed confirmed, masking out affected edges:"
	<< std::endl;
    std::vector<Morph::IndexType> box(4,0);
    if(bl_bleed && bl_bleed_trail){
      box[0] = 0;
      box[1] = ampx-1;
      box[2] = 0;
      box[3] = y00;
//       Morph::IndexType grow_row = box[3]+1;
//       int npix_grow = 1;
//       int nrows_grown = 0;
//       while(npix_grow > 0){
// 	npix_grow = 0;
// 	for(Morph::IndexType col = box[0]; col <= box[1];col++){
// 	  if(!(Inimage.DES()->mask[grow_row*Nx+col]&pixel_rejection_mask) &&
// 	     (temp_image[grow_row*Nx+col] < hole_level))
// 	    npix_grow++;
// 	}
// 	if(npix_grow > 10){
// 	  grow_row++;
// 	  nrows_grown++;
// 	}
// 	else
// 	  npix_grow = 0;
//       }
//       Out << "Grew edgebleed by " << nrows_grown << "(" << hole_level << ")." << std::endl;
//       box[3] += nrows_grown;

      int bleedsize = box[3] - box[2];
      if(bleedsize > 400){
	Out << "  Suspect overflagging with edgebleed of "
	    << bleedsize << " pixels wide."
	    << std::endl;
	Morph::StatType boxstats;
	Morph::IndexType npix_box = 0;
	Morph::IndexType niter = 0;
	if(Morph::GetSkyBox(&temp_image[0],Inimage.DES()->mask, box,Nx, Ny,300,
			    nbgit,ground_rejection_factor,1e-3,pixel_rejection_mask,0,
			    boxstats,npix_box,niter,NULL)){
	  Out << "Could not get good statistics in edgebleed region, cannot attempt correction."
	      << std::endl;
	}
        else if(npix_box > 300){								
	  Out << "Stats in edgebleed (sky,sigma) = (" << boxstats[Image::IMMEAN] << ","
	      << boxstats[Image::IMSIGMA] << ")" << std::endl;
	  double boxhole_level = boxstats[Image::IMMEAN] - 10.0*boxstats[Image::IMSIGMA];
	  Out << "Detection level for edgebleed = " << boxhole_level << std::endl;
	  Morph::IndexType nrow = box[3] - box[2] + 1;
	  std::vector<int> npix_low_row(nrow,0);
	  for(Morph::IndexType row = 0;row < nrow;row++){
	    Morph::IndexType base_index = (row+box[2])*Nx;
	    for(Morph::IndexType col = box[0];col <= box[1];col++){
	      Morph::IndexType image_index = base_index+col;
	      if(!(Inimage.DES()->mask[image_index]&pixel_rejection_mask)){
		if(temp_image[image_index] < boxhole_level)
		  npix_low_row[row]++;
	      }
	    }
	  }
	  std::vector<int>::iterator nlri = npix_low_row.begin();
	  Morph::IndexType min_row = 5000;
	  Morph::IndexType max_row = 0;
	  for(Morph::IndexType row = box[2];row <= box[3];row++){
	    Morph::IndexType npixlow = *nlri++;
	    if(npixlow > 10){
	      if(row < min_row) min_row = row;
	      if(row > max_row) max_row = row;
	    }
	  }
	  box[3] = max_row;
	  int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/10.0);
	  overflag++;
	  box[3] += overflag;
	} else {
	  Out << "Not enough valid pixels in edgebleed region to attempt correction." << std::endl;
	}
      } else {
	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
	overflag++;
	box[3] += overflag;
      }
//       if(bleedsize > 400){
// 	Out << "  Suspect overflagging with edgebleed of "
// 	    << bleedsize << " pixels wide."
// 	    << std::endl;
//       } else {
// 	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
// 	box[3] += overflag;
//       }
      Out << "  Edgebleed(" << box[0] << ":" << box[1]
	  << "," << box[2] << ":" << box[3] << ")" << std::endl;
      for(Morph::IndexType by = box[2];by <= box[3];by++){
	Morph::IndexType iminy = by*Nx;
	for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	  Morph::IndexType imin = iminy + bx;
	  Inimage.DES()->mask[imin] |= BADPIX_TRAIL;
	}
      } 
    }
    if(tl_bleed && tl_bleed_trail){
      box[0] = 0;
      box[1] = ampx-1;
      box[2] = y01;
      box[3] = Ny-1;

//       Morph::IndexType grow_row = box[2]-1;
//       int npix_grow = 1;
//       int nrows_grown = 0;
//       while(npix_grow > 0){
// 	npix_grow = 0;
// 	for(Morph::IndexType col = box[0]; col <= box[1];col++){
// 	  if(!(Inimage.DES()->mask[grow_row*Nx+col]&pixel_rejection_mask) &&
// 	     (temp_image[grow_row*Nx+col] < hole_level))
// 	    npix_grow++;
// 	}
// 	if(npix_grow > 10){
// 	  grow_row--;
// 	  nrows_grown++;
// 	}
// 	else
// 	  npix_grow = 0;
//       }
//       Out << "Grew edgebleed by " << nrows_grown << "(" << hole_level << ")." << std::endl;
//       box[2] -= nrows_grown;

      int bleedsize = box[3] - box[2];
      if(bleedsize > 400){
	Out << "  Suspect overflagging with edgebleed of "
	    << bleedsize << " pixels wide."
	    << std::endl;
	Morph::StatType boxstats;
	Morph::IndexType npix_box = 0;
	Morph::IndexType niter = 0;
	if(Morph::GetSkyBox(&temp_image[0],Inimage.DES()->mask, box,Nx, Ny,300,
			    nbgit,ground_rejection_factor,1e-3,pixel_rejection_mask,0,
			    boxstats,npix_box,niter,NULL)){
	  Out << "Could not get good statistics in edgebleed region, cannot attempt correction."
	      << std::endl;
	}
        else if(npix_box > 300){								
	  Out << "Stats in edgebleed (sky,sigma) = (" << boxstats[Image::IMMEAN] << ","
	      << boxstats[Image::IMSIGMA] << ")" << std::endl;
	  double boxhole_level = boxstats[Image::IMMEAN] - 10.0*boxstats[Image::IMSIGMA];
	  Out << "Detection level for edgebleed = " << boxhole_level << std::endl;
	  Morph::IndexType nrow = box[3] - box[2] + 1;
	  std::vector<int> npix_low_row(nrow,0);
	  for(Morph::IndexType row = 0;row < nrow;row++){
	    Morph::IndexType base_index = (row+box[2])*Nx;
	    for(Morph::IndexType col = box[0];col <= box[1];col++){
	      Morph::IndexType image_index = base_index+col;
	      if(!(Inimage.DES()->mask[image_index]&pixel_rejection_mask)){
		if(temp_image[image_index] < boxhole_level)
		  npix_low_row[row]++;
	      }
	    }
	  }
	  std::vector<int>::iterator nlri = npix_low_row.begin();
	  Morph::IndexType min_row = 5000;
	  Morph::IndexType max_row = 0;
	  for(Morph::IndexType row = box[2];row <= box[3];row++){
	    Morph::IndexType npixlow = *nlri++;
	    if(npixlow > 10){
	      if(row < min_row) min_row = row;
	      if(row > max_row) max_row = row;
	    }
	  }
	  box[2] = min_row;
	  int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/10.0);
	  overflag++;
	  box[2] -= overflag;
	} else {
	  Out << "Not enough valid pixels in edgebleed region to attempt correction." << std::endl;
	}
      } else {
	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
	overflag++;
	box[2] -= overflag;
      }
      Out << "  Edgebleed(" << box[0] << ":" << box[1]
	  << "," << box[2] << ":" << box[3] << ")" << std::endl;
      for(Morph::IndexType by = box[2];by <= box[3];by++){
	Morph::IndexType iminy = by*Nx;
	for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	  Morph::IndexType imin = iminy + bx;
	  Inimage.DES()->mask[imin] |= BADPIX_TRAIL;
	}
      } 
    }
    if(br_bleed && br_bleed_trail){
      box[0] = ampx;
      box[1] = Nx-1;
      box[2] = 0;
      box[3] = y10;

//       Morph::IndexType grow_row = box[3]+1;
//       int npix_grow = 1;
//       int nrows_grown = 0;
//       while(npix_grow > 0){
// 	npix_grow = 0;
// 	for(Morph::IndexType col = box[0]; col <= box[1];col++){
// 	  if(!(Inimage.DES()->mask[grow_row*Nx+col]&pixel_rejection_mask) &&
// 	     (temp_image[grow_row*Nx+col] < hole_level))
// 	    npix_grow++;
// 	}
// 	if(npix_grow > 10){
// 	  grow_row++;
// 	  nrows_grown++;
// 	}
// 	else
// 	  npix_grow = 0;
//       }
//       Out << "Grew edgebleed by " << nrows_grown << "(" << hole_level << ")." << std::endl;
//       box[3] += nrows_grown;

      int bleedsize = box[3] - box[2];
      if(bleedsize > 400){
	Out << "  Suspect overflagging with edgebleed of "
	    << bleedsize << " pixels wide."
	    << std::endl;
	Morph::StatType boxstats;
	Morph::IndexType npix_box = 0;
	Morph::IndexType niter = 0;
	if(Morph::GetSkyBox(&temp_image[0],Inimage.DES()->mask, box,Nx, Ny,300,
			    nbgit,ground_rejection_factor,1e-3,pixel_rejection_mask,0,
			    boxstats,npix_box,niter,NULL)){
	  Out << "Could not get good statistics in edgebleed region, cannot attempt correction."
	      << std::endl;
	}
        else if(npix_box > 300){								
	  Out << "Stats in edgebleed (sky,sigma) = (" << boxstats[Image::IMMEAN] << ","
	      << boxstats[Image::IMSIGMA] << ")" << std::endl;
	  double boxhole_level = boxstats[Image::IMMEAN] - 10.0*boxstats[Image::IMSIGMA];
	  Out << "Detection level for edgebleed = " << boxhole_level << std::endl;
	  Morph::IndexType nrow = box[3] - box[2] + 1;
	  std::vector<int> npix_low_row(nrow,0);
	  for(Morph::IndexType row = 0;row < nrow;row++){
	    Morph::IndexType base_index = (row+box[2])*Nx;
	    for(Morph::IndexType col = box[0];col <= box[1];col++){
	      Morph::IndexType image_index = base_index+col;
	      if(!(Inimage.DES()->mask[image_index]&pixel_rejection_mask)){
		if(temp_image[image_index] < boxhole_level)
		  npix_low_row[row]++;
	      }
	    }
	  }
	  std::vector<int>::iterator nlri = npix_low_row.begin();
	  Morph::IndexType min_row = 5000;
	  Morph::IndexType max_row = 0;
	  for(Morph::IndexType row = box[2];row <= box[3];row++){
	    Morph::IndexType npixlow = *nlri++;
	    if(npixlow > 10){
	      if(row < min_row) min_row = row;
	      if(row > max_row) max_row = row;
	    }
	  }
	  box[3] = max_row;
	  int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/10.0);
	  overflag++;
	  box[3] += overflag;
	} else {
	  Out << "Not enough valid pixels in edgebleed region to attempt correction." << std::endl;
	}
      } else {
	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
	overflag++;
	box[3] += overflag;
      }
//       if(bleedsize > 400){
// 	Out << "  Suspect overflagging with edgebleed of "
// 	    << bleedsize << " pixels wide."
// 	    << std::endl;
//       } else {
// 	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
// 	box[3] += overflag;
//       }
      Out << "  Edgebleed(" << box[0] << ":" << box[1]
	  << "," << box[2] << ":" << box[3] << ")" << std::endl;      
      for(Morph::IndexType by = box[2];by <= box[3];by++){
	Morph::IndexType iminy = by*Nx;
	for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	  Morph::IndexType imin = iminy + bx;
	  Inimage.DES()->mask[imin] |= BADPIX_TRAIL;
	}
      } 
    }
    if(tr_bleed && tr_bleed_trail){
      box[0] = ampx;
      box[1] = Nx-1;
      box[2] = y11;
      box[3] = Ny-1;
//       Morph::IndexType grow_row = box[2]-1;
//       int npix_grow = 1;
//       int nrows_grown = 0;
//       while(npix_grow > 0){
// 	npix_grow = 0;
// 	for(Morph::IndexType col = box[0]; col <= box[1];col++){
// 	  if(!(Inimage.DES()->mask[grow_row*Nx+col]&pixel_rejection_mask) &&
// 	     (temp_image[grow_row*Nx+col] < hole_level))
// 	    npix_grow++;
// 	}
// 	if(npix_grow > 10){
// 	  grow_row--;
// 	  nrows_grown++;
// 	}
// 	else
// 	  npix_grow = 0;
//       }
//       Out << "Grew edgebleed by " << nrows_grown << "(" << hole_level << ")." << std::endl;
//       box[2] -= nrows_grown;
      int bleedsize = box[3] - box[2];
      if(bleedsize > 400){
	Out << "  Suspect overflagging with edgebleed of "
	    << bleedsize << " pixels wide."
	    << std::endl;
	Morph::StatType boxstats;
	Morph::IndexType npix_box = 0;
	Morph::IndexType niter = 0;
	if(Morph::GetSkyBox(&temp_image[0],Inimage.DES()->mask, box,Nx, Ny,300,
			    nbgit,ground_rejection_factor,1e-3,pixel_rejection_mask,0,
			    boxstats,npix_box,niter,NULL)){
	  Out << "Could not get good statistics in edgebleed region, cannot attempt correction."
	      << std::endl;
	}
        else if(npix_box > 300){								
	  Out << "Stats in edgebleed (sky,sigma) = (" << boxstats[Image::IMMEAN] << ","
	      << boxstats[Image::IMSIGMA] << ")" << std::endl;
	  double boxhole_level = boxstats[Image::IMMEAN] - 10.0*boxstats[Image::IMSIGMA];
	  Out << "Detection level for edgebleed = " << boxhole_level << std::endl;
	  Morph::IndexType nrow = box[3] - box[2] + 1;
	  std::vector<int> npix_low_row(nrow,0);
	  for(Morph::IndexType row = 0;row < nrow;row++){
	    Morph::IndexType base_index = (row+box[2])*Nx;
	    for(Morph::IndexType col = box[0];col <= box[1];col++){
	      Morph::IndexType image_index = base_index+col;
	      if(!(Inimage.DES()->mask[image_index]&pixel_rejection_mask)){
		if(temp_image[image_index] < boxhole_level)
		  npix_low_row[row]++;
	      }
	    }
	  }
	  std::vector<int>::iterator nlri = npix_low_row.begin();
	  Morph::IndexType min_row = 5000;
	  Morph::IndexType max_row = 0;
	  for(Morph::IndexType row = box[2];row <= box[3];row++){
	    Morph::IndexType npixlow = *nlri++;
	    if(npixlow > 10){
	      if(row < min_row) min_row = row;
	      if(row > max_row) max_row = row;
	    }
	  }
	  box[2] = min_row;
	  int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/10.0);
	  overflag++;
	  box[2] -= overflag;
	} else {
	  Out << "Not enough valid pixels in edgebleed region to attempt correction." << std::endl;
	}
      } else {
	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
	overflag++;
	box[2] -= overflag;
      }
//       if(bleedsize > 400){
// 	Out << "  Suspect overflagging with edgebleed of "
// 	    << bleedsize << " pixels wide."
// 	    << std::endl;
//       } else {
// 	int overflag = static_cast<int>(static_cast<double>(box[3] - box[2])/8.0);
// 	box[2] -= overflag;
//       }
      Out << "  Edgebleed(" << box[0] << ":" << box[1]
	  << "," << box[2] << ":" << box[3] << ")" << std::endl;
      for(Morph::IndexType by = box[2];by <= box[3];by++){
	Morph::IndexType iminy = by*Nx;
	for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
	  Morph::IndexType imin = iminy + bx;
	  Inimage.DES()->mask[imin] |= BADPIX_TRAIL;
	}
      } 
    }  
    LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
  }
  // profiler.FunctionExit("EdgeBleed");
  
  
  // profiler.FunctionEntry("Stars");
  std::vector<Morph::BoxType>  star_boxes;
  std::vector<double> star_centers_x;
  std::vector<double> star_centers_y;
  std::vector<double> star_radii;
  std::vector<std::vector<double> > star_idf;
  // Loop through each BADPIX_STAR blob, get blob bounding box,
  // star center, and pixel intensity distribution function
  //
  // Step 0 - Initiate loop over blobs
  std::ofstream BlobFile;
  std::vector<Morph::BlobType>::iterator bbi   = saturated_blobs.begin();
  std::vector<int> nstars_pix(npix,0);
  while(bbi != saturated_blobs.end()){
    int blobno = bbi - saturated_blobs.begin() + 1;
    //      std::ostringstream BFOut;
    //      BFOut << "blob_" << blobno;
    //      BlobFile.open(BFOut.str().c_str());
    
    Morph::BlobType &blob = *bbi++;
    int nblobpix = blob.size();


    Morph::BoxType box;
    std::sort(blob.begin(),blob.end());
    Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
    star_boxes.push_back(box);
    
    Morph::IndexType nby = box[3] - box[2] + 1;
    Morph::IndexType nbx = box[1] - box[0] + 1;
    double cx = static_cast<double>(nbx);
    double cy = static_cast<double>(nby);
    cx /= 2.0;
    cy /= 2.0;
    double star_r = cx*cx;
    if(debug)
      std::cout << "Initial Stellar Radius(" << blobno << ") = " << std::sqrt(star_r) << std::endl;
    if(cy < cx) star_r = cy*cy;
    cx += (box[0]);
    cy += (box[2]);
    cx -= .5;
    cy -= .5;
    star_centers_x.push_back(cx);
    star_centers_y.push_back(cy);
    
    // **** NOTE: This parameter controls the star radius detection.  It
    // is the number of sigma above which star levels are detected.
    Morph::ImageDataType star_scalefactor2 = 3.0; 

    // profiler.FunctionEntry("StarRadius");
    ModifyStarR(Inimage.DES()->mask,Inimage.DES()->image,cx,cy,star_r,Nx,Ny,
		ground_rejection_mask,star_scalefactor2,image_stats);
    // profiler.FunctionExit("StarRadius");
    star_r *= rgf2;
    double srstarr = std::sqrt(star_r);
    if(debug)
      std::cout << "Modified Stellar Radius(" << blobno << "): " << srstarr << std::endl;

    star_radii.push_back(srstarr);

    if(do_star_interp){
      // profiler.FunctionEntry("StarInterp");
      long r = static_cast<long>((std::sqrt(star_r)) + 2.0);
      std::vector<double> idf(r,0);

      long center_x1 = static_cast<long>(cx);
      long center_y1 = static_cast<long>(cy);
      long center_x2 = static_cast<long>(cx+.5);
      long center_y2 = static_cast<long>(cy+.5);
      long center_y = center_y1;
      long npix_y = 0;
      if(center_y1 != center_y2){
	if((center_y2 < Ny) && (center_y2 >= 0)){
	  center_y = center_y2;
	  npix_y++;
	}
	if((center_y1 < Ny) && (center_y1 >= 0)){
	  npix_y++;
	  center_y = center_y1;
	}
      }
      else {
	if(center_y1 >= Ny) center_y1 = Ny;
	if(center_y1 < 0) center_y1 = 0;
	npix_y = 1;
      }

      std::vector<double>::iterator idfit = idf.begin();
      long pixx = 0;
      long npix = 0;
      while(idfit != idf.end()){
	Morph::IndexType vindex = idfit++ - idf.begin();
	if((center_x1 + pixx) < Nx){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x1 + pixx;
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if((center_x1 - pixx) >= 0){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x1 - pixx;
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if((center_x2 + pixx) < Nx){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x2 + pixx;
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if((center_x2 - pixx) >= 0){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x2 - pixx;
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if(npix > 0){
	  idf[vindex] /= static_cast<double>(npix);
	}
	npix = 0;
	pixx++;
      }
      //       idfit = idf.begin();
      //       while(idfit != idf.end()){
      // 	idfit++;
      //       }
      
      // Loop over blob pixels and set values according to linearly
      // interpolated values from the idf.
      Morph::BlobType::iterator sbi = blob.begin();
      while(sbi != blob.end()){
	Morph::IndexType pixel_index = *sbi++;
	if(Inimage.DES()->mask[pixel_index]&BADPIX_TRAIL){
	  Morph::IndexType spixy = pixel_index/Nx;
	  Morph::IndexType spixx = pixel_index%Nx;
	  double x_distance = static_cast<double>(spixx) - cx;
	  double y_distance = static_cast<double>(spixy) - cy;
	  double distance = std::sqrt(x_distance*x_distance + y_distance*y_distance);
	  Morph::IndexType dist_index = static_cast<Morph::IndexType>(distance);
	  Morph::IndexType other_dist_index = static_cast<Morph::IndexType>(distance+.5);
	  if(dist_index == other_dist_index)
	    temp_image[pixel_index] = idf[dist_index];
	  else{
	    double weight_1 = 1.0 - std::abs(distance - dist_index);
	    double weight_2 = 1.0 - std::abs(distance - other_dist_index);
	    if(weight_1 < 1e-5){
	      temp_image[pixel_index] = idf[other_dist_index];
	    }
	    else if(weight_2 < 1e-5){
	      temp_image[pixel_index] = idf[dist_index];
	    }
	    else
	      temp_image[pixel_index] = weight_1*idf[dist_index] + 
		weight_2*idf[other_dist_index];
	  }
	}
	if(!do_starmask)
	  Inimage.DES()->mask[pixel_index] ^= BADPIX_STAR;
      }
      // profiler.FunctionExit("StarInterp");
    }
    if(do_starmask){
      // profiler.FunctionEntry("StarMask");
      double sx1 = cx - static_cast<double>(nby)/2.0 - 1.0;
      double sx2 = cx + static_cast<double>(nby)/2.0 + 1.0;
      box[0] = static_cast<Morph::IndexType>(sx1);
      box[1] = static_cast<Morph::IndexType>(sx2);
      // Ensure boxes do not step outside image, and that they
      // exceed the detected star_r.
      double nxstar = (box[1] - box[0] + 1)/2.0;
      if(nxstar < srstarr){
	box[0] = static_cast<Morph::IndexType>(cx - srstarr - 4.0);
	box[1] = static_cast<Morph::IndexType>(cx + srstarr + 4.0);
	box[2] = static_cast<Morph::IndexType>(cy - srstarr - 4.0);
	box[3] = static_cast<Morph::IndexType>(cy + srstarr + 4.0);
      }
      if(box[0] < 0) box[0] = 0;
      if(box[1] >= Nx) box[1] = Nx-1;
      if(box[2] < 0) box[2] = 0;
      if(box[3] >= Ny) box[3] = Ny-1;
      for(Morph::IndexType byy=box[2];byy <= box[3];byy++){
	Morph::IndexType byindex = byy*Nx;
	for(Morph::IndexType bxx=box[0];bxx <= box[1];bxx++){
	  Morph::IndexType bindex = byindex + bxx;
	  double x_distance = static_cast<double>(bxx) - cx;
	  double y_distance = static_cast<double>(byy) - cy;
	  double r2 = x_distance*x_distance + y_distance*y_distance;
	  if(r2 <= star_r){
	    Inimage.DES()->mask[bindex] |= BADPIX_STAR;
	    nstars_pix[bindex]++;
	  }
	  else if((Inimage.DES()->mask[bindex]&BADPIX_STAR) && !nstars_pix[bindex])
	    Inimage.DES()->mask[bindex] ^= BADPIX_STAR;
	}
      }
      // profiler.FunctionExit("StarMask");
    } 
  }

  // profiler.FunctionEntry("RejectObscured");
  // Reject stars obscured by bleeds 
  blob_image.resize(npix,0);
  saturated_blobs.resize(0);
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,BADPIX_STAR,blob_image,saturated_blobs);
  std::vector<Morph::BlobType>::iterator blobit = saturated_blobs.begin();
  std::vector<double> new_center_x;
  std::vector<double> new_center_y;
  std::vector<double> new_radius;
  std::vector<double>::iterator cxi  = star_centers_x.begin();
  std::vector<double>::iterator cyi  = star_centers_y.begin();
  std::vector<double>::iterator sri  = star_radii.begin();
  Morph::IndexType ncovered = 0;
  while(blobit != saturated_blobs.end()){
    double center_x = *cxi++;
    double center_y = *cyi++;
    double srad = *sri++;
    Morph::BlobType &blob = *blobit++;
    Morph::BlobType::iterator bi = blob.begin();
    bool covered_by_trail = true;
    while(bi != blob.end() && covered_by_trail){
      Morph::IndexType pixel_index = *bi++;
      if(!(Inimage.DES()->mask[pixel_index] & BADPIX_TRAIL))
	covered_by_trail = false;
    }
    if(covered_by_trail){
      ncovered++;
      bi = blob.begin();
      while(bi != blob.end() && covered_by_trail){
	Morph::IndexType pixel_index = *bi++;
	Inimage.DES()->mask[pixel_index] &= ~BADPIX_STAR;
      }
    }
    else{
      new_center_x.push_back(center_x);
      new_center_y.push_back(center_y);
      new_radius.push_back(srad);
    }
  }
  star_centers_x.resize(0);
  star_centers_y.resize(0);
  star_radii.resize(0);
  cxi = new_center_x.begin();
  cyi = new_center_y.begin();
  sri = new_radius.begin();
  while(sri != new_radius.end()){
    star_centers_x.push_back(*cxi++);
    star_centers_y.push_back(*cyi++);
    star_radii.push_back(*sri++);
  }
  if(ncovered > 0){
    Out.str("");
    Out << "Rejected " << ncovered << " stars due to being completely obscured by bleeds.";
    LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
  }
  // profiler.FunctionExit("RejectObscured");
  // profiler.FunctionExit("Stars");
  // profiler.FunctionEntry("Output");
  // profiler.FunctionEntry("WCS");
  // Check for a WCS soln for this image
  bool get_wcs               = true;
  bool write_radii_in_pixels = false;
  double convert_to_arcsecs  = 1.0;
  wcsstruct *wcs_in          = NULL;
  catstruct *cat             = NULL;
  tabstruct *tab             = NULL;
  std::vector<double> imjacobian(4,0);
  if(FitsTools::HeaderKeyExists(Inimage.ImageHeader(),"SCAMPFLG")){
    if(FitsTools::GetHeaderValue<int>(Inimage.ImageHeader(),"SCAMPFLG") > 0) get_wcs = false;
  } else if(!FitsTools::HeaderKeyExists(Inimage.ImageHeader(),"SCAMPFLG")){
    get_wcs = false;
    Out.str("");
    Out << "Failed to find SCAMPFLG keyword in header of " << Inimage.DES()->name 
	<< ". Writing tables in image coordinates.";
    LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
  }
  if(get_wcs) {
    if (!(cat = read_cat(Inimage.DES()->name))){
      Out.str("");
      Out << "Failed to read catalog from " << Inimage.DES()->name 
	  << ". Writing tables in image coordinates.";
      LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
      get_wcs = false;
    } else {
      tab=cat->tab;
      wcs_in=read_wcs(tab);
      // Get coordinate transformation (if it exists)
      if(FitsTools::HeaderKeyExists(Inimage.ImageHeader(),"CD1_1")){
	imjacobian[0] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD1_1");
	imjacobian[1] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD1_2");
	imjacobian[2] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD2_1");
	imjacobian[3] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD2_2");
	convert_to_arcsecs = std::sqrt(imjacobian[0]*imjacobian[0] + imjacobian[1]*imjacobian[1]);
	convert_to_arcsecs += std::sqrt(imjacobian[2]*imjacobian[2] + imjacobian[3]*imjacobian[3]);
	convert_to_arcsecs *= 1800.0;
      } else {
	Out.str("");
	Out << "Failed to read coordinate transformation from " << Inimage.DES()->name 
	    << ". Writing object radii in pixels." << std::endl;
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	write_radii_in_pixels = true;
      }
    }
  }
  // profiler.FunctionExit("WCS");
  // profiler.FunctionEntry("ObjectTable");
  // Write the bright object table.
  if(!SaturatedObjectFileName.empty() && !star_centers_x.empty()){
    std::vector<std::string> field_names;
    std::vector<std::string> field_units(3,"pixels");
    std::vector<std::string> field_types(3,"D");
    field_names.push_back("X_IMAGE");
    field_names.push_back("Y_IMAGE");
    field_names.push_back("Radius");
    if(get_wcs){
      // replace field names if wcs will be used
      field_names[0].assign("RA");
      field_names[1].assign("DEC");
      field_units[0].assign("Degrees");
      field_units[1].assign("Degrees");
      if(!write_radii_in_pixels)
	field_units[2].assign("Arcsecs");
    }
    // Copy either the original, or the WCS converted values into
    // temporary arrays for table building
    double rawpos[2],wcspos[2];
    std::vector<double> star_centers_wcs_x(star_centers_x.size(),0);
    std::vector<double> star_centers_wcs_y(star_centers_y.size(),0);
    std::vector<double> star_radii_wcs(star_centers_x.size(),0);
    std::vector<double>::iterator wcsx_i  = star_centers_wcs_x.begin();
    std::vector<double>::iterator wcsy_i  = star_centers_wcs_y.begin();
    std::vector<double>::iterator wcsr_i  = star_radii_wcs.begin();
    std::vector<double>::iterator scx_i   = star_centers_x.begin();
    std::vector<double>::iterator scy_i   = star_centers_y.begin();
    std::vector<double>::iterator srad_i  = star_radii.begin();
    while(wcsx_i != star_centers_wcs_x.end()){
      if(get_wcs){
	rawpos[0] = *scx_i+1.0;
	rawpos[1] = *scy_i+1.0;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	*wcsx_i = wcspos[0];
	*wcsy_i = wcspos[1];
      } else {
	*wcsx_i = *scx_i+1.0;
	*wcsy_i = *scy_i+1.0;
      }
      *wcsr_i++ = *srad_i++*convert_to_arcsecs;
      wcsx_i++;
      wcsy_i++;
      scx_i++;
      scy_i++;
    }
    // Create the object table
    FitsTools::FitsImage StarTableOut;
    StarTableOut.SetOutStream(Out);
    StarTableOut.CreateFile(SaturatedObjectFileName,true,flag_verbose);
    StarTableOut.MakeTable(3,field_names,field_types,field_units,"TABLE",flag_verbose);

    // write center coordinate 1 to object table.
    if(StarTableOut.WriteTableColumn(1,star_centers_x.size(),TDOUBLE,
				     (void *)&star_centers_wcs_x[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "RA " : "x coordindate ")
	  << " of star centers in " << SaturatedObjectFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    // write center coordinate 2 to object table.
    if(StarTableOut.WriteTableColumn(2,star_centers_y.size(),TDOUBLE,
				     (void *)&star_centers_wcs_y[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "DEC " : "y coordindate ")
	  << " of star centers in " << SaturatedObjectFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
    // write radii to object table
    if(StarTableOut.WriteTableColumn(3,star_radii.size(),TDOUBLE,(void *)&star_radii_wcs[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert star radii in " 
	  << SaturatedObjectFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
    StarTableOut.Close();
  }
  // profiler.FunctionExit("ObjectTable");
  // profiler.FunctionEntry("TrailTable");

  // Make table for trail blob bounding boxes (if enabled)
  if(!TrailBoxesFileName.empty() && !trail_boxes.empty()){
    std::vector<std::string> field_names;
    std::vector<std::string> field_units(8,"Degrees");
    std::vector<std::string> field_types(8,"D");
    if(get_wcs){
      field_names.push_back("RA_1");
      field_names.push_back("RA_2");
      field_names.push_back("RA_3");
      field_names.push_back("RA_4");
      field_names.push_back("DEC_1");
      field_names.push_back("DEC_2");
      field_names.push_back("DEC_3");
      field_names.push_back("DEC_4");
    } else {
      field_names.push_back("X_1");
      field_names.push_back("X_2");
      field_names.push_back("X_3");
      field_names.push_back("X_4");
      field_names.push_back("Y_1");
      field_names.push_back("Y_2");
      field_names.push_back("Y_3");
      field_names.push_back("Y_4");
      field_units[0].assign("pixels");
      field_units[1].assign("pixels");
      field_units[2].assign("pixels");
      field_units[3].assign("pixels");
      field_units[4].assign("pixels");
      field_units[5].assign("pixels");
      field_units[6].assign("pixels");
      field_units[7].assign("pixels");
    }
    std::vector<double> LocX1;
    std::vector<double> LocX2;
    std::vector<double> LocX3;
    std::vector<double> LocX4;
    std::vector<double> LocY1;
    std::vector<double> LocY2;
    std::vector<double> LocY3;
    std::vector<double> LocY4;
    // Append the edgebleed boxes to the end of the trail boxes 
    std::vector<Morph::BoxType>::iterator ebi = edgebleed_boxes.begin();
    while(ebi != edgebleed_boxes.end())
      trail_boxes.push_back(*ebi++);
    
    std::vector<Morph::BoxType>::iterator tbi = trail_boxes.begin();
    while(tbi != trail_boxes.end()){
      unsigned int box_no = tbi-trail_boxes.begin() + 1;
      Morph::BoxType &box = *tbi++;
      if(get_wcs){
	double rawpos[2],wcspos[2];
	rawpos[0] = static_cast<double>(box[0])+.5;
	rawpos[1] = static_cast<double>(box[2])+.5;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX1.push_back(wcspos[0]);
	LocY1.push_back(wcspos[1]);

	rawpos[0] = static_cast<double>(box[1])+1.5;
	rawpos[1] = static_cast<double>(box[2])+.5;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX2.push_back(wcspos[0]);
	LocY2.push_back(wcspos[1]);

	rawpos[0] = static_cast<double>(box[1])+1.5;
	rawpos[1] = static_cast<double>(box[3])+1.5;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX3.push_back(wcspos[0]);
	LocY3.push_back(wcspos[1]);

	rawpos[0] = static_cast<double>(box[0])+.5;
	rawpos[1] = static_cast<double>(box[3])+1.5;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX4.push_back(wcspos[0]);
	LocY4.push_back(wcspos[1]);
      } else {
	LocX1.push_back(static_cast<double>(box[0])+.5);
	LocY1.push_back(static_cast<double>(box[2])+.5);
	LocX2.push_back(static_cast<double>(box[1])+1.5);
	LocY2.push_back(static_cast<double>(box[2])+.5);
	LocX3.push_back(static_cast<double>(box[1])+1.5);
	LocY3.push_back(static_cast<double>(box[3])+1.5);
	LocX4.push_back(static_cast<double>(box[0])+.5);
	LocY4.push_back(static_cast<double>(box[3])+1.5);
      }
    }
    // Write the FITS table for the blob boxes
    FitsTools::FitsImage TrailBoxesOut;
    TrailBoxesOut.SetOutStream(Out);
    TrailBoxesOut.CreateFile(TrailBoxesFileName,true,flag_verbose);
    TrailBoxesOut.MakeTable(8,field_names,field_types,field_units,"TABLE",flag_verbose);
    if(TrailBoxesOut.WriteTableColumn(1,trail_boxes.size(),TDOUBLE,(void *)&LocX1[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "RA " : "X coordinate ")
	  << " of trailbox corner 1 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(2,trail_boxes.size(),TDOUBLE,(void *)&LocX2[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "RA " : "X coordinate ")
	  << " of trailbox corner 2 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(3,trail_boxes.size(),TDOUBLE,(void *)&LocX3[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "RA " : "X coordinate ")
	  << " of trailbox corner 3 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(4,trail_boxes.size(),TDOUBLE,(void *)&LocX4[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "RA " : "X coordinate ")
	  << " of trailbox corner 4 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(5,trail_boxes.size(),TDOUBLE,(void *)&LocY1[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "DEC " : "Y coordinate ")
	  << " of trailbox corner 1 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(6,trail_boxes.size(),TDOUBLE,(void *)&LocY2[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "DEC " : "Y coordinate ")
	  << " of trailbox corner 2 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(7,trail_boxes.size(),TDOUBLE,(void *)&LocY3[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "DEC " : "Y coordinate ")
	  << " of trailbox corner 3 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    if(TrailBoxesOut.WriteTableColumn(8,trail_boxes.size(),TDOUBLE,(void *)&LocY4[0],flag_verbose)){
      Out << program_name << "::Error: Could not insert " << (get_wcs ? "DEC " : "Y coordinate ")
	  << " of trailbox corner 4 in " << TrailBoxesFileName;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
      
    }
    TrailBoxesOut.Close();
  }
  // profiler.FunctionExit("TrailTable");
  if(!do_starmask){
    for(int mindex = 0;mindex < npix;mindex++){
      if(Inimage.DES()->mask[mindex]&BADPIX_STAR)
	Inimage.DES()->mask[mindex] ^= BADPIX_STAR;
    }
  }
  if(do_interp || do_star_interp){
    std::vector<Morph::ImageDataType>::iterator tii = temp_image.begin();
    while(tii != temp_image.end()){
      Morph::IndexType indy = tii-temp_image.begin();
      Inimage.DES()->image[indy] = *tii++;
    }
  }
  if(Inimage.DES()->varim){
    for(int mindex = 0;mindex < npix;mindex++){
      if(Inimage.DES()->mask[mindex]&BADPIX_TRAIL)
	Inimage.DES()->varim[mindex] = 0.0;
    }
    if(do_zero){
      for(int mindex = 0;mindex < npix;mindex++){
	if(Inimage.DES()->mask[mindex]&BADPIX_STAR)
	  Inimage.DES()->varim[mindex] = 0.0;
      }
    }
  }

  UpdateImageHeader(Inimage,argv,bleed_status,do_interp,do_star_interp,
		    do_starmask,bleedmask_status);


  return(WriteOutputImage(Inimage,ofilename,flag_verbose));

}

// Starts with an initial guess at the star radius (star_r) and looks 
// in a disk of 5*star_r.  Integrates pixel energy over annular bins 
// and searches inward from disk edge until star level is detected.  
// Returns the resulting radius in star_r.
void ModifyStarR(Morph::MaskDataType *mask,Morph::ImageDataType *image,
		 double cx,double cy,double &star_r,Morph::IndexType Nx,
		 Morph::IndexType Ny,Morph::MaskDataType rejection_mask,
		 double star_scalefactor,Morph::StatType &stats)
{
  double starval = stats[Image::IMMEAN] + star_scalefactor*stats[Image::IMSIGMA];
  double temp_r = 10.0*std::sqrt(star_r);
  Morph::IndexType nbins = static_cast<Morph::IndexType>(temp_r + 1);  
  double binsize = temp_r/static_cast<double>(nbins);
  std::vector<std::vector<double> > starbins(nbins+1);
  Morph::IndexType x0 = static_cast<Morph::IndexType>(cx - temp_r - 1);
  Morph::IndexType x1 = static_cast<Morph::IndexType>(cx + temp_r + 1);
  Morph::IndexType y0 = static_cast<Morph::IndexType>(cy - temp_r - 1);
  Morph::IndexType y1 = static_cast<Morph::IndexType>(cy + temp_r + 1);
  if(x0 < 0) x0 = 0;
  if(x1 >= Nx) x1 = Nx-1;
  if(y0 < 0) y0 = 0;
  if(y1 >= Ny) y1 = Ny-1;
  for(Morph::IndexType pixx = x0;pixx <= x1;pixx++){
    for(Morph::IndexType pixy = y0;pixy <= y1;pixy++){
      Morph::IndexType index = pixy*Nx + pixx;
      if(!(mask[index]&rejection_mask)){
	double pixdist = std::sqrt((pixx-cx)*(pixx-cx) + (pixy-cy)*(pixy-cy));
	if(pixdist <= temp_r){
	  Morph::ImageDataType pixval = image[index];
	  Morph::IndexType bindex = static_cast<Morph::IndexType>(pixdist/binsize);
	  assert((bindex>=0) && (bindex <= nbins)); 
	  starbins[bindex].push_back(pixval);
	}
      }
    }
  }
  bool done = false;
  for(Morph::IndexType bindex = 0;((bindex < nbins) && !done);bindex++){
    double median = 0.0;
    if(!starbins[bindex].empty()){
      Morph::IndexType nval = starbins[bindex].size();
      std::sort(starbins[bindex].begin(),starbins[bindex].end());
      if(!(nval%2)){
	median = (starbins[bindex][nval/2]+starbins[bindex][(nval/2)-1])/2.0; 
      } else {
	median = starbins[bindex][nval/2];
      }
      if(median < starval){
	star_r = (bindex+1)*binsize;
	done = true;
      } 
    }
  }
  if(done)
    star_r *= star_r;
}


int UpdateImageHeader(FitsTools::FitsImage &Inimage,const char *argv[],
		      int nbleedpix,bool do_interp,bool do_star_interp,
		      bool do_starmask,int bleedmask_status)
{
  // Leave processing history in output header
  //
  
  // Get date and time
  time_t tm = time(NULL);
  std::string dateandtime(asctime(localtime(&tm)));
  
  // For some oddball reason, the dateandtime string appears to 
  // have a newline in it.  Work around this by reading in each
  // token as a separate string.
  std::istringstream DStr(dateandtime);
  std::string day,month,date,stime,year;
  DStr >> day >> month >> date >> stime >> year;
  
  // Form the processing messages and append them to the header
  std::ostringstream Hstr;
  Hstr << "DESBLEED= \'" << day << " " << month << " " << date 
       << " " << stime << " " << year 
       << (bleedmask_status ? "\' / bleed trail masking (FAILED)" : "\' / bleed trail masking");
  Inimage.AppendImageHeader(Hstr.str());
  Hstr.str("");
  Hstr << "NBLEED  = " << nbleedpix << "  / Number of bleed trail pixels.";
  Inimage.AppendImageHeader(Hstr.str());
  if(do_interp && !bleedmask_status){
    Hstr.str("");
    Hstr << "BLDINTRP= \'" << day << " " << month << " " << date
	 << " " << stime << " " << year << "\' / bleed trail interpolation";
    Inimage.AppendImageHeader(Hstr.str());
  }
  if(do_star_interp && !bleedmask_status){
    Hstr.str("");
    Hstr << "STRINTRP= \'" << day << " " << month << " " << date
	 << " " << stime << " " << year << "\' / saturated star interpolation";
    Inimage.AppendImageHeader(Hstr.str());
  }
  if(do_starmask && !bleedmask_status){
    Hstr.str("");
    Hstr << "STARMASK= \'" << day << " " << month << " " << date
	 << " " << stime << " " << year << "\' / created bright star mask";
    Inimage.AppendImageHeader(Hstr.str());
  }
  Hstr.str("");
  Hstr << "HISTORY DESDM: ";
  int argn = 0;
  while(argv[argn])
    Hstr << argv[argn++] << " ";
  Inimage.AppendImageHeader(Hstr.str());
  return(0);
}

int WriteOutputImage(FitsTools::FitsImage &Inimage,const std::string &ofilename,int flag_verbose)
{
  // profiler.FunctionEntry("ImageWrite");
  // Now actually write the output FITS file
  std::ostringstream Out;
  if(flag_verbose==3){
    Out.str("");
    Out << "Writing bleedtrail masked image into " << ofilename << ".";
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  if(Inimage.Write(ofilename,true,flag_verbose)){
    LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
    Out << "Failed write output image into  " << ofilename << ".";
    Inimage.Close();
    LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
    return(1);
  }
  // Close the output image.
  Inimage.Close();
  // profiler.FunctionExit("ImageWrite");
  if(flag_verbose==3){
    Out.str("");
    Out << "Finished writing " << ofilename;
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  return(0);
}
//  LocalWords:  LocY


int main(int argc,char *argv[])
{
  return(MakeBleedMask((const char **)argv));
}


