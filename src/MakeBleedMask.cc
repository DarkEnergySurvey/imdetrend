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
mkbleedmask [-aghimz] [-b <factor> -f <value> -l <value> -n <npix> -o <filename> -r <factor> -s <n> -t <npixels> 
            -v <level> -w <value> -x <filename> ] <infile> <outfile> 

  Input data
    <infile>   = FITS image (hopefully with astrometric soln)
    <outfile>  = FITS image output 
  Options:
  -g,--global_only:         Do not use local statistics.
  -s,--bgiters <n>:         Number of iterations for background estimation. (1)
  -b,--bgreject <factor>:   Use specified <factor> as scalefactor for background rejection. (5.0)
  -r,--growbox <factor>:    Factor by which to grow blob boxes for determining background level. (2)

  -t,--trail_length <npixels>: Use specified <npixels> for trail structuring element in Y. (20)
  -n,--numtrail <npix>:     Number of contiguous pixels required for trail detection. (trail_length/2)
  -f,--scalefactor <value>: Use specified <value> as scalefactor on background to detect bleedtrails. (10.0)

  -i,--interpolate:         Interpolate over bleedtrails.

  -l,--starlevel <value>:   Use specified <value> as scalefactor on background to detect stars. (5.0)

  -w,--growrad <value>:     Factor by which to grow detected star radii. (1.0)
  -a,--interpolate-stars:   Do radial interpolation over saturated stars for which bleedtrails were detected.
  -m,--starmask:            Create a mask for the detected bright objects. (No)
  -z,--zerostarweights:     Set weights to zero for masked stars.

  -o,--saturated_objects <filename>: Filename for saturated star table.(None)
  -x,--trailboxes <filename>: Filename for trail box table.(None)

  -h,--help:                Prints this long version of help.
  -v,--verb <level>:        Verbosity level [0-3], default = 1.

Summary:
  Detects bleedtrails in Y for incoming FITS image.
  Sets BADPIX_TRAIL for detected bleedtrails,
  and BADPIX_INTERP if interpolation is enabled,
  and BADPIX_STAR if star masking is enabled.
  Creates object table with WCS positions and radii.

Detailed Description:

Known Issues:


*/

#include "ComLine.hh"
#include "FitsTools.hh"
//#include "ImageMorphology.hh"
#include "BleedTrails.hh"
extern "C" {
#undef   OFF_T
#include "fitscat_mine.h"
#include "fitswcs.h"

static const char *svn_id = "$Id$";

}
// Local Function
void ModifyStarR(Morph::MaskDataType *mask,Morph::ImageDataType *image,
		 double cx,double cy,double &star_r,Morph::IndexType Nx,
		 Morph::IndexType Ny,Morph::MaskDataType rejection_mask,
		 double star_scalefactor,Morph::StatType &stats);

///
/// \brief ComLineObject for mkbleedmask app
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
    AddOption('g',"global_only");
    AddOption('h',"help");
    AddOption('i',"interpolate");
    AddOption('m',"starmask");
    AddOption('b',"bgreject",2,"factor");
    AddOption('f',"scalefactor",2,"value");
    AddOption('l',"starlevel",2,"value");
    AddOption('n',"numtrail",2,"npix");
    AddOption('o',"saturated_objects",2,"filename");
    AddOption('r',"growbox",2,"factor");
    AddOption('s',"bgiters",2,"n");
    AddOption('t',"trail_length",2,"npixels");
    AddOption('v',"verbose",2,"level");
    AddOption('e',"version");
    AddOption('w',"growrad",2,"value");
    AddOption('x',"trailboxes",2,"filename");
    AddOption('z',"zerostarweights");
    AddArgument("infile",1);
    AddArgument("outfile",1);
    AddHelp("interpolate-stars","Do radial interpolation over saturated stars for which bleedtrails were detected.");
    AddHelp("global_only","Do not use local statistics.");
    AddHelp("help","Prints this long version of help.");
    AddHelp("interpolate","Interpolate over bleedtrails.");
    AddHelp("starmask","Create a mask for the detected bright objects. (No)");
    AddHelp("bgreject","Use specified <factor> as scalefactor for background rejection. (5.0)");
    AddHelp("scalefactor","Use specified <value> as scalefactor on background to detect bleedtrails. (10.0)");
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



int MakeBleedMask(const char *argv[])
{
  // ------- Setup and Command Line Parsing --------
  
  // Initialization of command line object
  MakeBleedMaskComLine comline;
  int comline_status = comline.ProcessCommandLine((const char **)argv);

  std::string program_name(comline.ProgramName());
  
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
  bool do_star_interp    = !comline.GetOption("interpolate-stars").empty();
  std::string stw        =  comline.GetOption("trailwidth");
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
      Out << program_name << ": Invalid size for trail structuring element, " << trail_width << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
    }
  }

  double scalefactor = 5.0;
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

  // Step through the file lists and read in the FITS images

  // Declare and initialize FitsImage object which
  // will handle most of the FITS-specific interactions
  // under-the-covers.
  FitsTools::FitsImage Inimage;
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
  if(flag_verbose==3){
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

  Morph::IndexType Nx = Inimage.DES()->axes[0];
  Morph::IndexType Ny = Inimage.DES()->axes[1];
  Morph::IndexType npix = Nx*Ny;

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
  if(trail_width == 0)
    trail_width = 20;

  if(trail_width%2)
    trail_width++;
  
  int trail_arm = trail_width/2;
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

  unsigned int trail_pix = trail_width/2;
  std::vector<long> trail_structure(trail_width,0);
  std::vector<long>::iterator trailit = trail_structure.begin();
  for(unsigned int i = trail_pix;i > 0;i--)
    *trailit++ = -static_cast<long>(Nx*i);
  for(unsigned int i = 1;i <= trail_pix;i++)
    *trailit++ = i*Nx;
  
  short trail_rejection_mask  = BADPIX_CRAY | BADPIX_BPM;
  short trail_exception_mask  = 0;
  short trail_mask = BADPIX_TRAIL;
  short ground_rejection_mask = trail_rejection_mask | trail_mask | BADPIX_SATURATE; 
  int bleed_status = 0; 
  int total_trail_pixels = 0;
  int method = 1;
  

  // Get the "blobs" of BADPIX_SATURATE pixels. The number of blobs is the 
  // initial number of saturated objects before bleedtrail rejection.
  //
  // - blob_image is an integer image wherein the value indicates which
  //   blob to which the pixel belongs.
  // - blobs contains the pixel indices for each blob
  //
  std::vector<Morph::IndexType> blob_image(npix,0);
  std::vector<std::vector<Morph::IndexType> > blobs;
  // This call populates the above two data structures
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,BADPIX_SATURATE,blob_image,blobs);
  
  if(flag_verbose){
    Out.str("");
    Out << "Found " << blobs.size() << " saturated objects initially."; 
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  
    
  // THESE data structures hold the bounding boxes for
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
  
  // Get global image statistics
  Morph::BoxStats(Inimage.DES()->image,Inimage.DES()->mask,Nx,Ny,image_box,
		  BADPIX_SATURATE | BADPIX_CRAY |            // reject pixels with
		  BADPIX_BPM | BADPIX_STAR | BADPIX_TRAIL,   // these bits set
		  BADPIX_INTERP,                             // accept these
		  image_stats,image_npix_box);
  int bbb = 0;
  bool converged = false;
  double last_mean = image_stats[Image::IMMEAN];
  while((bbb < nbgit) && !converged){ 
    bbb++;
    Morph::ImageDataType ground_rejection_level = image_stats[Image::IMMEAN] + 
      ground_rejection_factor*image_stats[Image::IMSIGMA];
    // Get stats again, but reject high pixels
    Morph::BoxStats(Inimage.DES()->image,Inimage.DES()->mask,Nx,Ny,image_box,
		    BADPIX_SATURATE | BADPIX_CRAY |            // reject pixels with
		    BADPIX_BPM | BADPIX_STAR | BADPIX_TRAIL,   // these bits set
		    BADPIX_INTERP,                             // accept these
		    0,ground_rejection_level,                  // rejection levels
		    image_stats,image_npix_box);
    double residual = std::abs(image_stats[Image::IMMEAN]-last_mean);
    if(residual < 1e-3){
      converged = true;
      if(flag_verbose){
	Out.str("");
	Out << "Global image statistics converged after NITER=" << bbb 
	    << (bbb==1 ? " iteration." : " iterations.") << std::endl; 
	LX::ReportMessage(flag_verbose,QA,1,Out.str());
	Out.str("");
	Out << Util::stripdirs(filename) << "   BKGD=" << image_stats[Image::IMMEAN] << ", "
	    << "BKGD_SIGMA=" << image_stats[Image::IMSIGMA] << std::endl;
	LX::ReportMessage(flag_verbose,QA,1,Out.str());
      }
    }
    last_mean = image_stats[Image::IMMEAN];
  }
  if(!converged){
    Out.str("");
    Out << "Global image statistics did not converge after " << nbgit << " iterations.";
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }

  // Experimental - push saturated,interpolated pixels back up to detectable levels
  double trail_level = image_stats[Image::IMMEAN] + 2.0*scalefactor*image_stats[Image::IMSIGMA];
  for(int pixi = 0;pixi < Nx*Ny;pixi++){
    if((Inimage.DES()->mask[pixi]&BADPIX_SATURATE)&&(Inimage.DES()->mask[pixi]&BADPIX_INTERP)){
      Inimage.DES()->mask[pixi]^=BADPIX_INTERP;
      Inimage.DES()->image[pixi] = trail_level;
    }
  }

  std::vector<Morph::MaskDataType> temp_mask(npix,0);
  if(!Inimage.DES()->mask){
    method = 0;
    Inimage.DES()->mask = &temp_mask[0];
    if(flag_verbose){
      Out.str("");
      Out << program_name << "::Warning: Incoming image had no mask, using temporary.";
      LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
    }
  }
  else {
    std::vector<Morph::MaskDataType>::iterator  tmi = temp_mask.begin();
    while(tmi != temp_mask.end()){
      Morph::IndexType index = tmi - temp_mask.begin();
      *tmi++ = Inimage.DES()->mask[index];
    }
  }
    
  // Copy the input image into temp_image
  std::vector<Morph::ImageDataType> temp_image(npix,0);
  std::vector<Morph::ImageDataType>::iterator tii = temp_image.begin();
  while(tii != temp_image.end()){
    Morph::IndexType index = tii - temp_image.begin();
    *tii++ = Inimage.DES()->image[index];
  }
  
  // Loop through each blob, get image stats inside extended 
  // blob bounding box.
  //
  // Step 0 - Initiate loop over blobs
  std::vector<Morph::BlobType>::iterator bbbi   = blobs.begin();
  while(bbbi != blobs.end()){
    int blobno = bbbi - blobs.begin() + 1;
    Morph::BlobType &blob = *bbbi++;
    int nblobpix = blob.size();
    // Step 1 - Get blob bounding box
    std::vector<Morph::IndexType> box;
    std::sort(blob.begin(),blob.end());
    Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
    //  adds the box to our list of blob boxes
    blob_boxes.push_back(box);

    // Step 2 - Expand bounding boxes
    int bleed_over = (trail_width/2) - trail_arm + 1;
    assert(bleed_over > 0);
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
//     std::cout << "BackGround Box: [" << box[0] << ":" << box[1] << ","
// 	      << box[2] << ":" << box[3] << "]" << std::endl;
    
    // Step 3 - Get image statistics inside box
    Morph::IndexType npix_box;
    Morph::StatType stats;
    bool use_image_stats = false;
    //    if((nblobpix < trail_arm) || ){
    if(global_stats_only){
      use_image_stats = true;
    }
    else{
      Morph::BoxStats(Inimage.DES()->image,Inimage.DES()->mask,Nx,Ny,box,
		      BADPIX_SATURATE | BADPIX_CRAY |            // reject pixels with
		      BADPIX_BPM | BADPIX_STAR | BADPIX_TRAIL,   // these bits set
		      BADPIX_INTERP,                             // accept these
		      stats,npix_box);
      if(npix_box < 300)
	use_image_stats = true;
      else{
	//	std::cout << "Box had " << npix_box << " valid image pixels for statistics." << std::endl;
	//	std::cout << "Initial Stats: (" << stats[Image::IMMIN] << "," << stats[Image::IMMAX]
	//		  << "," << stats[Image::IMMEAN] << "," << stats[Image::IMSIGMA] << std::endl;
	bool stats_converged = false;
	bbb = 0;
	last_mean = stats[Image::IMMEAN];
	while((bbb < nbgit) && !stats_converged && !use_image_stats){
	  bbb++;
	  Morph::ImageDataType local_ground_rejection_level = stats[Image::IMMEAN] + 
	    ground_rejection_factor*stats[Image::IMSIGMA];
	  // Get stats again, but reject high pixels
	  Morph::BoxStats(Inimage.DES()->image,Inimage.DES()->mask,Nx,Ny,box,
			  BADPIX_SATURATE | BADPIX_CRAY |                  // reject pixels with
			  BADPIX_BPM | BADPIX_STAR | BADPIX_TRAIL,         // these bits set
			  BADPIX_INTERP,                                   // accept these
			  0,local_ground_rejection_level,                  // rejection levels
			  stats,npix_box);
	  if(npix_box < 300)
	    use_image_stats = true;
	  double residual = std::abs(stats[Image::IMMEAN] - last_mean);
	  if(residual < 1e-3){
	    stats_converged = true;
	    if(flag_verbose >= 5){
	      Out.str("");
	      Out << "Local image statistics converged after " << bbb << " iteration"
		  << (bbb==1 ? "." : "s."); 
	      LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	    }
	  }
	  last_mean = stats[Image::IMMEAN];
	  //	  std::cout << "Stats After Rejection(" << bbb << "): (" << stats[Image::IMMIN] << "," << stats[Image::IMMAX]
	  //		    << "," << stats[Image::IMMEAN] << "," << stats[Image::IMSIGMA] << std::endl;
	}
      }
    }
    if(use_image_stats){
      stats = image_stats;
      if(flag_verbose){
	Out.str("");
	Out << "Global image statistics used for blob " << blobno << ".";
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
      }
    }
    // Record statistics for this box in persistent data structure
    //    box_stats.push_back(stats);
    box_stats.push_back(stats);
  }

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
  

  // PASS 1 - DETECT EVERYTHING BRIGHT
  bleed_status = DetectBleedTrailsInBoxes(Inimage.DES()->image,
					  NULL,
					  Inimage.DES()->mask,
					  &temp_mask[0],
					  Nx,Ny,candidate_trail_boxes,box_stats,
					  trail_structure,
					  trail_rejection_mask,trail_exception_mask,
					  ground_rejection_mask,trail_exception_mask,
					  trail_mask,
					  BADPIX_INTERP,
					  BADPIX_STAR,
					  trail_arm,scalefactor,
					  star_scalefactor,
					  false,std::cout);
  // PASS 2 - DETECT BLEED TRAILS
  bleed_status = DetectBleedTrailsInBoxes(Inimage.DES()->image,
					  &temp_image[0],
					  &temp_mask[0],
					  Inimage.DES()->mask,
					  Nx,Ny,candidate_trail_boxes,box_stats,
					  trail_structure,
					  trail_rejection_mask,trail_exception_mask,
					  ground_rejection_mask,trail_exception_mask,
					  trail_mask,
					  BADPIX_INTERP,
					  BADPIX_STAR,
					  trail_arm,scalefactor,
					  star_scalefactor,
					  do_interp,std::cout);
  
  total_trail_pixels = bleed_status;

  if(bleed_status < 0)
    { 
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      Out.str("");
      Out << "DetectBleedTrails failed on pixel " << bleed_status << std::endl;
      LX::ReportMessage(flag_verbose,STATUS,5,Out.str());
      return(1);
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
  // with BADPIX_STAR as well.  We want to do "star interpolation" for
  // the pixels marked with BADPIX_STAR, instead of the run-of-the-mill
  // linear interpolation in the row for the normal bleed trail pixels.
  
  // Step 1 - Detect Star Blobs
  // Get the "blobs" of BADPIX_STAR pixels. The number of blobs is the 
  // new number of saturated objects after bleedtrail rejection.
  //
  // - blob_image is an integer image wherein the value indicates which
  //   blob to which the pixel belongs.
  // - blobs contains the pixel indices for each blob
  //
  blob_image.resize(npix,0);
  blobs.resize(0);
  // This call populates the above two data structures
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,BADPIX_STAR,blob_image,blobs);

  if(flag_verbose){
    Out.str("");
    Out << "Found " << blobs.size() << " saturated stars."; 
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }
  
    
  // Get boxes just around actual trails
  std::vector<Morph::BoxType>  trail_boxes;
  std::vector<std::vector<Morph::IndexType> > trail_blobs;
  blob_image.resize(npix,0);
  trail_blobs.resize(0);
  Morph::GetBlobs(Inimage.DES()->mask,Nx,Ny,trail_mask,blob_image,trail_blobs);
  bbbi   = trail_blobs.begin();
  while(bbbi != trail_blobs.end()){
    int blobno = bbbi - trail_blobs.begin() + 1;
    Morph::BlobType &blob = *bbbi++;
    int nblobpix = blob.size();
    // Step 1 - Get blob bounding box
    std::vector<Morph::IndexType> box;
    std::sort(blob.begin(),blob.end());
    Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
    //  adds the box to our list of blob boxes
    trail_boxes.push_back(box);
  }

  std::vector<Morph::BoxType>  star_boxes;
  std::vector<double> star_centers_x;
  std::vector<double> star_centers_y;
  std::vector<double> star_radii;
  std::vector<std::vector<double> > star_idf;
  // Loop through each blob, get blob bounding box,
  // star center, and pixel intensity distribution function
  //
  // Step 0 - Initiate loop over blobs
  std::ofstream BlobFile;
  std::vector<Morph::BlobType>::iterator bbi   = blobs.begin();
  while(bbi != blobs.end()){
    int blobno = bbi - blobs.begin() + 1;
    //      std::ostringstream BFOut;
    //      BFOut << "blob_" << blobno;
    //      BlobFile.open(BFOut.str().c_str());
    
    Morph::BlobType &blob = *bbi++;
    int nblobpix = blob.size();
    // Step 1 - Get blob bounding box
    Morph::BoxType box;
    std::sort(blob.begin(),blob.end());
    Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
    //  adds the box to our list of blob boxes
    star_boxes.push_back(box);
    
    
    // Step 2 - Get Center and number of pixels for intensity
    // distribution function
    Morph::IndexType nby = box[3] - box[2] + 1;
    Morph::IndexType nbx = box[1] - box[0] + 1;
    double cx = static_cast<double>(nbx);
    double cy = static_cast<double>(nby);
    cx /= 2.0;
    cy /= 2.0;
    double star_r = cy*cy;
    long r = static_cast<long>((std::sqrt(cx*cx + cy*cy)) + 2.0);
    std::vector<double> idf(r,0);
    //    std::vector<Morph::IndexType> npix(r,0);
    cx += (box[0]);
    cy += (box[2]);
    cx -= .5;
    cy -= .5;
    star_centers_x.push_back(cx);
    star_centers_y.push_back(cy);
    
    if(do_star_interp){
      //      BlobFile << "# Box: [" << box[0] << ":" << box[1] << ","
      //	       << box[2] << ":" << box[3] << "]" << std::endl;
      
      // Step 3 - Populate intensity distribution function
      // 
      // Go to the star's center
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
      //      BlobFile << "# Center (x1:x2,y1:y2) = (" << center_x1 << ":" << center_x2
      //	       << "," << center_y1 << ":" << center_y2 << ")" << std::endl;
      std::vector<double>::iterator idfit = idf.begin();
      long pixx = 0;
      long npix = 0;
      while(idfit != idf.end()){
	Morph::IndexType vindex = idfit++ - idf.begin();
	//	BlobFile << "# idf[" << vindex << "] = ";
	//	double &value = *idfit++;
	if((center_x1 + pixx) < Nx){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x1 + pixx;
	    //	    BlobFile << "+(" << center_x1+pixx << "," << center_y + yp << ","
	    //		     << cpindex << "," << Inimage.DES()->image[cpindex] << ")";
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if((center_x1 - pixx) >= 0){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x1 - pixx;
	    //	    BlobFile << "+(" << center_x1-pixx << "," << center_y + yp << ","
	    //		     << cpindex << "," << Inimage.DES()->image[cpindex] << ")";
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if((center_x2 + pixx) < Nx){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x2 + pixx;
	    //    BlobFile << "+(" << center_x2+pixx << "," << center_y + yp << ","
	    //		     << cpindex << "," << Inimage.DES()->image[cpindex] << ")";
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if((center_x2 - pixx) >= 0){
	  for(int yp = 0;yp < npix_y;yp++){
	    Morph::IndexType cpindex = (center_y + yp)*Nx + center_x2 - pixx;
	    //	    BlobFile << "+(" << center_x2-pixx << "," << center_y + yp << ","
	    //		     << cpindex << "," << Inimage.DES()->image[cpindex] << ")";
	    idf[vindex] += Inimage.DES()->image[cpindex];
	    npix++;
	  }
	}
	if(npix > 0){
	  idf[vindex] /= static_cast<double>(npix);
	  //	  BlobFile << "=" << idf[vindex] << std::endl;
	}
	npix = 0;
	pixx++;
      }
      idfit = idf.begin();
      while(idfit != idf.end()){
	//BlobFile << idfit - idf.begin() << " " << *idfit << std::endl;
	idfit++;
      }
      //      BlobFile.close();

      // Loop over blob pixels and set values according to linearly
      // interpolated values from the idf.
      Morph::BlobType::iterator sbi = blob.begin();
      while(sbi != blob.end()){
	Morph::IndexType pixel_index = *sbi++;
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
	if(!do_starmask)
	  Inimage.DES()->mask[pixel_index] ^= BADPIX_STAR;
      }
    }
    ModifyStarR(Inimage.DES()->mask,Inimage.DES()->image,cx,cy,star_r,Nx,Ny,
		ground_rejection_mask,star_scalefactor,image_stats);
    star_r *= (radius_growth_factor*radius_growth_factor);
    double srstarr = std::sqrt(star_r);
    star_radii.push_back(srstarr);
    if(do_starmask){
      //      star_radii[star_radii.size()-1] = std::sqrt(star_r);
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
	  if(r2 <= star_r)
	    Inimage.DES()->mask[bindex] |= BADPIX_STAR;
	  else if(Inimage.DES()->mask[bindex]&BADPIX_STAR)
	    Inimage.DES()->mask[bindex] ^= BADPIX_STAR;
	}
      }
    }  
  }
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
  }
  else if(!FitsTools::HeaderKeyExists(Inimage.ImageHeader(),"SCAMPFLG")){
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
    }
    else{
      tab=cat->tab;
      wcs_in=read_wcs(tab);
      // Get coordinate transformation (if it exists)
      if(FitsTools::HeaderKeyExists(Inimage.ImageHeader(),"CD1_1")){
	imjacobian[0] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD1_1");
	imjacobian[1] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD1_2");
	imjacobian[2] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD2_1");
	imjacobian[3] = FitsTools::GetHeaderValue<double>(Inimage.ImageHeader(),"CD3_1");
	convert_to_arcsecs = std::sqrt(imjacobian[0]*imjacobian[0] + imjacobian[1]*imjacobian[1]);
	convert_to_arcsecs += std::sqrt(imjacobian[2]*imjacobian[2] + imjacobian[3]*imjacobian[3]);
	convert_to_arcsecs *= 1800.0;
      }
      else{
	Out.str("");
	Out << "Failed to read coordinate transformation from " << Inimage.DES()->name 
	    << ". Writing object radii in pixels." << std::endl;
	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
	write_radii_in_pixels = true;
      }
    }
  }
  // Write the bright object table.
  if(!SaturatedObjectFileName.empty()){
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
      }
      else{
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

  // Make table for trail blob bounding boxes (if enabled)
  if(!TrailBoxesFileName.empty()){
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
    }
    else{
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
    std::vector<Morph::BoxType>::iterator tbi = trail_boxes.begin();
    while(tbi != trail_boxes.end()){
      unsigned int box_no = tbi-trail_boxes.begin() + 1;
      Morph::BoxType &box = *tbi++;
      if(get_wcs){
	double rawpos[2],wcspos[2];
	rawpos[0] = static_cast<double>(box[0])+1.0;
	rawpos[1] = static_cast<double>(box[2])+1.0;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX1.push_back(wcspos[0]);
	LocY1.push_back(wcspos[1]);

	rawpos[0] = static_cast<double>(box[1])+1.0;
	rawpos[1] = static_cast<double>(box[2])+1.0;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX2.push_back(wcspos[0]);
	LocY2.push_back(wcspos[1]);

	rawpos[0] = static_cast<double>(box[1])+1.0;
	rawpos[1] = static_cast<double>(box[3])+1.0;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX3.push_back(wcspos[0]);
	LocY3.push_back(wcspos[1]);

	rawpos[0] = static_cast<double>(box[0])+1.0;
	rawpos[1] = static_cast<double>(box[3])+1.0;
	raw_to_wcs(wcs_in,rawpos,wcspos);
	LocX4.push_back(wcspos[0]);
	LocY4.push_back(wcspos[1]);
      }
      else {
	LocX1.push_back(static_cast<double>(box[0])+1.0);
	LocY1.push_back(static_cast<double>(box[2])+1.0);
	LocX2.push_back(static_cast<double>(box[1])+1.0);
	LocY2.push_back(static_cast<double>(box[2])+1.0);
	LocX3.push_back(static_cast<double>(box[1])+1.0);
	LocY3.push_back(static_cast<double>(box[3])+1.0);
	LocX4.push_back(static_cast<double>(box[0])+1.0);
	LocY4.push_back(static_cast<double>(box[3])+1.0);
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

  // - - Output into FITS - - 
  
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
       << " " << stime << " " << year << "\' / bleed trail masking";
  Inimage.AppendImageHeader(Hstr.str());
  Hstr.str("");
  Hstr << "NBLEED  = " << bleed_status << "  / Number of bleed trail pixels.";
  Inimage.AppendImageHeader(Hstr.str());
  if(do_interp){
    Hstr.str("");
    Hstr << "BLDINTRP= \'" << day << " " << month << " " << date
	 << " " << stime << " " << year << "\' / bleed trail interpolation";
    Inimage.AppendImageHeader(Hstr.str());
  }
  if(do_star_interp){
    Hstr.str("");
    Hstr << "STRINTRP= \'" << day << " " << month << " " << date
	 << " " << stime << " " << year << "\' / saturated star interpolation";
    Inimage.AppendImageHeader(Hstr.str());
  }
  if(do_starmask){
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

  // Now actually write the output FITS file
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

  if(flag_verbose==3){
    Out.str("");
    Out << "Finished writing " << ofilename;
    LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
  }

  return(0);
}


int main(int argc,char *argv[])
{
  return(MakeBleedMask((const char **)argv));
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
  double temp_r = 5.0*std::sqrt(star_r);
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
      }
      else{
	median = starbins[bindex][nval/2];
      }
      if(median < starval){
	star_r = (bindex+1)*binsize;
	done = true;
      } 
    }
  }
  star_r *= star_r;
}

//  LocalWords:  LocY