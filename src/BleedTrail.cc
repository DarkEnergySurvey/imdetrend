/// 
/// @file
/// @brief Morphological bleed trail detection routines
/// @author Mike Campbell (mtcampbe@illinois.edu)
///
#include "BleedTrails.hh"

inline Morph::IndexType GetTrailWidth(Morph::ImageDataType *image,Morph::MaskDataType *mask,
				      Morph::IndexType Nx,Morph::IndexType Ny,Morph::IndexType index,
				      Morph::MaskDataType RejectMask,Morph::MaskDataType AcceptMask)
{
  Morph::IndexType y = index/Nx;
  Morph::IndexType x = index%Nx;
  Morph::IndexType sx = x;
  Morph::IndexType sy = y;
  Morph::IndexType sdist = 0;
  bool search_done = false;
  while(!search_done){
    sx--;
    sdist++;
    if(sx < 0)
      search_done = true;
    else{
      Morph::IndexType sindex = y*Nx+sx;
      if((!(mask[sindex]&RejectMask)) || (mask[sindex]&AcceptMask))
	search_done = true;
    }
  }
  search_done = false;
  sx = x;
  while(!search_done){
    sx++;
    sdist++;
    if(sx >= Nx)
      search_done = true;
    else{
      Morph::IndexType sindex = y*Nx+sx;
      if((!(mask[sindex]&RejectMask)) || (mask[sindex]&AcceptMask))
	search_done = true;
    }
  }
  return(sdist);
}



///
/// @brief Basic structuring element-based bleed trail detection
/// @param input_image Pointer to image data
/// @param output_image Pointer to output image data (only used if interpolation enabled)
/// @param input_mask  Pointer to input mask data
/// @param output_mask  Pointer to output mask data
/// @param Nx    Input number of columns in image
/// @param Ny    Input number of rows in image
/// @param blobs Input rectangular regions in which to detect bleedtrails
/// @param stats Input image statistics for each blob
/// @param trail_struct  Input structuring element for trail direction
/// @param TrailRejectionMask Input bit mask for pixels that should not be used in trail S.E.
/// @param TrailExceptionMask Input bit mask for rejection override in trail S.E.
/// @param ImageRejectionMask Input bit mask for pixels that should not be used in ground S.E.
/// @param ImageExceptionMask Input bit mask for rejection override in ground S.E.
/// @param TrailMask     Input bit mask for setting a bleedtrail pixel
/// @param InterpMask    Input bit mask for identifying interpolated pixel
/// @param StarMask      Input bit mask for identifying star pixel
/// @param trail_pix_limit Indicates how many bright pixels in a row determines trail pixels
/// @param scalefactor  Indicates how many sigma away from mean indicates bright
/// @param interpolate   Input flag for whether to interpolate over bleedtrails
/// @param Ostr Stream object for output
///
/// This function does the following:\n
/// * Loop over every pixel in the input blobs
/// * Evaluate the trail structuring element
/// * If trail_pix_limit bright pixels in a row are found 
///   in the structuring element, then mark current pixel
///   as a bleed trail with TrailMask
/// * If it's a bleed trail, check the valid image values
///   around the pixel to see if they are of star level, if
///   so mark the pixel as a star with StarMask.
/// * If interpolation is enabled, use the local mean (if it exists)
///   plus some randomized portion of sigma as the new image value
///
int DetectBleedTrailsInBlobs(Morph::ImageDataType *input_image,
			     Morph::ImageDataType *output_image,
			     Morph::MaskDataType *input_mask,
			     Morph::MaskDataType *output_mask,
			     Morph::IndexType Nx,
			     Morph::IndexType Ny,
			     std::vector<std::vector<Morph::IndexType> >  &blobs,
			     std::vector<Morph::StatType> &blob_stats,
			     std::vector<Morph::IndexType> &trail_struct,
			     Morph::MaskDataType TrailRejectionMask,
			     Morph::MaskDataType TrailExceptionMask,
			     Morph::MaskDataType ImageRejectionMask,
			     Morph::MaskDataType ImageExceptionMask,
			     Morph::MaskDataType TrailMask,
			     Morph::MaskDataType InterpMask,
			     Morph::MaskDataType StarMask,
			     Morph::IndexType trail_pix_limit,
			     double scalefactor,
			     double star_scalefactor,
			     bool interpolate,
			     std::ostream &Ostr)
{
  Morph::IndexType npix            = Nx*Ny;
  Morph::IndexType npix_trail      = 0;
  Morph::IndexType npix_bright     = 0;
  Morph::IndexType ntrail          = 0;
  Morph::IndexType tborder_y_minus = 0;
  Morph::IndexType tborder_y_plus  = 0;
  Morph::IndexType tborder_x_minus = 0;
  Morph::IndexType tborder_x_plus  = 0;
  Morph::ImageDataType *imptr = input_image;
  Morph::ImageDataType *omptr = output_image;
  
  long ranseed = -17;
  if(interpolate){
    if(output_image){
      for(Morph::IndexType i = 0;i < npix;i++)
	*omptr++ = *imptr++;
    } else {
      Ostr << "DetectTrailsInBlobs::Error: No output image supplied for interpolation.";
      return(-1);
    }
  }
			     
  //  std::vector<Morph::MaskDataType> temp_mask(&input_mask[0],&input_mask[npix-1]);
  Morph::GetSEAttributes(trail_struct,Nx,npix_trail,tborder_y_minus,tborder_y_plus,
			 tborder_x_minus,tborder_x_plus);
  
  Morph::IndexType yborder_minus = tborder_y_minus;
  Morph::IndexType yborder_plus  = tborder_y_plus;
  Morph::IndexType xborder_minus = tborder_x_minus;
  Morph::IndexType xborder_plus  = tborder_x_plus;
  Morph::IndexType ylimit        = Ny - yborder_plus;
  Morph::IndexType xlimit        = Nx - xborder_plus;
  
  std::vector<std::vector<Morph::IndexType> >::iterator bbi = blobs.begin();
  std::vector<Morph::StatType>::iterator bsi = blob_stats.begin();
  while(bbi != blobs.end()){
    Morph::IndexType blob_index = bbi - blobs.begin();
    std::vector<Morph::IndexType> &blob = *bbi++;
    Morph::StatType &stats = *bsi++;
    double trail_level = stats[Image::IMMEAN] + scalefactor*stats[Image::IMSIGMA];
    double star_level  = stats[Image::IMMEAN] + star_scalefactor*stats[Image::IMSIGMA];
    std::vector<Morph::IndexType>::iterator bi = blob.begin();
    while(bi != blob.end()){
      Morph::IndexType pixel_index = *bi++;
      // Make sure pixel at index isn't in the rejection mask before doing next set of things
      if(!(input_mask[pixel_index]&TrailRejectionMask) || 
	 (input_mask[pixel_index]&TrailExceptionMask)){
	Morph::IndexType pixx = pixel_index%Nx;
	Morph::IndexType pixy = pixel_index/Nx;
	bool need_trim = ((pixx < xborder_minus) || (pixx > xlimit) ||
			  (pixy < yborder_minus) || (pixy > ylimit));
	Morph::ImageDataType trail_mean = 0;
	Morph::ImageDataType trail_min = 0;
	Morph::ImageDataType trail_max = 0;
	Morph::ImageDataType trail_sigma = 0;
	Morph::IndexType trail_npix = 0;
	std::vector<Morph::IndexType> trimmed_trail_struct;
	std::vector<Morph::IndexType> *se_ptr = &trail_struct;
	if(need_trim){
	  Morph::TrimStructuringElement(pixel_index,Nx,Ny,trail_struct,trimmed_trail_struct);
	  se_ptr = &trimmed_trail_struct;
	}
	std::vector<Morph::IndexType> &tse(*se_ptr);
	std::vector<Morph::IndexType>::iterator selIt = tse.begin();
	npix_bright = 0;
	while(selIt != tse.end() && (npix_bright < trail_pix_limit)){
	  Morph::IndexType sindex = pixel_index + *selIt++;
	  if(!(input_mask[sindex]&TrailRejectionMask) || (input_mask[sindex]&TrailExceptionMask)){
	    if(input_image[sindex] > trail_level){
	      npix_bright++;
	    } else { 
	      npix_bright = 0;
	    }
	  }
	}
	if(npix_bright >= trail_pix_limit){
	  output_mask[pixel_index] |= TrailMask;
	  ntrail++;
	  std::vector<Morph::ImageDataType> vals(4,0);
	  std::vector<Morph::IndexType>     dist(4,0);
	  Morph::IndexType trailwidth = GetTrailWidth(input_image,input_mask,Nx,Ny,pixel_index,
						      ImageRejectionMask|TrailMask,ImageExceptionMask);
	  Morph::BoxType bgbox(4,0);
	  bgbox[0] = pixx - (trailwidth+2);
	  bgbox[1] = pixx + (trailwidth+2);
	  bgbox[2] = pixy - 1;
	  bgbox[3] = pixy + 1;
	  if(bgbox[0] < 0) bgbox[0] = 0;
	  if(bgbox[1] >= Nx) bgbox[1] = Nx-1;
	  if(bgbox[2] < 0) bgbox[2] = 0;
	  if(bgbox[3] >= Ny) bgbox[3] = Ny-1;
	  Morph::StatType interpstats;
	  Morph::IndexType npix_interp = 0;
	  Morph::BoxStats(input_image,input_mask,Nx,Ny,bgbox,
			  ImageRejectionMask|TrailMask,ImageExceptionMask,
			  interpstats,npix_interp);
	  Morph::ImageDataType baseval = interpstats[Image::IMMEAN];
	  if(!npix_interp){
	    bgbox[0] = pixx - 2*(trailwidth+2);
	    bgbox[1] = pixx + 2*(trailwidth+2);
	    bgbox[2] = pixy - 2;
	    bgbox[3] = pixy + 2;
	    if(bgbox[0] < 0) bgbox[0] = 0;
	    if(bgbox[1] >= Nx) bgbox[1] = Nx-1;
	    if(bgbox[2] < 0) bgbox[2] = 0;
	    if(bgbox[3] >= Ny) bgbox[3] = Ny-1;
	    Morph::BoxStats(input_image,input_mask,Nx,Ny,bgbox,
			    ImageRejectionMask|TrailMask,ImageExceptionMask,
			    interpstats,npix_interp);
	    if(!npix_interp)
	      baseval = stats[Image::IMMEAN];
	    else
	      baseval = interpstats[Image::IMMEAN];
	  }
	  if(baseval > star_level)
	    output_mask[pixel_index] |= StarMask;
	  if(interpolate){
	    output_image[pixel_index] = baseval + mygasdev(&ranseed)*stats[Image::IMSIGMA];
	    output_mask[pixel_index] |= InterpMask;
	  }
	}
      }
    }
  }
  return(ntrail);
}

int BleedTrailDetectionFilter::DetectBleedsOnBlobs(Morph::BlobsType &blobs,std::vector<Morph::ImageDataType> &high_levels)
{
  Morph::IndexType ntrailpix = 0;
  Initialize();
  unsigned int trail_detection_threshold = GetTrailDetectionThreshold();
  // Loop through all blobs
  Morph::BlobsType::iterator blobs_visitor = blobs.begin();
  std::vector<Morph::ImageDataType>::iterator level_iterator = high_levels.begin();
  while(blobs_visitor != blobs.end()){
    Morph::BlobType &blob(*blobs_visitor++);
    Morph::ImageDataType trail_level = *level_iterator++;
    SetTrailDetectionLevel(trail_level);
    // Loop through pixels in blob
    Morph::BlobType::iterator blob_pixel_visitor = blob.begin();
    while(blob_pixel_visitor != blob.end()){
      Morph::IndexType blob_pixel_index = *blob_pixel_visitor++;
      //      std::cout << "blob_pixel_index = " << blob_pixel_index << std::endl;
      if(ReferenceIsFlagged(blob_pixel_index)){
	  MarkAsTrailPixel(blob_pixel_index);
	  ntrailpix++;
      } else if(PixelShouldBeProcessed(blob_pixel_index)){
	// Get the structuring element for the filter this way because it automatically
	// detects intersections with the input image boundary and trims the structuring
	// element (if necessary).
	Morph::StructuringElementType &structuringelement(GetStructuringElement(blob_pixel_index));
	// Loop through structuring element pixels 
	Morph::StructuringElementType::iterator structure_offset_visitor = structuringelement.begin();
	unsigned int number_of_adjacent_bright_pixels = 0;
	while(structure_offset_visitor != structuringelement.end() && 
	      (number_of_adjacent_bright_pixels < trail_detection_threshold)){
	  Morph::IndexType considered_pixel_index = blob_pixel_index + 
	    ResolveAbsoluteOffset(*structure_offset_visitor++);
	  if(PixelDataIsValid(considered_pixel_index)){
	    if(PixelValueIsTrailLevel(considered_pixel_index)) {
	      number_of_adjacent_bright_pixels++;
	    } else {
	      number_of_adjacent_bright_pixels = 0;
	    }
	  }
	}
	if(number_of_adjacent_bright_pixels >= trail_detection_threshold){
	  MarkAsTrailPixel(blob_pixel_index);
	  ntrailpix++;
	}
      }
    }
  }
  return(ntrailpix);
};

int BleedTrailDetectionFilter::DetectStarsOnBlobs(Morph::BlobsType &blobs,std::vector<Morph::ImageDataType> &high_levels)
{
  Morph::IndexType nstarpixel = 0;
  Initialize();
  std::vector<bool> PixelHasNotBeenProcessed(npix,true);
  bool blob_on_left  = false;
  bool blob_on_right = false;
  // Loop through all blobs
  Morph::BlobsType::iterator blobs_visitor = blobs.begin();
  std::vector<Morph::ImageDataType>::iterator level_iterator = high_levels.begin();
  while(blobs_visitor != blobs.end()){
    Morph::BlobType &blob(*blobs_visitor++);
    Morph::ImageDataType star_level = *level_iterator++;
    SetStarDetectionLevel(star_level);
    // Loop through pixels in blob
    Morph::BlobType::iterator blob_pixel_visitor = blob.begin();
    while(blob_pixel_visitor != blob.end()){
      Morph::IndexType blob_pixel_index = *blob_pixel_visitor++;
      if(PixelShouldBeProcessed(blob_pixel_index) &&
	 PixelIsFlagged(blob_pixel_index,TrailMask) &&
	 PixelHasNotBeenProcessed[blob_pixel_index]){
	PixelHasNotBeenProcessed[blob_pixel_index] = false;
	Morph::ImageDataType value_minus = 0;
	Morph::ImageDataType value_plus = 0;
	Morph::IndexType pixcolumn = blob_pixel_index%nx;
	Morph::IndexType column = pixcolumn-1;
	Morph::IndexType row_pixel_index_minus = blob_pixel_index - 1;
	blob_on_left  = true;
	blob_on_right = true;
	while(PixelIsFlagged(row_pixel_index_minus,DataRejectionMask) && column >= 0){
	  row_pixel_index_minus--;
	  column--;
	}
	if(column >= 0){
	  blob_on_left = false;
	  value_minus = PixelValue(row_pixel_index_minus);
	}
	column = pixcolumn+1;
	Morph::IndexType row_pixel_index_plus = blob_pixel_index + 1;
	while(PixelIsFlagged(row_pixel_index_plus,DataRejectionMask) && column < nx){
	  row_pixel_index_plus++;
	  column++;
	}
	if(column < nx){
	  blob_on_right = false;
	  value_plus = PixelValue(row_pixel_index_plus);
	}
	if((blob_on_left && (value_plus > star_level)) || (blob_on_right && (value_minus > star_level)) || 
	   ((value_minus > star_level) && (value_plus > star_level))){
	  for(Morph::IndexType i = row_pixel_index_minus+1;i < row_pixel_index_plus;i++){
	    MarkAsStarPixel(i);
	    PixelHasNotBeenProcessed[i] = false;
	    nstarpixel++;
	  }
	}
      }
    }
  }
  return(nstarpixel);
};

int BleedTrailDetectionFilter::LinearInterpolateOverTrailsInBlobs(Morph::BlobsType &blobs)
{
  Initialize();
  int ninterppixel = 0;
  std::vector<bool> PixelHasNotBeenProcessed(npix,true);
  // Loop through all blobs
  Morph::BlobsType::iterator blobs_visitor = blobs.begin();
  while(blobs_visitor != blobs.end()){
    Morph::BlobType &blob(*blobs_visitor++);
    // Loop through pixels in blob
    Morph::BlobType::iterator blob_pixel_visitor = blob.begin();
    while(blob_pixel_visitor != blob.end()){
      Morph::IndexType blob_pixel_index = *blob_pixel_visitor++;
      if(PixelShouldBeProcessed(blob_pixel_index) &&
	 PixelIsFlagged(blob_pixel_index,TrailMask) &&
	 PixelHasNotBeenProcessed[blob_pixel_index]){
	PixelHasNotBeenProcessed[blob_pixel_index] = false;
	Morph::ImageDataType value_minus = 0;
	Morph::ImageDataType value_plus = 0;
	Morph::IndexType pixcolumn = blob_pixel_index%nx;
	Morph::IndexType column = pixcolumn-1;
	Morph::IndexType row_pixel_index_minus = blob_pixel_index - 1;
	bool found_minus = false;
	bool found_plus = false;
	while(PixelIsFlagged(row_pixel_index_minus,TrailMask|DataRejectionMask) && column >= 0){
	  row_pixel_index_minus--;
	  column--;
	}
	if(column >= 0){
	  found_minus = true;
	  value_minus = PixelValue(row_pixel_index_minus);
	}
	column = pixcolumn+1;
	Morph::IndexType row_pixel_index_plus = blob_pixel_index + 1;
	while(PixelIsFlagged(row_pixel_index_plus,TrailMask|DataRejectionMask) && column < nx){
	  row_pixel_index_plus++;
	  column++;
	}
	if(column < nx){
	  found_plus = true;
	  value_plus = PixelValue(row_pixel_index_plus);
	}
	if(value_minus == 0)
	  value_minus = value_plus;
	if(value_plus == 0)
	  value_plus = value_minus;
// 	if(!found_minus && !found_plus){
// 	  std::cout << "Warning: Couldn't find either plus or minus" << std::endl;
// 	} else if(!found_minus || !found_plus) {
// 	  std::cout << "Warning: Couldn't find trail bounds." << std::endl;
// 	}
	Morph::ImageDataType dval = value_plus - value_minus;
	Morph::ImageDataType dx   = row_pixel_index_plus - row_pixel_index_minus;
	Morph::ImageDataType slope = dval/dx;
	for(Morph::IndexType i = 1;i < dx;i++){
	  Morph::ImageDataType data = value_minus + slope*i;
	  Morph::IndexType iindex = row_pixel_index_minus + i;
	  if(PixelIsFlagged(iindex,StarMask))
	    data += StarNoise();
	  else
	    data += BackgroundNoise();
	  SetOutputPixelValue(iindex,data);
	  MarkAsInterpolatedPixel(iindex);
	  PixelHasNotBeenProcessed[i] = false;
	  ninterppixel++;
	}
      }
    }
  }
  return(ninterppixel);

};

// int EdgeBleedDetectionFilter::Run()
// {
//   // ******** Edge Bleed Detection ************
//   //
//   // This vector will store the bounding boxes
//   // for the edgebleeds (if any).
//   std::vector<Morph::BoxType> edgebleed_boxes;
  
//   // -- These configuration values are used for  --
//   // -- edgebleed detection and processing
//   // These y values will indicate the extent of
//   // of the edgebleed damage (if any).
//   std::vector<Morph::BoxType> edgebleed_detection_boxes; 
//   Morph::IndexType y00 = 0;
//   Morph::IndexType y01 = Ny;
//   Morph::IndexType y10 = 0;
//   Morph::IndexType y11 = Ny;
//   Morph::IndexType ampx = Nx/2;
//   Morph::IndexType nrows_detection = 20;
//   double initial_edgebleed_detection_level = 5.0;
//   int NPIX_EDGEBLEED    = 500;    // how large should the "low" section be to get caught
//   int NPIXMAX_EDGEBLEED = 200000; // larger than this? probably not an edgebleed
//   int XSIZE_EDGEBLEED   = 100;    // must be extended at least this far in X to be caught
//   // Using arbitrary number (200) here to be the max size
//   // of an edgebleed bounding box's extent in Y.  
//   Morph::IndexType ytop = Ny - (npix_edge+nrows_detection);
//   Morph::IndexType ybottom = (npix_edge+nrows_detection);
//   Morph::IndexType edgebleed_rows_left = 0;
//   Morph::IndexType edgebleed_rows_right = 0;
//   // For indicating whether edgebleed was 
//   // actually detected on a given edge.
//   bool bl_bleed                  = false;
//   bool tl_bleed                  = false;
//   bool br_bleed                  = false;
//   bool tr_bleed                  = false;
//   bool suspect_edge_bleed        = false;
//   bool suspect_edge_bleed_bottom = false;
//   bool suspect_edge_bleed_top    = false;
//   bool ultra_edgebleed_left      = false;
//   bool ultra_edgebleed_right     = false;
//   bool definite_edgebleed_left   = false;
//   bool definite_edgebleed_right  = false;
//   bool collision_on_left         = false;
//   bool collision_on_right        = false;
//   bool huge_edgebleed            = false;
//   bool small_edgebleed           = false;
//   bool small_left                = false;
//   bool small_right               = false;
//   // --------------------------

//   if(image_stats[Image::IMSIGMA] >= 200.0){
//     Out.str("");
//     Out << "Background sigma (" << image_stats[Image::IMSIGMA] << "). Detecting edgebleed at 1.0.";
//     LX::ReportMessage(flag_verbose,STATUS,1,Out.str());  
//     Out.str("");
//     initial_edgebleed_detection_level = 1.0;
//   }

//   // define the levels at which the threshold the image
//   // for detecting low pixels
//   double hole_level = image_stats[Image::IMMEAN] - initial_edgebleed_detection_level*image_stats[Image::IMSIGMA];
//   double high_level = std::numeric_limits<Morph::ImageDataType>::max();
//   if(hole_level < 0)
//     hole_level = 0;

//   // Get initial blobs of saturated pixels
//   Morph::BlobsType SatBlobs;
//   std::vector<Morph::IndexType> SatBlobImage(npix,0);
//   Morph::GetBlobsWithRejection(&temp_mask[0],Nx,Ny,BADPIX_SATURATE,BLOB_REJECTION_MASK,0,SatBlobImage,SatBlobs);
//   SatBlobImage.resize(0);

//   // This block just checks to see if we should expect there might be an edgebleed, which means:
//   // (0) A saturated blob collides with the read registers -or-
//   // (1) There is a blob of saturated pixels that runs the full 
//   //     extent of the image in the read direction -or-
//   // (2) There is a HUGE (> 150) pixel wide saturated blob that
//   //     extends more than 1500 pixels in the read direction 
//   Morph::BlobsType::iterator satblobit = SatBlobs.begin();
//   std::vector<unsigned int> full_length_blob_indices;
//   while(satblobit != SatBlobs.end()){
//     unsigned int bindex = satblobit - SatBlobs.begin();
//     Morph::BlobType &satblob(*satblobit++);
//     Morph::BoxType box;
//     std::sort(satblob.begin(),satblob.end());
//     Morph::GetBlobBoundingBox(satblob,Nx,Ny,box);
//     Morph::IndexType box_xsize = box[1] - box[0];
//     Morph::IndexType box_ysize = box[3] - box[2];
//     double xsizeb2 = static_cast<double>(box_xsize)/2.0;
//     Morph::IndexType approach = 20;
//     double xpos = static_cast<double>(box[1] + box[0])/2.0;
//     Morph::IndexType boxpos = static_cast<unsigned int>(xpos);
//     bool by_side = ((xpos < (xsizeb2+npix_edge)) || (std::abs(xpos-ampx) < xsizeb2) || 
// 		    ((Nx-xpos) < (xsizeb2+npix_edge)));
//     //    if(box_ysize > 1500 || by_side) {
//     if(box_xsize > 150 && box_ysize > 1500){
//       huge_edgebleed = true;
//       Out.str("");
//       approach = 200;
//       Out << "Suspect HUGE edgebleed case with box (" << box[0] << ":" 
// 	  << box[1] << "," << box[2] << ":" << box[3] << ")." << std::endl;
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str()); 
//     }
//     bool collision_with_bottom = false;
//     bool collision_with_top = false;
//     if((box[2] <= approach)){
//       if(read_register_on_bottom)
// 	suspect_edge_bleed_bottom = true;
//       if(box[2] <= npix_edge+2){
// 	if(boxpos < ampx)
// 	  collision_on_left = true;
// 	else
// 	  collision_on_right = true;
// 	collision_with_bottom = true;
//       }
//     }
//     if((box[3] >= (Ny - approach))){
//       if(read_register_on_top)
// 	suspect_edge_bleed_top = true;
//       if(box[3] >= (Ny-npix_edge-3)){
// 	if(boxpos < ampx)
// 	  collision_on_left = true;
// 	else
// 	  collision_on_right = true;
// 	collision_with_top = true;
//       }
//     }
//     if(collision_with_bottom && collision_with_top){
//       Out.str("");
//       Out << "Full length saturated blob(" << bindex << ") at (" << 
// 	box[0] << ":" << box[1] << "," << box[2] << ":" << box[3] << ")"
// 	  << " will be checked for edgebleed.";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
//       full_length_blob_indices.push_back(bindex);
//     }
//     //  }
//   }
//   suspect_edge_bleed = (suspect_edge_bleed_bottom || suspect_edge_bleed_top); 
//   if(suspect_edge_bleed){
//     // Give a quick status message to indicate potential edgebleed.
//     if(flag_verbose){
//       Out.str("");
//       Out << "Checking for possible edgebleed on "
// 	  << ((collision_on_left && collision_on_right) ? "both sides " : 
// 	      (collision_on_left ? "left side " : "right side "))
// 	  << "because of saturated blob proximity to read registers.";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
//     }
//     // Look through the rows by the read registers and find out
//     // if there are some rows with large numbers of negative or 
//     // anomalously low pixel values.
//     std::vector<Morph::IndexType> ndip(4,0);
//     Morph::IndexType nrows_check = 20;
//     if(huge_edgebleed) // For HUGE edgebleeds, check more rows
//       nrows_check = 1024;
//     Morph::IndexType bottom_y       = npix_edge;
//     Morph::IndexType top_y          = Ny - npix_edge - 1;
//     Morph::IndexType bottom_row     = (read_register_on_bottom ? 
// 				       (read_register_on_top ? top_y : bottom_y) : 
// 				       (top_y-nrows_check));
//     Morph::IndexType top_row        = (read_register_on_top ? 
// 				       (read_register_on_bottom ? bottom_y : top_y) : 
// 				       (bottom_y+nrows_check));
//     Morph::IndexType edgebleed_size = 0;

//     for(Morph::IndexType rowin = bottom_row;rowin <= top_row;rowin++){
//       Morph::IndexType index_base = rowin*Nx;
//       for(Morph::IndexType pixin = 0;pixin < ampx;pixin++){
// 	if(!(Inimage.DES()->mask[index_base+pixin]&(BADPIX_BPM|BLOB_REJECTION_MASK))){
// 	  Morph::ImageDataType left_value = Inimage.DES()->image[index_base+pixin];
// 	  if(left_value < hole_level){
// 	    ndip[0]++;
// 	    if(left_value < 0)
// 	      ndip[2]++;
// 	  }
// 	}
// 	if((!Inimage.DES()->mask[index_base+pixin+ampx]&(BADPIX_BPM|BLOB_REJECTION_MASK))){
// 	  Morph::ImageDataType right_value = Inimage.DES()->image[index_base+pixin+ampx];
// 	  if(right_value < hole_level){
// 	    ndip[1]++;
// 	    if(right_value < 0)
// 	      ndip[3]++;
// 	  }
// 	}
//       }
//     }
//     bool detected_dip_left     = ((ndip[0] > 30*nrows_check)&&collision_on_left);
//     bool detected_dip_right    = ((ndip[1] > 30*nrows_check)&&collision_on_right);
//     bool strong_dip_left       = ((ndip[0] > 60*nrows_check)&&collision_on_left);
//     bool strong_dip_right      = ((ndip[1] > 60*nrows_check)&&collision_on_right);
//     ultra_edgebleed_left       = ((ndip[0] > 500*nrows_check)&&collision_on_left);
//     ultra_edgebleed_right      = ((ndip[1] > 500*nrows_check)&&collision_on_right);
//     definite_edgebleed_left    = ((ndip[2] > 5)&&collision_on_left);
//     definite_edgebleed_right   = ((ndip[3] > 5)&&collision_on_right);
//     if((ndip[2] > 0 || ndip[3] > 0) && !(definite_edgebleed_left || definite_edgebleed_right)){
//       Out.str("");
//       Out << "Detected negative pixels near read registers: (" << ndip[2] << "," << ndip[3] 
// 	  << "), but edgebleed not definite.";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str()); 
//     }
//     if(definite_edgebleed_left || definite_edgebleed_right){
//       Out.str("");
//       Out << "Possible edgebleed with (" << ndip[2] << "," << ndip[3] << ") negative pixels near read registers."
// 	  << std::endl;
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str()); 
//       Morph::IndexType rrowrec = 0;
//       Morph::IndexType lrowrec = 0;
//       if(read_register_on_bottom){
// 	Morph::IndexType rowin = npix_edge;
// 	Morph::IndexType rowmax = Ny - npix_edge;
// 	while(rowin < rowmax){
// 	  Morph::IndexType index_base = rowin*Nx;
// 	  Morph::IndexType leftvals = 0;
// 	  Morph::IndexType rightvals = 0;
// 	  for(Morph::IndexType colin = 0; colin < ampx;colin++){
// 	    if(!(Inimage.DES()->mask[index_base+colin]&(BADPIX_BPM|BLOB_REJECTION_MASK)))
// 	      if(Inimage.DES()->image[index_base+colin] < 0) leftvals++;
// 	    if(!(Inimage.DES()->mask[index_base+colin+ampx]&(BADPIX_BPM|BLOB_REJECTION_MASK)))
// 	      if(Inimage.DES()->image[index_base+colin+ampx] < 0) rightvals++;
// 	  }	  
// 	  if(leftvals){
// 	    edgebleed_rows_left++;
// 	    lrowrec = rowin;
// 	  }
// 	  if(rightvals > 0){
// 	    edgebleed_rows_right++;
// 	    rrowrec = rowin;
// 	  }
// 	  rowin++;
// 	}
// 	if(edgebleed_rows_left > 0)
// 	  edgebleed_rows_left  = lrowrec;
// 	if(edgebleed_rows_right > 0)
// 	  edgebleed_rows_right = rrowrec;
//       }
//       if(read_register_on_top){
// 	Morph::IndexType rowin = Ny - npix_edge -1;
// 	while((rowin >=  npix_edge)){
// 	  Morph::IndexType index_base = rowin*Nx;
// 	  Morph::IndexType leftvals = 0;
// 	  Morph::IndexType rightvals = 0;
// 	  for(Morph::IndexType colin = 0; colin < ampx;colin++){
// 	    if(!(Inimage.DES()->mask[index_base+colin]&(BADPIX_BPM|BLOB_REJECTION_MASK)))
// 	      if(Inimage.DES()->image[index_base+colin] < 0) leftvals++;
// 	    if(!(Inimage.DES()->mask[index_base+colin+ampx]&(BADPIX_BPM|BLOB_REJECTION_MASK)))
// 	      if(Inimage.DES()->image[index_base+colin+ampx] < 0) rightvals++;
// 	  }	  
// 	  if(leftvals){
// 	    edgebleed_rows_left++;
// 	    lrowrec = rowin;
// 	  }
// 	  if(rightvals){
// 	    rrowrec = rowin;
// 	    edgebleed_rows_right++;
// 	  }
// 	  rowin--;
// 	}   
// 	if(edgebleed_rows_left > 0)
// 	  edgebleed_rows_left  = 4096 - lrowrec + 1;
// 	if(edgebleed_rows_right > 0)
// 	  edgebleed_rows_right = 4096 - rrowrec + 1;
//       }
//       Out.str("");
//       Out << "Possible edgebleed - number of rows with negative pixels: (" << edgebleed_rows_left << "," 
// 	  << edgebleed_rows_right << ")";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str());        
//     }
//     if(detected_dip_left || detected_dip_right){
//       Out.str("");
//       Out << "Number of Low Pixels near read registers: (" << ndip[0] << "," 
// 	  << ndip[1] << ")";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str());  
//     } 
//     if(!(strong_dip_left || strong_dip_right)){
//       Out.str("");
//       LX::ReportMessage(flag_verbose,STATUS,1,"No strong edgebleed dips detected near read registers.");  
//       if(!full_length_blob_indices.empty()){
// 	Out.str("");
// 	std::vector<unsigned int>::iterator flti = full_length_blob_indices.begin();
// 	while(flti != full_length_blob_indices.end()){
// 	  unsigned int bindex = *flti++;
// 	  Morph::BlobType &satblob(SatBlobs[bindex]);
// 	  Morph::BlobType::iterator blob_pixel_index = satblob.begin();
// 	  Morph::BoxType box;
// 	  Morph::GetBlobBoundingBox(satblob,Nx,Ny,box);
// 	  Out.str("");
// 	  Out << "Full length saturated blob(" << bindex << ") at (" << 
// 	    box[0] << ":" << box[1] << "," << box[2] << ":" << box[3] << ")"
// 	      << " being checked for small, otherwise undetected edgebleed.";
// 	  LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
// 	  Out.str("");
// 	  unsigned int xpos = (box[1] + box[0])/2;
// 	  unsigned int npix_detect = 0;
// 	  Morph::IndexType npix_blob = satblob.size();
// 	  while(blob_pixel_index != satblob.end()){
// 	    unsigned int pixind = *blob_pixel_index++;
// 	    unsigned int pixcol = pixind%Nx;
// 	    unsigned int pixrow = pixind/Nx;
// 	    pixind = pixrow*Nx + (2*ampx - pixcol - 1);
// 	    if(!(Inimage.DES()->mask[pixind]&(BADPIX_BPM|BLOB_REJECTION_MASK))){
// 	      Morph::ImageDataType pixval = Inimage.DES()->image[pixind];
// 	      if(pixval < hole_level)
// 		npix_detect++;
// 	    }
// 	  }
// 	  double fraction = (1.0*std::fabs(npix_blob - npix_detect))/(1.0*npix_blob);
// 	  Out << "Fraction of mirror pixels below background for potential edgebleed: " << fraction << ".";
// 	  LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
// 	  if(fraction > .5){
// 	    small_edgebleed = true;
// 	    if(xpos < ampx){
// 	      small_left = true;
// 	      LX::ReportMessage(flag_verbose,STATUS,1,"Suspect small edgebleed on left.");
// 	    }
// 	    else {
// 	      small_right = true;
// 	      LX::ReportMessage(flag_verbose,STATUS,1,"Suspect small edgebleed on right.");
// 	    }
// 	  }
// 	}
//       }
//     }
//     if(!small_edgebleed && !((definite_edgebleed_right || definite_edgebleed_left)&&huge_edgebleed)){
//       LX::ReportMessage(flag_verbose,STATUS,1,"No evidence for edgebleed detected. Pushing further detection level to 5.0.");
//       hole_detection_level *= 5.0;
//     }
//     if((definite_edgebleed_right || definite_edgebleed_left) && huge_edgebleed){
//       hole_detection_level = 1.0;
//       Out.str("");
//       Out << "Possible anomalously large edgebleed detected, resetting hole detection to "
// 	  << hole_detection_level << ".";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str());  
//     }
//     hole_level = image_stats[Image::IMMEAN] - hole_detection_level*image_stats[Image::IMSIGMA];
//     if(hole_level < 0){
//       hole_level = 0;
//       Out.str("");
//       Out << "Hole detection level for edgebleed should not be < 0, resetting to 0.";
//       LX::ReportMessage(flag_verbose,STATUS,1,Out.str());
//       Out.str("");
//     }
//     // This sets BADPIX_LOW (in a temporary mask) for the low pixels
//     std::vector<Morph::MaskDataType> hole_mask(npix,0);
//     std::vector<Morph::IndexType> hole_image(npix,0);
//     std::vector<std::vector<Morph::IndexType> > hole_blobs;
//     Morph::ThresholdImage(Inimage.DES()->image,&hole_mask[0],Nx,Ny,
// 			  hole_level,high_level,BADPIX_LOW,BADPIX_SATURATE);  

//     // Erosion followed by Dilation will removed straggly connections 
//     // but preserve the effective extent of the bounding box for the large
//     // part of the "blob" - hopefully help to clean up these vastly overflagged
//     // edgebleeds.
//     Morph::ErodeMask(&hole_mask[0],Nx,Ny,structuring_element,BADPIX_LOW);
//     Morph::DilateMask(&hole_mask[0],Nx,Ny,structuring_element,BADPIX_LOW);

//     // Update the temporary mask so that it knows about the image's BPM and EDGE
//     Morph::MaskDataType bits_to_set = BADPIX_BPM|BADPIX_EDGE;
//     for(int jj = 0;jj < npix;jj++)
//       hole_mask[jj] |= (Inimage.DES()->mask[jj]&bits_to_set);
      
//     // Get the blobs of BADPIX_LOW while respecing the BPM and/or EDGE
//     Morph::GetBlobsWithRejection(&hole_mask[0],Nx,Ny,BADPIX_LOW,BADPIX_SATURATE|BLOB_REJECTION_MASK,
// 				 0,hole_image,hole_blobs);
  
//     // Loop through the blobs and if a large (i.e. larger than NPIX_EDGEBLEED)
//     // blob is found, then check whether it collides with the image boundary.
//     // If so, it is a positive early edgebleed detection
//     std::vector<Morph::BlobType>::iterator hbi   = hole_blobs.begin();
//     Morph::IndexType edgebleed_y_lower = 0;
//     Morph::IndexType edgebleed_y_upper = 0;
//     while(hbi != hole_blobs.end()){
//       int blobno = hbi - hole_blobs.begin() + 1;
//       Morph::BlobType &blob = *hbi++;
//       int nblobpix = blob.size();
//       std::vector<Morph::IndexType> box;
//       // Sort the blob before bounding box determination
//       std::sort(blob.begin(),blob.end());
//       Morph::GetBlobBoundingBox(blob,Nx,Ny,box);
//       int xsize = box[1] - box[0];
//       if((((nblobpix > NPIX_EDGEBLEED) && (xsize > XSIZE_EDGEBLEED)) || (detected_dip_left || detected_dip_right))
// 	 && (nblobpix < NPIXMAX_EDGEBLEED)){ // blob is of proper size
// 	if(box[0] < ampx){ // "left" side of chip
// 	  // box collides with bottom image boundary (and not the top)
// 	  if((box[2] <= ybottom) && (box[3] < ytop) && read_register_on_bottom){
// 	    bl_bleed = true;
// 	    if(box[3] > y00)
// 	      y00 = box[3];
// 	    // box collides with top image boundary (and not the bottom)
// 	  } else if((box[3] >= ytop) && (box[2] > ybottom) && read_register_on_top){ 
// 	    tl_bleed = true;
// 	    if(box[2] < y01)
// 	      y01 = box[2];
// 	  } 
// 	} 
// 	if(box[1] >= ampx) { // "right" side of chip
// 	  // box collides with bottom image boundary (and not the top)
// 	  if((box[2]  <= ybottom) && (box[3] < ytop) && read_register_on_bottom){
// 	    br_bleed = true;
// 	    if(box[3] > y10)
// 	      y10 = box[3];
// 	    // box collides with top image boundary (and not the bottom)
// 	  } else if((box[3] >= ytop) && (box[2] > ybottom) && read_register_on_top){
// 	    tr_bleed = true;
// 	    if(box[2] < y11)
// 	      y11 = box[2];
// 	  } 
// 	}
//       }
//     }
//     // ******** End Edge Bleed Detection ************
  
//     // If possible edge bleed is detected, then recalculate image statistics
//     // with edge bleed masked out - and then make sure to use global
//     // image statistics and *not* local statistics for the rest of the
//     // processing
//     if(bl_bleed || tl_bleed || br_bleed || tr_bleed || definite_edgebleed_left || definite_edgebleed_right){
//       // Dilate twice again if edge bleed.
//       //    Morph::DilateMask(&temp_mask[0],Nx,Ny,structuring_element,BADPIX_SATURATE);
//       //    Morph::DilateMask(&temp_mask[0],Nx,Ny,structuring_element,BADPIX_SATURATE);
//       global_stats_only = true;
//       std::vector<Morph::MaskDataType> temp_mask2(temp_mask.begin(),temp_mask.end());
//       if(bl_bleed){
// 	small_left = false;
// 	Out.str("");
// 	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (0:" 
// 	    << ampx-1 << ",0:" << y00 << ")";
// 	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
// 	std::vector<Morph::IndexType> box(4,0);
// 	box[0] = 0;
// 	box[1] = ampx-1;
// 	box[2] = 0;
// 	box[3] = y00;
// 	edgebleed_boxes.push_back(box);
// 	for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	  Morph::IndexType iminy = by*Nx;
// 	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	    Morph::IndexType imin = iminy + bx;
// 	    temp_mask2[imin] |= BADPIX_LOW;
// 	  }
// 	} 
//       }
//       if(tl_bleed){
// 	small_left = false;
// 	Out.str("");
// 	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (0:" << ampx-1 
// 	    << "," << y01 << ":" << Ny-1 << ").";
// 	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
// 	std::vector<Morph::IndexType> box(4,0);
// 	box[0] = 0;
// 	box[1] = ampx-1;
// 	box[2] = y01;
// 	box[3] = Ny-1;
// 	edgebleed_boxes.push_back(box);      
// 	for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	  Morph::IndexType iminy = by*Nx;
// 	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	    Morph::IndexType imin = iminy + bx;
// 	    temp_mask2[imin] |= BADPIX_LOW;
// 	  }
// 	} 
//       }
//       if(br_bleed){
// 	small_right = false;
// 	Out.str("");
// 	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (" << ampx << ":" 
// 	    << Nx-1 << ",0:" << y10 << ").";
// 	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());    
// 	std::vector<Morph::IndexType> box(4,0);
// 	box[0] = ampx;
// 	box[1] = Nx-1;
// 	box[2] = 0;
// 	box[3] = y10;
// 	edgebleed_boxes.push_back(box);
// 	for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	  Morph::IndexType iminy = by*Nx;
// 	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	    Morph::IndexType imin = iminy + bx;
// 	    temp_mask2[imin] |= BADPIX_LOW;
// 	  }
// 	} 
//       }
//       if(tr_bleed){
// 	small_right = false;
// 	Out.str("");
// 	Out << "Warning: Detected possible edgebleed with low blobs near read registers at (" << ampx << ":" 
// 	    << Nx-1 << "," << y11 << ":" << Ny - 1 << ")."; 
// 	LX::ReportMessage(flag_verbose,STATUS,3,Out.str());    
// 	std::vector<Morph::IndexType> box(4,0);
// 	box[0] = ampx;
// 	box[1] = Nx-1;
// 	box[2] = y11;
// 	box[3] = Ny-1;
// 	edgebleed_boxes.push_back(box);      
// 	for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	  Morph::IndexType iminy = by*Nx;
// 	  for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	    Morph::IndexType imin = iminy + bx;
// 	    temp_mask2[imin] |= BADPIX_LOW;
// 	  }
// 	} 
//       }
//       if(definite_edgebleed_left){
// 	small_left = false;
// 	if(read_register_on_bottom){
// 	  bl_bleed = true;
// 	  Out.str("");
// 	  Out << "Warning: Detected definite edgebleed with low blobs near read registers at (0:" 
// 	      << ampx-1 << ",0:" << edgebleed_rows_left << ")";
// 	  LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
// 	  std::vector<Morph::IndexType> box(4,0);
// 	  box[0] = 0;
// 	  box[1] = ampx-1;
// 	  box[2] = 0;
// 	  box[3] = edgebleed_rows_left;
// 	  edgebleed_boxes.push_back(box);
// 	  for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	    Morph::IndexType iminy = by*Nx;
// 	    for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	      Morph::IndexType imin = iminy + bx;
// 	      temp_mask2[imin] |= BADPIX_LOW;
// 	    }
// 	  } 
// 	} 
// 	if(read_register_on_top){
// 	  tl_bleed = true;
// 	  Out.str("");
// 	  Out << "Warning: Detected definite edgebleed with low blobs near read registers at (0:" << ampx-1 
// 	      << "," << 4096-edgebleed_rows_left-1 << ":" << Ny-1 << ").";
// 	  LX::ReportMessage(flag_verbose,STATUS,3,Out.str());
// 	  std::vector<Morph::IndexType> box(4,0);
// 	  box[0] = 0;
// 	  box[1] = ampx-1;
// 	  box[2] = 4096 - edgebleed_rows_left - 1;
// 	  box[3] = Ny-1;
// 	  edgebleed_boxes.push_back(box);      
// 	  for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	    Morph::IndexType iminy = by*Nx;
// 	    for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	      Morph::IndexType imin = iminy + bx;
// 	      temp_mask2[imin] |= BADPIX_LOW;
// 	    }
// 	  } 
// 	}
//       }
//       if(definite_edgebleed_right){
// 	small_right = false;
// 	if(read_register_on_bottom){
// 	  br_bleed = true;
// 	  Out.str("");
// 	  Out << "Warning: Detected definite edgebleed with low blobs near read registers at (" << ampx << ":" 
// 	      << Nx-1 << ",0:" << edgebleed_rows_right << ").";
// 	  LX::ReportMessage(flag_verbose,STATUS,3,Out.str());    
// 	  std::vector<Morph::IndexType> box(4,0);
// 	  box[0] = ampx;
// 	  box[1] = Nx-1;
// 	  box[2] = 0;
// 	  box[3] = edgebleed_rows_right;
// 	  edgebleed_boxes.push_back(box);
// 	  for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	    Morph::IndexType iminy = by*Nx;
// 	    for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	      Morph::IndexType imin = iminy + bx;
// 	      temp_mask2[imin] |= BADPIX_LOW;
// 	    }
// 	  } 
// 	}
// 	if(read_register_on_top){
// 	  tr_bleed = true;
// 	  Out.str("");
// 	  Out << "Warning: Detected definite edgebleed with low blobs near read registers at (" << ampx << ":" 
// 	      << Nx-1 << "," << 4096-edgebleed_rows_right-1 << ":" << Ny - 1 << ")."; 
// 	  LX::ReportMessage(flag_verbose,STATUS,3,Out.str());    
// 	  std::vector<Morph::IndexType> box(4,0);
// 	  box[0] = ampx;
// 	  box[1] = Nx-1;
// 	  box[2] = 4096-edgebleed_rows_right-1;
// 	  box[3] = Ny-1;
// 	  edgebleed_boxes.push_back(box);      
// 	  for(Morph::IndexType by = box[2];by <= box[3];by++){
// 	    Morph::IndexType iminy = by*Nx;
// 	    for(Morph::IndexType bx = box[0];bx <= box[1];bx++){
// 	      Morph::IndexType imin = iminy + bx;
// 	      temp_mask2[imin] |= BADPIX_LOW;
// 	    }
// 	  } 
// 	}
//       }
//       Out.str("");
//       if(Morph::GetSky(Inimage.DES()->image,&temp_mask2[0],Nx,Ny,
// 		       minpix,nbgit,ground_rejection_factor,1e-3,
// 		       ground_rejection_mask,BADPIX_INTERP,
// 		       image_stats,image_npix_box,niter,&Out)){
// 	std::ostringstream GSO;
// 	GSO << "Global image statistics failed:" << std::endl
// 	    << Out.str() << std::endl;
// 	LX::ReportMessage(flag_verbose,STATUS,4,GSO.str());
// 	bleedmask_status = 1;
// 	UpdateImageHeader(Inimage,argv,0,do_interp,do_star_interp,
// 			  do_starmask,bleedmask_status);
// 	return(WriteOutputImage(Inimage,ofilename,flag_verbose));
//       } else if(flag_verbose){
// 	Out << "Global image statistics (recalculated due to strongly suspected edgebleed)"
// 	    << " converged (again) after NITER=" << niter 
// 	    << (niter==1 ? " iteration." : " iterations.") << std::endl; 
// 	LX::ReportMessage(flag_verbose,QA,1,Out.str());
// 	Out.str("");
// 	Out << Util::stripdirs(filename) << "   BKGD=" << image_stats[Image::IMMEAN] << ", "
// 	    << "BKGD_SIGMA=" << image_stats[Image::IMSIGMA] << std::endl;
// 	LX::ReportMessage(flag_verbose,QA,1,Out.str());
//       }
//     } 
//     if (small_edgebleed){
//       LX::ReportMessage(flag_verbose,STATUS,1,
// 			"Strongly suspect small edgebleed due to intrachip xtalk, but otherwise undetected.");
//       if(small_left){
// 	if(read_register_on_bottom){
// 	  bl_bleed = true;
// 	  unsigned int test = 10+npix_edge;
// 	  if(y00 < test)
// 	    y00 = test;
// 	  LX::ReportMessage(flag_verbose,STATUS,1,"Flagging 10 edgebleed rows on bottom left to be safe.");
// 	} else {
// 	  tl_bleed = true;
// 	  unsigned int test = Ny - 1 - npix_edge - 10;
// 	  if(y01 > test)
// 	    y01 = test;
// 	  LX::ReportMessage(flag_verbose,STATUS,1,"Flagging 10 edgebleed rows on top left to be safe.");
// 	}
//       }
//       if(small_right){
// 	if(read_register_on_bottom){
// 	  br_bleed = true;
// 	  unsigned int test = 10+npix_edge;
// 	  if(y10 < test)
// 	    y10 = test;
// 	  LX::ReportMessage(flag_verbose,STATUS,1,"Flagging 10 edgebleed rows on bottom right to be safe.");
// 	} else {
// 	  tr_bleed = true;
// 	  LX::ReportMessage(flag_verbose,STATUS,1,"Flagging 10 edgebleed rows on top right to be safe.");
// 	  unsigned int test = Ny - 1 - npix_edge - 10;
// 	  if(y11 > test) 
// 	    y11 = test;
// 	}
//       }
//     }
//   } // suspect_edge_bleed
// }
