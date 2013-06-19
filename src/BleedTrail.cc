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
      if(PixelShouldBeProcessed(blob_pixel_index)){
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
	while(PixelIsFlagged(row_pixel_index_minus,DataRejectionMask) && column >= 0){
	  row_pixel_index_minus--;
	  column--;
	}
	if(column >= 0){
	  value_minus = PixelValue(row_pixel_index_minus);
	}
	column = pixcolumn+1;
	Morph::IndexType row_pixel_index_plus = blob_pixel_index + 1;
	while(PixelIsFlagged(row_pixel_index_plus,DataRejectionMask) && column < nx){
	  row_pixel_index_plus++;
	  column++;
	}
	if(column < nx){
	  value_plus = PixelValue(row_pixel_index_plus);
	}
	if(value_minus >= star_level && value_plus >= star_level){
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
	while(PixelIsFlagged(row_pixel_index_minus,TrailMask) && column >= 0){
	  row_pixel_index_minus--;
	  column--;
	}
	if(column >= 0){
	  found_minus = true;
	  value_minus = PixelValue(row_pixel_index_minus);
	}
	column = pixcolumn+1;
	Morph::IndexType row_pixel_index_plus = blob_pixel_index + 1;
	while(PixelIsFlagged(row_pixel_index_plus,TrailMask) && column < nx){
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
	if(!found_minus && !found_plus){
	  std::cout << "Warning: Couldn't find either plus or minus" << std::endl;
	} else if(!found_minus || !found_plus) {
	  std::cout << "Warning: Couldn't find trail bounds." << std::endl;
	}
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
