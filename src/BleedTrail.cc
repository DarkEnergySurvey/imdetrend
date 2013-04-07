/// 
/// \file
/// \brief Morphological bleed trail detection routines
/// \author Mike Campbell (mtcampbe@illinois.edu)
///
#include "ImageMorphology.hh"
//#include "NRFuncs.hh"
//#include <iostream>

///
/// \brief Basic structuring element-based bleed trail detection
/// \param image Pointer to image data
/// \param output_image Pointer to output image data (only used if interpolation enabled)
/// \param input_mask  Pointer to input mask data
/// \param output_mask  Pointer to output mask data
/// \param Nx    Input number of columns in image
/// \param Ny    Input number of rows in image
/// \param trail_struct  Input structuring element for trail direction
/// \param ground_struct Input struct element for background direction
/// \param TrailRejectionMask Input bit mask for pixels that should not be used in trail S.E.
/// \param TrailExceptionMask Input bit mask for rejection override in trail S.E.
/// \param GroundRejectionMask Input bit mask for pixels that should not be used in ground S.E.
/// \param GroundExceptionMask Input bit mask for rejection override in ground S.E.
/// \param TrailMask     Input bit mask for setting a bleedtrail pixel
/// \param InterpMask    Input bit mask for identifying interpolated pixel
/// \param trail_pix_limit Indicates how many bright pixels in a row determines trail pixels
/// \param scalefactor  Indicates how many sigma away from mean indicates bright
/// \param interpolate   Input flag for whether to interpolate over bleedtrails
/// \param Ostr Stream object for output
///
/// This function does the following:\n
/// - Loop over every pixel in the input image
///   - Evaluate the structuring elements:
///     - trail_struct (shaped like the object you are trying to detect)
///     - ground_struct (shaped to sample local background)
///   - Compare the mean of the two structuring elements
///     - if (\f$\bar{p_t} > (\bar{p_g} + {\sigma}_g)\f$) then set BADPIX_TRAIL
///
/// \note The actual implementation of the procedure breaks the image into 
///       chunks in order to do special treatment of image regions where the
///       structuring element collides with the image boundaries.
///
int DetectBleedTrails(Morph::ImageDataType *image,
		      Morph::ImageDataType *output_image,
		      Morph::MaskDataType *input_mask,
		      Morph::MaskDataType *output_mask,
		      Morph::IndexType Nx,
		      Morph::IndexType Ny,
		      std::vector<Morph::IndexType> &trail_struct,
		      std::vector<Morph::IndexType> &ground_struct,
		      Morph::MaskDataType TrailRejectionMask,
		      Morph::MaskDataType TrailExceptionMask,
		      Morph::MaskDataType GroundRejectionMask,
		      Morph::MaskDataType GroundExceptionMask,
		      Morph::MaskDataType TrailMask,
		      Morph::MaskDataType InterpMask,
		      Morph::IndexType trail_pix_limit,
		      double scalefactor,
		      bool interpolate,
		      std::ostream &Ostr)
{
  Morph::IndexType npix = Nx*Ny;
  Morph::IndexType npix_processed  = 0;
  Morph::IndexType npix_trail      = 0;
  Morph::IndexType npix_bright     = 0;
  Morph::IndexType ntrail          = 0;
  Morph::IndexType npix_ground     = 0;
  Morph::IndexType tborder_y_minus = 0;
  Morph::IndexType tborder_y_plus  = 0;
  Morph::IndexType tborder_x_minus = 0;
  Morph::IndexType tborder_x_plus  = 0;
  Morph::IndexType gborder_y_minus = 0;
  Morph::IndexType gborder_y_plus  = 0;
  Morph::IndexType gborder_x_minus = 0;
  Morph::IndexType gborder_x_plus  = 0;
  Morph::ImageDataType *imptr = image;
  Morph::ImageDataType *omptr = output_image;
  Morph::ImageDataType ground_rejection_factor = 2.0;

  long ranseed = -17;
  if(interpolate)
    for(Morph::IndexType i = 0;i < npix;i++)
      *omptr++ = *imptr++;
  
    
  //  Morph::MaskDataType  GroundRejectionMask = RejectionMask|TrailMask|BADPIX_SATURATE;
  //  Morph::MaskDataType  GroundRejectionMask = RejectionMask|TrailMask;
  std::vector<Morph::MaskDataType> temp_mask(&input_mask[0],&input_mask[npix-1]);
  Morph::GetSEAttributes(trail_struct,Nx,npix_trail,tborder_y_minus,tborder_y_plus,
			 tborder_x_minus,tborder_x_plus);
  Morph::GetSEAttributes(ground_struct,Nx,npix_ground,gborder_y_minus,gborder_y_plus,
			 gborder_x_minus,gborder_x_plus);
  Morph::IndexType yborder_minus = (tborder_y_minus > gborder_y_minus ? tborder_y_minus :
			gborder_y_minus);
  Morph::IndexType yborder_plus  = (tborder_y_plus > gborder_y_plus ? tborder_y_plus :
			gborder_y_plus);
  Morph::IndexType xborder_minus = (tborder_x_minus > gborder_x_minus ? tborder_x_minus :
			gborder_x_minus);
  Morph::IndexType xborder_plus  = (tborder_x_plus > gborder_x_plus ? tborder_x_plus :
			gborder_x_plus);
  Morph::IndexType ylimit = Ny - yborder_plus;
  Morph::IndexType xlimit = Nx - xborder_plus;
  
  std::string messageHdr("DetectBleedTrails: ");
  Ostr << messageHdr << "Morphological Boundary (x-:x+,y-:y+) = ("
       << xborder_minus << "," << xborder_plus << "," << yborder_minus
       << "," << yborder_plus << ")" << std::endl;
  // Do the slow parts (i.e. parts near image borders)
  //
  // Slow Part 1 [1:Nx,1:border_y] (INSIDE THE LEFT Y BORDER)
  Ostr << messageHdr << "Processing Y- Border" << std::endl;
  for(Morph::IndexType y = 0;y < yborder_minus;y++){
    for(Morph::IndexType x = 0;x < Nx;x++){
      npix_processed++;
      Morph::IndexType index = y*Nx + x;
      // Make sure pixel at index isn't in the rejection mask before doing next set of things
      if(!(temp_mask[index]&TrailRejectionMask) || (temp_mask[index]&TrailExceptionMask)){
	Morph::ImageDataType trail_mean = 0;
	Morph::ImageDataType trail_min = 0;
	Morph::ImageDataType trail_max = 0;
	Morph::ImageDataType trail_sigma = 0;
	Morph::IndexType trail_npix = 0;
	Morph::ImageDataType ground_mean = 0;
	Morph::ImageDataType ground_min = 0;
	Morph::ImageDataType ground_max = 0;
	Morph::ImageDataType ground_sigma = 0;
	Morph::IndexType ground_npix = 0;
	std::vector<Morph::IndexType> trimmed_trail_struct;
	std::vector<Morph::IndexType> trimmed_ground_struct;
	Morph::TrimStructuringElement(index,Nx,Ny,ground_struct,trimmed_ground_struct);
	Morph::TrimStructuringElement(index,Nx,Ny,trail_struct,trimmed_trail_struct);
	//	Morph::TrimStructuringElement(index,Nx,Ny,ground_struct,trimmed_ground_struct);
	//	Morph::TrimStructuringElement(index,Nx,Ny,trail_struct,trimmed_trail_struct);
	//	Morph::SEStats(image,&temp_mask[0],index,trimmed_trail_struct,TrailRejectionMask,
	//		       TrailExceptionMask,trail_min,trail_max,trail_mean,trail_sigma,trail_npix);
	Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix);
	Morph::ImageDataType ground_rejection_level = ground_mean + ground_rejection_factor*ground_sigma;
	Morph::IndexType ground_npix_reject = 0;
	Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix,
		       ground_rejection_level,ground_npix_reject);	
	std::vector<Morph::IndexType>::iterator selIt = trimmed_trail_struct.begin();
	npix_bright = 0;
	while(selIt != trimmed_trail_struct.end() && (npix_bright < trail_pix_limit)){
	  Morph::IndexType sindex = index + *selIt++;
	  if(!(temp_mask[sindex]&TrailRejectionMask) || (temp_mask[sindex]&TrailExceptionMask)){
	    if(image[sindex] > (ground_mean + scalefactor*ground_sigma))
	      npix_bright++;
	    else if(npix_bright) npix_bright = 0;
	  }
	}
	if(npix_bright >= trail_pix_limit){
	 output_mask[index] |= TrailMask;
	 if(interpolate){
	   output_image[index] = ground_mean + mygasdev(&ranseed)*ground_sigma;
	   output_mask[index] |= InterpMask;
	 }
	 ntrail++;
	}
	//	else
	//	  output_mask[index] ^= (temp_mask[index]&TrailMask);
	// 	if(!((trail_npix > 0) && (ground_npix > 0))){
	// 	  Ostr << messageHdr << "trail element size = " 
	// 	       << trimmed_trail_struct.size() << std::endl
	// 	       << messageHdr << "ground element size  = " 
	// 	       << trimmed_ground_struct.size() << std::endl
	// 	       << messageHdr << "trail_npix = " << trail_npix 
	// 	       << ", ground_npix = " << ground_npix << std::endl
	// 	       << messageHdr << "index = " << index << std::endl
	// 	       << messageHdr << "temp_mask[index] = " << temp_mask[index] 
	// 	       << std::endl;
	// 	  return(-npix_processed);
	// 	}
	// Mark it as a bleed trail if it's outside where the background should be
	// 	if(trail_mean > (ground_mean + scalefactor*ground_sigma)){ 
	// 	 output_mask[index] |= TrailMask;
	// 	  ntrail++;
	// 	}
      }
    }
  }

  // Do the X-borders for the parts of the image outside of Y borders
  if(xborder_minus > 0 || xborder_plus > 0){
    Ostr << messageHdr << "Processing X-/X+ borders " << std::endl;
    // LEFT X-BORDER
    for(Morph::IndexType y = yborder_minus;y < ylimit;y++){
      for(Morph::IndexType x = 0;x < xborder_minus;x++){
	npix_processed++;
	Morph::IndexType index = y*Nx + x;
	// Make sure pixel at index isn't in the rejection mask before doing next set of things
	if(!(temp_mask[index]&TrailRejectionMask) || (temp_mask[index]&TrailExceptionMask)){
	  Morph::ImageDataType trail_mean = 0;
	  Morph::ImageDataType trail_min = 0;
	  Morph::ImageDataType trail_max = 0;
	  Morph::ImageDataType trail_sigma = 0;
	  Morph::IndexType trail_npix = 0;
	  Morph::ImageDataType ground_mean = 0;
	  Morph::ImageDataType ground_min = 0;
	  Morph::ImageDataType ground_max = 0;
	  Morph::ImageDataType ground_sigma = 0;
	  Morph::IndexType ground_npix = 0;
	  std::vector<Morph::IndexType> trimmed_trail_struct;
	  std::vector<Morph::IndexType> trimmed_ground_struct;
	  Morph::TrimStructuringElement(index,Nx,Ny,ground_struct,trimmed_ground_struct);
	  Morph::TrimStructuringElement(index,Nx,Ny,trail_struct,trimmed_trail_struct);
	//	  Morph::TrimStructuringElement(index,npix,ground_struct,trimmed_ground_struct);
	//	  Morph::TrimStructuringElement(index,npix,trail_struct,trimmed_trail_struct);
	  //	  assert(trimmed_ground_struct.size() > 0 && trimmed_trail_struct.size() > 0);
// 	  Morph::SEStats(image,&temp_mask[0],index,trimmed_trail_struct,TrailRejectionMask,
// 			 TrailExceptionMask,trail_min,trail_max,trail_mean,trail_sigma,trail_npix);
	  Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
			 GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix);
	Morph::ImageDataType ground_rejection_level = ground_mean + ground_rejection_factor*ground_sigma;
	Morph::IndexType ground_npix_reject = 0;
	Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix,
		       ground_rejection_level,ground_npix_reject);	
	  //	  assert((trail_npix > 0) && (ground_npix > 0));
	  std::vector<Morph::IndexType>::iterator selIt = trimmed_trail_struct.begin();
	  npix_bright = 0;
	  while(selIt != trimmed_trail_struct.end() && (npix_bright < trail_pix_limit)){
	    Morph::IndexType sindex = index + *selIt++;
	    if(!(temp_mask[sindex]&TrailRejectionMask) || (temp_mask[sindex]&TrailExceptionMask)){
	      //	    if(!(mask[sindex]&RejectionMask) || (mask[sindex]&ExceptionMask)){
	      if(image[sindex] > (ground_mean + scalefactor*ground_sigma))
		npix_bright++;
	      else if(npix_bright) npix_bright = 0;
	    }
	  }
	  if(npix_bright >= trail_pix_limit){
	   output_mask[index] |= TrailMask;
	   if(interpolate){
	     output_image[index] = ground_mean + mygasdev(&ranseed)*ground_sigma;
	     output_mask[index] |= InterpMask;
	   }
	   ntrail++;
	  }
	  //	else
	  //	  output_mask[index] ^= (temp_mask[index]&TrailMask);
// 	  if(!((trail_npix > 0) && (ground_npix > 0)))
// 	    return(-npix_processed);
// 	  // Mark it as a bleed trail if it's outside where the background should be
// 	  if(trail_mean > (ground_mean + scalefactor*ground_sigma)){
// 	   output_mask[index] |= TrailMask;
// 	    ntrail++;
// 	  }
	}
      }
      // RIGHT X-BORDER
      for(Morph::IndexType x = xlimit;x < Nx;x++){
	npix_processed++;
	Morph::IndexType index = y*Nx + x;
	// Make sure pixel at index isn't in the rejection mask before doing next set of things
	if(!(temp_mask[index]&TrailRejectionMask) || (temp_mask[index]&TrailExceptionMask)){
	  Morph::ImageDataType trail_mean = 0;
	  Morph::ImageDataType trail_min = 0;
	  Morph::ImageDataType trail_max = 0;
	  Morph::ImageDataType trail_sigma = 0;
	  Morph::IndexType trail_npix = 0;
	  Morph::ImageDataType ground_mean = 0;
	  Morph::ImageDataType ground_min = 0;
	  Morph::ImageDataType ground_max = 0;
	  Morph::ImageDataType ground_sigma = 0;
	  Morph::IndexType ground_npix = 0;
	  std::vector<Morph::IndexType> trimmed_trail_struct;
	  std::vector<Morph::IndexType> trimmed_ground_struct;
	  Morph::TrimStructuringElement(index,Nx,Ny,ground_struct,trimmed_ground_struct);
	  Morph::TrimStructuringElement(index,Nx,Ny,trail_struct,trimmed_trail_struct);
	  //	  Morph::TrimStructuringElement(index,npix,ground_struct,trimmed_ground_struct);
	  //	  Morph::TrimStructuringElement(index,npix,trail_struct,trimmed_trail_struct);
	  //	  assert(trimmed_ground_struct.size() > 0 && trimmed_trail_struct.size() > 0);
	  // 	  Morph::SEStats(image,&temp_mask[0],index,trimmed_trail_struct,TrailRejectionMask,
	  // 			 TrailExceptionMask,trail_min,trail_max,trail_mean,trail_sigma,trail_npix);
	  Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
			 GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix);
	  //	  assert((trail_npix > 0) && (ground_npix > 0));
	Morph::ImageDataType ground_rejection_level = ground_mean + ground_rejection_factor*ground_sigma;
	Morph::IndexType ground_npix_reject = 0;
	Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix,
		       ground_rejection_level,ground_npix_reject);	
	  std::vector<Morph::IndexType>::iterator selIt = trimmed_trail_struct.begin();
	  npix_bright = 0;
	  while(selIt != trimmed_trail_struct.end() && (npix_bright < trail_pix_limit)){
	    Morph::IndexType sindex = index + *selIt++;
	    if(!(temp_mask[sindex]&TrailRejectionMask) || (temp_mask[sindex]&TrailExceptionMask)){
	    //	    if(!(mask[sindex]&RejectionMask) || (mask[sindex]&ExceptionMask)){
	      if(image[sindex] > (ground_mean + scalefactor*ground_sigma))
		npix_bright++;
	      else if(npix_bright) npix_bright = 0;
	    }
	  }
	  if(npix_bright >= trail_pix_limit){
	   output_mask[index] |= TrailMask;
	   if(interpolate){
	     output_image[index] = ground_mean + mygasdev(&ranseed)*ground_sigma;
	     output_mask[index] |= InterpMask;
	   }
	   ntrail++;
	  }
	  //	else
	  //	  output_mask[index] ^= (temp_mask[index]&TrailMask);
// 	  if(!((trail_npix > 0) && (ground_npix > 0)))
// 	    return(-npix_processed);
// 	  // Mark it as a bleed trail if it's outside where the background should be
// 	  if(trail_mean > (ground_mean + scalefactor*ground_sigma)){
// 	   output_mask[index] |= TrailMask;
// 	    ntrail++;
// 	  }
	}
      }
    }
  }

  Ostr << messageHdr << "Processing Y+ Border " << std::endl;
  // Now do the RIGHT Y-border (if it exists) 
  for(Morph::IndexType y = ylimit;y < Ny;y++){
    for(Morph::IndexType x = 0;x < Nx;x++){
      npix_processed++;
      Morph::IndexType index = y*Nx + x;
      // Make sure pixel at index isn't in the rejection mask before doing next set of things
      if(!(temp_mask[index]&TrailRejectionMask) || (temp_mask[index]&TrailExceptionMask)){
	Morph::ImageDataType trail_mean = 0;
	Morph::ImageDataType trail_min = 0;
	Morph::ImageDataType trail_max = 0;
	Morph::ImageDataType trail_sigma = 0;
	Morph::IndexType trail_npix = 0;
	Morph::ImageDataType ground_mean = 0;
	Morph::ImageDataType ground_min = 0;
	Morph::ImageDataType ground_max = 0;
	Morph::ImageDataType ground_sigma = 0;
	Morph::IndexType ground_npix = 0;
	std::vector<Morph::IndexType> trimmed_trail_struct;
	std::vector<Morph::IndexType> trimmed_ground_struct;
	Morph::TrimStructuringElement(index,Nx,Ny,ground_struct,trimmed_ground_struct);
	Morph::TrimStructuringElement(index,Nx,Ny,trail_struct,trimmed_trail_struct);
	//	Morph::TrimStructuringElement(index,npix,ground_struct,trimmed_ground_struct);
	//	Morph::TrimStructuringElement(index,npix,trail_struct,trimmed_trail_struct);
	//	assert(trimmed_ground_struct.size() > 0 && trimmed_trail_struct.size() > 0);
	// 	Morph::SEStats(image,&temp_mask[0],index,trimmed_trail_struct,TrailRejectionMask,
	// 		       TrailExceptionMask,trail_min,trail_max,trail_mean,trail_sigma,trail_npix);
	Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix);
	//	assert((trail_npix > 0) && (ground_npix > 0));
	Morph::ImageDataType ground_rejection_level = ground_mean + ground_rejection_factor*ground_sigma;
	Morph::IndexType ground_npix_reject = 0;
	Morph::SEStats(image,&temp_mask[0],index,trimmed_ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix,
		       ground_rejection_level,ground_npix_reject);	
	std::vector<Morph::IndexType>::iterator selIt = trimmed_trail_struct.begin();
	npix_bright = 0;
	while(selIt != trimmed_trail_struct.end() && (npix_bright < trail_pix_limit)){
	  Morph::IndexType sindex = index + *selIt++;
	  if(!(temp_mask[sindex]&TrailRejectionMask) || (temp_mask[sindex]&TrailExceptionMask)){
	    //	  if(!(mask[sindex]&RejectionMask) || (mask[sindex]&ExceptionMask)){
	    if(image[sindex] > (ground_mean + scalefactor*ground_sigma))
	      npix_bright++;
	    else if(npix_bright) npix_bright = 0;
	  }
	}
	if(npix_bright >= trail_pix_limit){
	  if(interpolate){
	    output_image[index] = ground_mean + mygasdev(&ranseed)*ground_sigma;
	    output_mask[index] |= InterpMask;
	  }
	  output_mask[index] |= TrailMask;
	  ntrail++;
	}
	//	else
	//	  output_mask[index] ^= (temp_mask[index]&TrailMask);
// 	if(!((trail_npix > 0) && (ground_npix > 0)))
// 	  return(-npix_processed);
// 	// Mark it as a bleed trail if it's outside where the background should be
// 	if(trail_mean > (ground_mean + scalefactor*ground_sigma)){
// 	 output_mask[index] |= TrailMask;
// 	  ntrail++;  
// 	}
      }
    }
  }
  
  Ostr << messageHdr << "Processing central part of image." << std::endl;
  // Now do the central part of the image 
  for(Morph::IndexType y = yborder_minus;y < ylimit;y++){
    for(Morph::IndexType x = xborder_minus;x < xlimit;x++){
      npix_processed++;
      Morph::IndexType index = y*Nx + x;
      // Make sure pixel at index isn't in the rejection mask before doing next set of things
      if(!(temp_mask[index]&TrailRejectionMask) || (temp_mask[index]&TrailExceptionMask)){
	Morph::ImageDataType trail_mean = 0;
	Morph::ImageDataType trail_min = 0;
	Morph::ImageDataType trail_max = 0;
	Morph::ImageDataType trail_sigma = 0;
	Morph::IndexType trail_npix = 0;
	Morph::ImageDataType ground_mean = 0;
	Morph::ImageDataType ground_min = 0;
	Morph::ImageDataType ground_max = 0;
	Morph::ImageDataType ground_sigma = 0;
	Morph::IndexType ground_npix = 0;
// 	Morph::SEStats(image,&temp_mask[0],index,trail_struct,TrailRejectionMask,
// 		       TrailExceptionMask,trail_min,trail_max,trail_mean,trail_sigma,trail_npix);
	Morph::SEStats(image,&temp_mask[0],index,ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix);
	//	assert((trail_npix > 0) && (ground_npix > 0));
	Morph::ImageDataType ground_rejection_level = ground_mean + ground_rejection_factor*ground_sigma;
	Morph::IndexType ground_npix_reject = 0;
	Morph::SEStats(image,&temp_mask[0],index,ground_struct,GroundRejectionMask,
		       GroundExceptionMask,ground_min,ground_max,ground_mean,ground_sigma,ground_npix,
		       ground_rejection_level,ground_npix_reject);	
	std::vector<Morph::IndexType>::iterator selIt = trail_struct.begin();
	npix_bright = 0;
	while(selIt != trail_struct.end() && (npix_bright < trail_pix_limit)){
	  Morph::IndexType sindex = index + *selIt++;
	  if(!(temp_mask[sindex]&TrailRejectionMask) || (temp_mask[sindex]&TrailExceptionMask)){
	    //	  if(!(mask[sindex]&RejectionMask) || (mask[sindex]&ExceptionMask)){
	    if(image[sindex] > (ground_mean + scalefactor*ground_sigma))
	      npix_bright++;
	    else if(npix_bright) npix_bright = 0;
	  }
	}
	if(npix_bright >= trail_pix_limit){
	  output_mask[index] |= TrailMask;
	  if(interpolate){
	    output_mask[index] |= InterpMask;
	    output_image[index] = ground_mean + mygasdev(&ranseed)*ground_sigma;
	  }
	  ntrail++;
	}
	//	else
	//	  output_mask[index] ^= (temp_mask[index]&TrailMask);
// 	if(!((trail_npix > 0) && (ground_npix > 0)))
// 	    return(-npix_processed);
// 	// Mark it as a bleed trail if it's outside where the background should be
// 	  if(trail_mean > (ground_mean + scalefactor*ground_sigma)){
// 	   output_mask[index] |= TrailMask;
// 	    ntrail++;
// 	  }
      }
    }
  }
  //  assert(npix_processed == npix);
  if(npix_processed != npix)
    return(-npix_processed);
  return(ntrail);
}


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
/// \brief Basic structuring element-based bleed trail detection
/// \param input_image Pointer to image data
/// \param output_image Pointer to output image data (only used if interpolation enabled)
/// \param input_mask  Pointer to input mask data
/// \param output_mask  Pointer to output mask data
/// \param Nx    Input number of columns in image
/// \param Ny    Input number of rows in image
/// \param boxes Input rectangular regions in which to detect bleedtrails
/// \param box_stats Input image statistics for boxes
/// \param trail_struct  Input structuring element for trail direction
/// \param TrailRejectionMask Input bit mask for pixels that should not be used in trail S.E.
/// \param TrailExceptionMask Input bit mask for rejection override in trail S.E.
/// \param ImageRejectionMask Input bit mask for pixels that should not be used in ground S.E.
/// \param ImageExceptionMask Input bit mask for rejection override in ground S.E.
/// \param TrailMask     Input bit mask for setting a bleedtrail pixel
/// \param InterpMask    Input bit mask for identifying interpolated pixel
/// \param StarMask      Input bit mask for identifying star pixel
/// \param trail_pix_limit Indicates how many bright pixels in a row determines trail pixels
/// \param scalefactor  Indicates how many sigma away from mean indicates bright
/// \param interpolate   Input flag for whether to interpolate over bleedtrails
/// \param Ostr Stream object for output
///
/// This function does the following:\n
/// - Loop over every pixel in the input image
///   - Evaluate the structuring elements:
///     - trail_struct (shaped like the object you are trying to detect)
///     - ground_struct (shaped to sample local background)
///   - Compare the mean of the two structuring elements
///     - if (\f$\bar{p_t} > (\bar{p_g} + {\sigma}_g)\f$) then set BADPIX_TRAIL
///
/// \note The actual implementation of the procedure breaks the image into 
///       chunks in order to do special treatment of image regions where the
///       structuring element collides with the image boundaries.
///
int DetectBleedTrailsInBoxes(Morph::ImageDataType *input_image,
			     Morph::ImageDataType *output_image,
			     Morph::MaskDataType *input_mask,
			     Morph::MaskDataType *output_mask,
			     Morph::IndexType Nx,
			     Morph::IndexType Ny,
			     std::vector<Morph::BoxType>   &boxes,
			     std::vector<Morph::StatType> &box_stats,
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
  if(interpolate)
    for(Morph::IndexType i = 0;i < npix;i++)
      *omptr++ = *imptr++;
  
  //  std::vector<Morph::MaskDataType> temp_mask(&input_mask[0],&input_mask[npix-1]);
  Morph::GetSEAttributes(trail_struct,Nx,npix_trail,tborder_y_minus,tborder_y_plus,
			 tborder_x_minus,tborder_x_plus);

  Morph::IndexType yborder_minus = tborder_y_minus;
  Morph::IndexType yborder_plus  = tborder_y_plus;
  Morph::IndexType xborder_minus = tborder_x_minus;
  Morph::IndexType xborder_plus  = tborder_x_plus;
  Morph::IndexType ylimit        = Ny - yborder_plus;
  Morph::IndexType xlimit        = Nx - xborder_plus;
  

  std::vector<Morph::BoxType>::iterator bi = boxes.begin();
  std::vector<Morph::StatType>::iterator bsi = box_stats.begin();
  while(bi != boxes.end()){
    Morph::StatType &stats = *bsi++;
    Morph::BoxType &box = *bi++;
    bool need_trim = ((box[0] < xborder_minus) || (box[1] > xlimit) ||
		      (box[2] < yborder_minus) || (box[3] > ylimit));
    double trail_level = stats[Image::IMMEAN] + scalefactor*stats[Image::IMSIGMA];
    double star_level  = stats[Image::IMMEAN] + star_scalefactor*stats[Image::IMSIGMA];
    if(need_trim){
      for(Morph::IndexType y = box[2];y <= box[3];y++){
	for(Morph::IndexType x = box[0];x <= box[1];x++){
	  Morph::IndexType index = y*Nx + x;
	  // Make sure pixel at index isn't in the rejection mask before doing next set of things
	  if(!(input_mask[index]&TrailRejectionMask) || (input_mask[index]&TrailExceptionMask)){
	    Morph::ImageDataType trail_mean = 0;
	    Morph::ImageDataType trail_min = 0;
	    Morph::ImageDataType trail_max = 0;
	    Morph::ImageDataType trail_sigma = 0;
	    Morph::IndexType trail_npix = 0;
	    std::vector<Morph::IndexType> trimmed_trail_struct;
	    Morph::TrimStructuringElement(index,Nx,Ny,trail_struct,trimmed_trail_struct);
	    std::vector<Morph::IndexType>::iterator selIt = trimmed_trail_struct.begin();
	    npix_bright = 0;
	    while(selIt != trimmed_trail_struct.end() && (npix_bright < trail_pix_limit)){
	      Morph::IndexType sindex = index + *selIt++;
	      if(!(input_mask[sindex]&TrailRejectionMask) || (input_mask[sindex]&TrailExceptionMask)){
		if(input_image[sindex] > trail_level)
		  npix_bright++;
		else 
		  npix_bright = 0;
	      }
	    }
	    if(npix_bright >= trail_pix_limit){
	      output_mask[index] |= TrailMask;
	      if(interpolate){
		std::vector<Morph::ImageDataType> vals(4,0);
		std::vector<Morph::IndexType>     dist(4,0);
		Morph::IndexType trailwidth = GetTrailWidth(input_image,input_mask,Nx,Ny,index,
							    ImageRejectionMask|TrailMask,ImageExceptionMask);
		//		assert(trailwidth < 20);
		Morph::BoxType bgbox(4,0);
		bgbox[0] = x - (trailwidth+4);
		bgbox[1] = x + (trailwidth+4);
		bgbox[2] = y - 1;
		bgbox[3] = y + 1;
		if(bgbox[0] < 0) bgbox[0] = 0;
		if(bgbox[1] >= Nx) bgbox[1] = Nx-1;
		if(bgbox[2] < 0) bgbox[2] = 0;
		if(bgbox[3] >= Ny) bgbox[3] = Ny-1;
		Morph::StatType interpstats;
		Morph::IndexType npix_interp;
		Morph::BoxStats(input_image,input_mask,Nx,Ny,bgbox,
				ImageRejectionMask|TrailMask,ImageExceptionMask,
				interpstats,npix_interp);
// 		assert(npix_interp > 0);
// 		Morph::ClosestValidValues(input_image,input_mask,Nx,Ny,index,
// 					  ImageRejectionMask,ImageExceptionMask,
// 					  vals,dist);
// 		if(dist[0] > 0 && dist[1] > 0){
// 		  Morph::ImageDataType total_dist = dist[0] + dist[1];
// 		  Morph::ImageDataType distm1 = 1.0/total_dist;
// 		  Morph::ImageDataType weight_1 = (total_dist - dist[0]);
// 		  Morph::ImageDataType weight_2 = (total_dist - dist[1]);
// 		  if((vals[0] > star_level) || (vals[1] > star_level)){
// 		    output_mask[index] |= StarMask;
// 		  }
// 		  else{
// 		    output_mask[index] |= InterpMask;
// 		    output_image[index] = (weight_1*vals[0] + weight_2*vals[1])*distm1 + mygasdev(&ranseed)*stats[Image::IMSIGMA];
// 		  }

// 		}
		Morph::ImageDataType baseval = interpstats[Image::IMMEAN];
		output_image[index] = baseval + mygasdev(&ranseed)*stats[Image::IMSIGMA];
		output_mask[index] |= InterpMask;
		if(baseval > star_level)
		  output_mask[index] |= StarMask;
		//		  else {
		//		  }
		//		}
// 		else if(dist[0] > 0){
// 		  if(vals[0] > star_level)
// 		    output_mask[index] |= StarMask;
// 		  else{
// 		    output_mask[index] |= InterpMask;
// 		    output_image[index] = vals[0] + mygasdev(&ranseed)*stats[Image::IMSIGMA];
// 		  }
// 		}
// 		else if(dist[1] > 0){
// 		  if(vals[1] > star_level)
// 		    output_mask[index] |= StarMask;
// 		  else{
// 		    output_image[index] = vals[1] + mygasdev(&ranseed)*stats[Image::IMSIGMA];
// 		    output_mask[index] |= InterpMask;
// 		  }
// 		}
	      }
	      ntrail++;
	    }
	  }
	}
      }
    }
    else {
      for(Morph::IndexType y = box[2];y <= box[3];y++){
	for(Morph::IndexType x = box[0];x <= box[1];x++){
	  Morph::IndexType index = y*Nx + x;
	  // Make sure pixel at index isn't in the rejection mask before doing next set of things
	  if(!(input_mask[index]&TrailRejectionMask) || (input_mask[index]&TrailExceptionMask)){
	    Morph::ImageDataType trail_mean = 0;
	    Morph::ImageDataType trail_min = 0;
	    Morph::ImageDataType trail_max = 0;
	    Morph::ImageDataType trail_sigma = 0;
	    Morph::IndexType trail_npix = 0;
	    std::vector<Morph::IndexType>::iterator selIt = trail_struct.begin();
	    npix_bright = 0;
	    while(selIt != trail_struct.end() && (npix_bright < trail_pix_limit)){
	      Morph::IndexType sindex = index + *selIt++;
	      if(!(input_mask[sindex]&TrailRejectionMask) || (input_mask[sindex]&TrailExceptionMask)){
		if(input_image[sindex] > trail_level)
		  npix_bright++;
		else 
		  npix_bright = 0;
	      }
	    }
	    if(npix_bright >= trail_pix_limit){
	      output_mask[index] |= TrailMask;
	      if(interpolate){
		std::vector<Morph::ImageDataType> vals(4,0);
		std::vector<Morph::IndexType>     dist(4,0);
		Morph::IndexType trailwidth = GetTrailWidth(input_image,input_mask,Nx,Ny,index,
							    ImageRejectionMask|TrailMask,ImageExceptionMask);
		//		assert(trailwidth < 20);
		Morph::BoxType bgbox(4,0);
		bgbox[0] = x - (trailwidth+4);
		bgbox[1] = x + (trailwidth+4);
		bgbox[2] = y - 1;
		bgbox[3] = y + 1;
		if(bgbox[0] < 0) bgbox[0] = 0;
		if(bgbox[1] >= Nx) bgbox[1] = Nx-1;
		if(bgbox[2] < 0) bgbox[2] = 0;
		if(bgbox[3] >= Ny) bgbox[3] = Ny-1;
		Morph::StatType interpstats;
		Morph::IndexType npix_interp;
		Morph::BoxStats(input_image,input_mask,Nx,Ny,bgbox,
				ImageRejectionMask|TrailMask,ImageExceptionMask,
				interpstats,npix_interp);
		//		assert(npix_interp > 0);
		//		Morph::ClosestValidValues(input_image,input_mask,Nx,Ny,index,
		//					  ImageRejectionMask,ImageExceptionMask,
		//					  vals,dist);
		//		if(dist[0] > 0 && dist[1] > 0){
		//		    Morph::ImageDataType total_dist = dist[0] + dist[1];
		//		    Morph::ImageDataType distm1     = 1.0/total_dist;
		//		    Morph::ImageDataType weight_1   = (total_dist - dist[0]);
		//		    Morph::ImageDataType weight_2   = (total_dist - dist[1]);
		//		    if((vals[0] > star_level) || (vals[1] > star_level)){
		//		output_image[index] = (weight_1*vals[0] + weight_2*vals[1])*distm1 + mygasdev(&ranseed)*stats[Image::IMSIGMA];
		Morph::ImageDataType baseval = interpstats[Image::IMMEAN];
		output_image[index] = baseval + mygasdev(&ranseed)*stats[Image::IMSIGMA];
		output_mask[index] |= InterpMask;
		if(baseval > star_level)
		  output_mask[index] |= StarMask;
		//		    }
		//		    else{
		//		    }
		//		}
		//		  else {

		//		  }
		//		}
		//		else if(dist[0] > 0){
		//		  if(vals[0] < star_level){
		//		    output_mask[index] |= InterpMask;
		//		    output_image[index] = vals[0] + mygasdev(&ranseed)*stats[Image::IMSIGMA];
		//		  }
		//		  else
		    //	    if(vals[0] > star_level)
		//		      output_mask[index] |= StarMask;
		//		}
		//		else if(dist[1] > 0){
		//		  if(vals[1] < star_level){
		//		    output_mask[index] |= InterpMask;
		//		    output_image[index] = vals[1] + mygasdev(&ranseed)*stats[Image::IMSIGMA];
		//		  }
		//		  else
		    //		  if(vals[1] > star_level)
		//		    output_mask[index] |= StarMask;
		//		}
	      }
	      ntrail++;
	    }
	  }
	}
      }
    }
  }
  return(ntrail);
}

///
/// \brief Basic structuring element-based bleed trail detection
/// \param input_image Pointer to image data
/// \param output_image Pointer to output image data (only used if interpolation enabled)
/// \param input_mask  Pointer to input mask data
/// \param output_mask  Pointer to output mask data
/// \param Nx    Input number of columns in image
/// \param Ny    Input number of rows in image
/// \param blobs Input rectangular regions in which to detect bleedtrails
/// \param stats Input image statistics for each blob
/// \param trail_struct  Input structuring element for trail direction
/// \param TrailRejectionMask Input bit mask for pixels that should not be used in trail S.E.
/// \param TrailExceptionMask Input bit mask for rejection override in trail S.E.
/// \param ImageRejectionMask Input bit mask for pixels that should not be used in ground S.E.
/// \param ImageExceptionMask Input bit mask for rejection override in ground S.E.
/// \param TrailMask     Input bit mask for setting a bleedtrail pixel
/// \param InterpMask    Input bit mask for identifying interpolated pixel
/// \param StarMask      Input bit mask for identifying star pixel
/// \param trail_pix_limit Indicates how many bright pixels in a row determines trail pixels
/// \param scalefactor  Indicates how many sigma away from mean indicates bright
/// \param interpolate   Input flag for whether to interpolate over bleedtrails
/// \param Ostr Stream object for output
///
/// This function does the following:\n
/// - Loop over every pixel in the input blobs
/// - Evaluate the trail structuring element
/// - If trail_pix_limit bright pixels in a row are found 
///   in the structuring element, then mark current pixel
///   as a bleed trail with TrailMask
/// - If it's a bleed trail, check the valid image values
///   around the pixel to see if they are of star level, if
///   so mark the pixel as a star with StarMask.
/// - If interpolation is enabled, use the local mean (if it exists)
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
    }
    else{
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
	    }
	    else{ 
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
	    bgbox[0] = pixx - 2.0*(trailwidth+2);
	    bgbox[1] = pixx + 2.0*(trailwidth+2);
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
