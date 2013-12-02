#include "ImageMorphology.hh"

int DetectBleedTrails(Morph::ImageDataType *input_image,
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
		      double scalefactor=1.0,
		      bool interpolate = false,
		      std::ostream &Ostr=std::cout);

int DetectBleedTrailsInBoxes(Morph::ImageDataType *image,
			     Morph::ImageDataType *output_image,
			     Morph::MaskDataType *input_mask,
			     Morph::MaskDataType *output_mask,
			     Morph::IndexType Nx,
			     Morph::IndexType Ny,
			     std::vector<Morph::BoxType>   &boxes,
			     std::vector<Morph::StatType>  &box_stats,
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
			     std::ostream &Ostr);

int DetectBleedTrailsInBlobs(Morph::ImageDataType *image,
			     Morph::ImageDataType *output_image,
			     Morph::MaskDataType *input_mask,
			     Morph::MaskDataType *output_mask,
			     Morph::IndexType Nx,
			     Morph::IndexType Ny,
			     std::vector<std::vector<Morph::IndexType> >   &blobs,
			     std::vector<Morph::StatType>  &blob_stats,
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
			     std::ostream &Ostr);


class BleedTrailDetectionFilter : public Morph::ImageFilter
{
protected:
  Morph::StructuringElementType structuring_element;
  Morph::StructuringElementType trimmed_structuring_element;
  Morph::IndexType number_of_pixels_for_trail_detection;
  Morph::ImageDataType trail_detection_level;
  Morph::ImageDataType star_detection_level;
  Morph::ImageDataType background_noise;
  Morph::ImageDataType star_noise;
  Morph::MaskDataType TrailMask;
  Morph::MaskDataType StarMask;
  Morph::MaskDataType InterpolatedMask;
private:
  int yborder_minus;
  int yborder_plus;
  int xborder_minus;
  int xborder_plus;
  long ranseed;
  Morph::IndexType ylimit;
  Morph::IndexType xlimit;
  Morph::IndexType pixx;
  Morph::IndexType pixy;
public:
  virtual void Initialize()
  {
    ImageFilter::Initialize();
    ranseed = -17;
    yborder_minus = structuring_element.MinYOffset();
    yborder_plus  = structuring_element.MaxYOffset();
    xborder_minus = structuring_element.MinXOffset();
    xborder_plus  = structuring_element.MaxXOffset();
    if(yborder_minus < 0)
      yborder_minus *= -1;
    else
      yborder_minus = 0;
    if(xborder_minus < 0)
      xborder_minus *= -1;
    else
      xborder_minus = 0;
    if(xborder_plus < 0)
      xborder_plus = 0;
    if(yborder_plus < 0)
      yborder_plus = 0;
    xlimit = nx - xborder_plus;
    ylimit = NumberOfPixelsInY() - yborder_plus;
  };
  void SetStructuringElement(Morph::StructuringElementType &input_structuring_element)
  {structuring_element = input_structuring_element;};
  void SetTrailDetectionLevel(Morph::ImageDataType level)  { trail_detection_level = level;                 }; 
  void SetStarDetectionLevel(Morph::ImageDataType level)   { star_detection_level = level;                  }; 
  void SetTrailDetectionThreshold(Morph::IndexType thresh) { number_of_pixels_for_trail_detection = thresh; };
  void SetTrailMask(Morph::MaskDataType mask)              { TrailMask = mask;                              };
  void SetStarMask(Morph::MaskDataType mask)               { StarMask  = mask;                              };
  void SetInterpolatedMask(Morph::MaskDataType mask)       { InterpolatedMask  = mask;                      };
  void SetStarNoise(Morph::ImageDataType noise)            { star_noise = noise;                            };
  void SetBackgroundNoise(Morph::ImageDataType noise)      { background_noise = noise;                      };
  inline Morph::ImageDataType StarNoise() const 
  {
    return(mygasdev(const_cast<long *>(&ranseed))*star_noise);
  }
  inline Morph::ImageDataType BackgroundNoise() const 
  {
    return(mygasdev(const_cast<long *>(&ranseed))*star_noise);
  }
  inline unsigned int GetTrailDetectionThreshold() { return(number_of_pixels_for_trail_detection); };
  inline void TrimStructuringElement()
  {
    trimmed_structuring_element.resize(0);
    Morph::StructuringElementType::iterator sei = structuring_element.begin();
    while(sei != structuring_element.end()){
      int xx = pixx + sei->first;
      int yy = pixy + sei->second;
      if(xx < nx && xx >= 0 && yy >= 0 && yy < ny)
	trimmed_structuring_element.push_back(*sei);
      sei++;
    }
  };
  inline bool StructuringElementIntersectsImageBoundary(Morph::IndexType pixel_index)
  {
    pixx = pixel_index%nx;
    pixy = pixel_index/nx;
    return((pixx < xborder_minus) || (pixx > xlimit) ||
	   (pixy < yborder_minus) || (pixy > ylimit));
  };
  inline Morph::StructuringElementType &GetTrimmedStructuringElement(Morph::IndexType pixel_index)
  {
    if(StructuringElementIntersectsImageBoundary(pixel_index)){
      TrimStructuringElement();
      return(trimmed_structuring_element);
    } else {
      return(structuring_element);
    }
  };
  inline void MarkAsTrailPixel(Morph::IndexType pixel_index)        { OutputMaskBuffer[pixel_index] |= TrailMask; };
  inline void MarkAsStarPixel(Morph::IndexType pixel_index)         { OutputMaskBuffer[pixel_index] |= StarMask; };
  inline void MarkAsInterpolatedPixel(Morph::IndexType pixel_index) { OutputMaskBuffer[pixel_index] |= InterpolatedMask; };
  inline Morph::StructuringElementType &GetStructuringElement(Morph::IndexType pixel_index)
  { return(GetTrimmedStructuringElement(pixel_index)); };
  inline bool PixelValueIsTrailLevel(Morph::IndexType pixel_index)
  { return(InputImageBuffer[pixel_index] >= trail_detection_level); };
  int DetectBleedsOnBlobs(Morph::BlobsType &blobs,std::vector<Morph::ImageDataType> &high_levels);
  int DetectStarsOnBlobs(Morph::BlobsType &blobs,std::vector<Morph::ImageDataType> &high_levels);
  int LinearInterpolateOverTrailsInBlobs(Morph::BlobsType &blobs);
};
