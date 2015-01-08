#ifndef utilities_hxx_
#define utilities_hxx_

// std
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <complex>
#include <algorithm>

// gnu
//#include <dirent.h>
//#include <sys/types.h>

// itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImportImageFilter.h"
#include "itkOrientImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
//#include "itkIdentityTransform.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include  "itkXorImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
//#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkAddImageFilter.h"

#include <csignal>

namespace LASeg
{
    template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName);
  /**
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName)
  {
    typedef itk::ImageFileReader< itkImage_t > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename itkImage_t::Pointer image;

    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl;
        std::cerr<< err << std::endl;
        raise(SIGABRT);
      }

    return image;
  }

  /**
   * readImage in an oriented way
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer ReadImageOriented(const char *fileName)
  {
      typedef itk::ImageFileReader< itkImage_t > ImageReaderType;
      typename itkImage_t::Pointer rval;

       itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
       typename ImageReaderType::Pointer fileReader = ImageReaderType::New();
//       fileReader->SetImageIO(io);
       fileReader->SetFileName(fileName);

       try {
           std::cout << "start reading...\n";
           fileReader->Update();
           rval = fileReader->GetOutput();

           typename itk::OrientImageFilter<itkImage_t,itkImage_t>::Pointer orienter =
                                                 itk::OrientImageFilter<itkImage_t,itkImage_t>::New();
           orienter->UseImageDirectionOn();
           orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
           orienter->SetInput(rval);
           orienter->Update();
           rval = orienter->GetOutput();
       }
        catch ( itk::ExceptionObject &err) {
            std::cerr<< "ExceptionObject caught !" << std::endl;
            std::cerr<< err << std::endl;
            raise(SIGABRT);
        }

       return rval;
  }


  /**
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName)
  {
    typedef itk::ImageFileWriter< itkImage_t > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        raise(SIGABRT);
      }
  }
template<typename ImSrcType>
float* ExtractCubicNonZeroRegion(typename ImSrcType::Pointer im, int& dim) {

    typedef itk::ImageMaskSpatialObject< 3 > ImageMaskSpatialObjectType;

    typename ImageMaskSpatialObjectType::Pointer imageMaskSpatialObject = ImageMaskSpatialObjectType::New();
    imageMaskSpatialObject->SetImage( im );

    typename ImSrcType::RegionType boundingBox;
    typename ImSrcType::SizeType size;
    typename ImSrcType::IndexType start,idx;

    boundingBox = imageMaskSpatialObject->GetAxisAlignedBoundingBoxRegion();
    size = boundingBox.GetSize();
    start = boundingBox.GetIndex();

    dim = std::max(std::max(size[0], size[1]), size[2]);
    long total, widLine, areaPlane;
    total = dim*dim*dim;
    float* pData = new float[total];
    std::fill(pData, pData+total, 0.0);

    unsigned int i,j,k;
    widLine = size[0];
    areaPlane = size[0]*size[1];
    for(i = 0; i < size[0]; i++) {
        for(j = 0; j < size[1]; j++) {
            for(k = 0; k < size[2]; k++) {
                idx[0] = i + start[0];
                idx[1] = j + start[1];
                idx[2] = k + start[2];
                pData[i + dim*(j + k*dim)] = im->GetPixel(idx);
            }
        }
    }

    return pData;
}

  /**
   * Read a series of images.
   */
  template< typename inputPixel_t, typename outputPixel_t >
  typename itk::Image<outputPixel_t, 3>::Pointer
  castItkImage( typename itk::Image<inputPixel_t, 3>::Pointer inputImage )
  {
    typedef itk::Image<inputPixel_t, 3> inputImage_t;
    typedef itk::Image<outputPixel_t, 3> outputImage_t;

    typedef itk::CastImageFilter< inputImage_t, outputImage_t > itkCastFilter_t;

    typename itkCastFilter_t::Pointer caster = itkCastFilter_t::New();
    caster->SetInput( inputImage );
    caster->Update();


    return caster->GetOutput();
  }


  /**
   * Read a series of images.
   */
  template< typename itkImage_t >
  std::vector< typename itkImage_t::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList )
  {
    typedef typename itkImage_t::Pointer itkImagePointer_t;
    typedef std::vector< itkImagePointer_t > itkImageList_t;
    typedef itk::ImageFileReader< itkImage_t > itkImageReader_t;


    itkImageList_t imageSeries;

    int n = imageNameList.size();
    for (int i = 0; i < n; ++i)
      {
        std::string thisName = imageNameList[i];

        typename itkImageReader_t::Pointer reader = itkImageReader_t::New();
        reader->SetFileName(thisName);

        itkImagePointer_t img;

        try
          {
            reader->Update();
            img = reader->GetOutput();
          }
        catch ( itk::ExceptionObject &err)
          {
            std::cerr<< "ExceptionObject caught !" << std::endl;
            std::cerr<< err << std::endl;
            raise(SIGABRT);
          }


        imageSeries.push_back(img);
      }

    return imageSeries;
  }

  /*============================================================
   * readTextLineToListOfString
   */
  template<typename TNull>
  std::vector< std::string > readTextLineToListOfString(const char* textFileName)
  {
    /* The file MUST end with an empty line, then each line will be
       stored as an element in the returned vector object. */


    // Here speed and space is not a big issue, so just let it copy and return stuff.
    std::vector< std::string > listOfStrings;

    std::ifstream f(textFileName);
    std::string thisLine;

    if (f.good())
      {
        while( std::getline(f, thisLine))
          {
            if(!thisLine.empty())
                listOfStrings.push_back(thisLine);
            else break;
          }
      }
    else
      {
        std::cerr<<"Error: can not open file:"<<textFileName<<std::endl;
        raise(SIGABRT);
      }

    f.close();

    return listOfStrings;
  }


   /**
        Compute XOR image for extracting myocardium wall
        NOTE: image resampling is needed in case the input images have different SIZEs!
    */
  template<typename inputType1, typename inputType2, typename outputType>
    typename outputType::Pointer
  ExtractXORImage(typename inputType1::Pointer input1, typename inputType2::Pointer input2) {

            typedef itk::XorImageFilter <outputType> XorImageFilterType;
            typename XorImageFilterType::Pointer xorFilter = XorImageFilterType::New();
            xorFilter->SetInput1(input1);
            xorFilter->SetInput2(input2);
            xorFilter->Update();

            return  xorFilter->GetOutput();
  }
  template<typename image_t>
  double getVol(typename image_t::Pointer img, typename image_t::PixelType thld = 0)
  {
    typedef itk::ImageRegionConstIterator<image_t> ImageRegionConstIterator_t;
    ImageRegionConstIterator_t it(img, img->GetLargestPossibleRegion() );

    double cell = (img->GetSpacing()[0])*(img->GetSpacing()[1])*(img->GetSpacing()[2]);

    double v = 0.0;

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        typename image_t::PixelType f = it.Get();

        v += (f>thld?cell:0.0);
      }

    return v;
  }



  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  thld3(typename input_image_t::Pointer input,                            \
        typename input_image_t::PixelType lowerT,                         \
        typename input_image_t::PixelType upperT, \
        typename output_image_t::PixelType insideValue = 1,               \
        typename output_image_t::PixelType outsideValue = 0)
  {
    /**
     * O(x) :=    I(x) \in [lowerT, upperT] ? insideValue : outsideValue
     */

    //tst
    //   std::cout<<lowerT<<std::endl;
    //   std::cout<<upperT<<std::endl;

    typedef itk::BinaryThresholdImageFilter<input_image_t, output_image_t> binaryThresholdImageFilter_t;

    typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
    thlder->SetInput(input);
    thlder->SetInsideValue(insideValue);
    thlder->SetOutsideValue(outsideValue);
    thlder->SetUpperThreshold(upperT);
    thlder->SetLowerThreshold(lowerT);
    thlder->Update();

    return thlder->GetOutput();
  }


  /**
   * Post process probabiliry map: Multiply by 200, then convert to
   * uchar image. This will make the final result to be easier
   * thresholded by Slicer.
   */
  template<typename image_t>
  typename itk::Image<unsigned char, 3>::Pointer
  postProcessProbabilityMap(typename image_t::Pointer probMap, typename image_t::PixelType c)
  {
    typedef itk::MultiplyImageFilter< image_t, typename image_t::PixelType, itk::Image<unsigned char, 3> > \
      itkMultiplyImageFilter_t;

    typename itkMultiplyImageFilter_t::Pointer mul = itkMultiplyImageFilter_t::New();
    mul->SetInput(probMap);
    mul->SetConstant(c);
    mul->Update();

    return mul->GetOutput();
  }

  template<typename InputImageType, typename OutputImageType>
  typename OutputImageType::Pointer
  ResampleImage(typename InputImageType::Pointer imIn, double* sampleSpacing, int imTypeFlag) {

      typename InputImageType::SizeType inSize;
      typename InputImageType::SpacingType inSpacing;
      typename OutputImageType::SizeType outSize;
      typename OutputImageType::SpacingType outSpacing;

      inSize = imIn->GetLargestPossibleRegion().GetSize();
      inSpacing = imIn->GetSpacing();
      outSpacing = sampleSpacing;

      for(unsigned int i = 0; i < 3; i++) {
          outSize[i] = std::floor(inSpacing[i]*inSize[i]/outSpacing[i]+0.5);
      }

      if(imTypeFlag) {
          typedef itk::ResampleImageFilter< InputImageType, OutputImageType >    ResampleFilterType;
          typedef itk::LinearInterpolateImageFunction< InputImageType, double >  InterpolatorType;

          typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
          typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

          resampler->SetInterpolator( interpolator );
          resampler->SetInput( imIn );
          resampler->SetSize( outSize);
          resampler->SetOutputOrigin(  imIn->GetOrigin() );
          resampler->SetOutputSpacing(sampleSpacing);
          resampler->SetOutputDirection( imIn->GetDirection() );
          resampler->SetDefaultPixelValue( 0 );
          resampler->Update();

          return resampler->GetOutput();
      }
      else {
            typedef itk::ResampleImageFilter< InputImageType, OutputImageType >    ResampleFilterType;
            typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, double >  InterpolatorType;
            typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
            typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

            resampler->SetInterpolator( interpolator );
            resampler->SetInput( imIn );
            resampler->SetSize( outSize);
            resampler->SetOutputOrigin(  imIn->GetOrigin() );
            resampler->SetOutputSpacing(sampleSpacing);
            resampler->SetOutputDirection( imIn->GetDirection() );
            resampler->SetDefaultPixelValue( 0 );
            resampler->Update();

            return resampler->GetOutput();

      }
  }

  template<typename ImType>
  void FindMaskCentroid(typename ImType::Pointer imMask, double& cx, double& cy, double& cz) {

      typedef typename ImType::IndexType TIndex;
      typename ImType::SizeType size;
      long ix, iy, iz, vol;

      size = imMask->GetLargestPossibleRegion().GetSize();

      cx = 0;
      cy = 0;
      cz = 0;
      vol = 0;
      for (ix = 0; ix < size[0]; ++ix) {
          for (iy = 0; iy < size[1]; ++iy) {
              for (iz = 0; iz < size[2]; ++iz) {

                  TIndex idx = {{ix, iy, iz}};
                  if(imMask->GetPixel(idx) > 0) {
                      cx += double(ix);
                      cy += double(iy);
                      cz += double(iz);
                      vol++;
                  }
               }
           }
     }

      if(vol > 0) {
          cx /= vol;
          cy /= vol;
          cz /= vol;
      }
  }

  template<typename ImSrcType>
  typename ImSrcType::Pointer
  AnisotropicSmoothing(typename ImSrcType::Pointer imSrc) {

      typedef itk::GradientAnisotropicDiffusionImageFilter< ImSrcType, ImSrcType> FilterType;
      typename FilterType::Pointer filter = FilterType::New();

      filter->SetInput(imSrc);
      filter->SetNumberOfIterations(1);
      filter->SetTimeStep( 0.0390625);
      filter->SetConductanceParameter(1.0);
      filter->Update();

      return filter->GetOutput();
  }

  template<typename ImSrcType, typename ImMaskType>
  void
  MaskOutVolume(typename ImSrcType::Pointer imSrc, typename ImMaskType::Pointer imLab, int maskVal = -1000) {

      typedef itk::ImageRegionConstIterator<ImMaskType> ItConst;
      typedef itk::ImageRegionIterator<ImSrcType> ItSrc;

      ItConst itConst(imLab, imLab->GetLargestPossibleRegion() );
      ItSrc itSrc(imSrc, imSrc->GetLargestPossibleRegion() );

      for(itConst.GoToBegin(), itSrc.GoToBegin(); !itConst.IsAtEnd(); ++itConst, ++itSrc) {

          if(itConst.Get() > 0) {
              itSrc.Set(maskVal);
          }
      }

      LASeg::writeImage<ImSrcType>(imSrc, "../Results/Test/maskout.nrrd");
  }

  template<typename ImSrcType>
  typename ImSrcType::Pointer
  HomogeneityEnhancement(typename ImSrcType::Pointer imSrc) {

      // Homogeneity smoothing
      typename ImSrcType::Pointer imHomo;
      imHomo = LASeg::AnisotropicSmoothing<ImSrcType>(imSrc);

    //  Find gradients
      double variance = 2.0;
      double upperThreshold = 0.0;
      double lowerThreshold = 0.0;
      typedef itk::CannyEdgeDetectionImageFilter <ImSrcType, ImSrcType>  CannyEdgeDetectionImageFilterType;

      typename CannyEdgeDetectionImageFilterType::Pointer cannyFilter  = CannyEdgeDetectionImageFilterType::New();
      cannyFilter->SetInput(imHomo);
      cannyFilter->SetVariance( variance );
      cannyFilter->SetUpperThreshold( upperThreshold );
      cannyFilter->SetLowerThreshold( lowerThreshold );

      cannyFilter->Update();

      // Add two image together
      typedef itk::AddImageFilter <ImSrcType, ImSrcType> AddImageFilterType;

      typename  AddImageFilterType::Pointer addFilter = AddImageFilterType::New ();
      addFilter->SetInput1(imHomo);
      addFilter->SetInput2(cannyFilter->GetOutput());
      addFilter->Update();

      return addFilter->GetOutput();
  }

}// douher

#endif
