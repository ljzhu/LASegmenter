#ifndef imageProcessing_txx_
#define imageProcessing_txx_


#include "imageProcessing.h"

#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkContinuousIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "itkSignedDanielssonDistanceMapImageFilter.h"


#include "vnl/vnl_det.h"


#include <csignal>


namespace douher
{
  /*==================================================*/
  template< typename T > 
  void gradientImage2(typename itk::Image< T, 2 >::Pointer image, \
                      typename itk::Image< T, 2 >::Pointer &gradX, \
                      typename itk::Image< T, 2 >::Pointer &gradY, \
                      unsigned char zeroBoundary)
  {
    typedef itk::Image< T, 2 > InputImageType;

    typedef itk::DerivativeImageFilter<InputImageType, InputImageType> DerivativeImageFilterType;

    typename DerivativeImageFilterType::Pointer filterX = DerivativeImageFilterType::New();
    filterX->SetOrder( 1 );
    filterX->SetDirection( 0 ); // 0 is X dir
    filterX->SetInput( image );
    filterX->Update();
    gradX = filterX->GetOutput();

    // std::cout<<"size is "<<gradX->GetLargestPossibleRegion().GetSize()<<std::endl;


    typename DerivativeImageFilterType::Pointer filterY = DerivativeImageFilterType::New();
    filterY->SetOrder( 1 );
    filterY->SetDirection( 1 ); // 0 is Y dir
    filterY->SetInput( image );
    filterY->Update();
    gradY = filterY->GetOutput();


    if (zeroBoundary == 1 )
      {
        typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorWithIdxType;
        IteratorWithIdxType xIter( gradX, gradX->GetLargestPossibleRegion() );
        IteratorWithIdxType yIter( gradY, gradY->GetLargestPossibleRegion() );

        typename InputImageType::SizeType imSize = image->GetLargestPossibleRegion().GetSize();

        for (xIter.GoToBegin(), yIter.GoToBegin(); !xIter.IsAtEnd(); ++xIter, ++yIter)
          {
            typename InputImageType::IndexType sub = xIter.GetIndex();
            long ix = static_cast< long >(sub[0]);
            long iy = static_cast< long >(sub[1]);

            if (ix == 0 || iy == 0 || ix == imSize[0] - 1 || iy == imSize[1] - 1)
              {
                xIter.Set(0);
                yIter.Set(0);
              }
          }
      }
  }


  /*==================================================*/
  template< typename T > 
  void gradientImage3(typename itk::Image< T, 3 >::Pointer image, \
                      typename itk::Image< T, 3 >::Pointer &gradX, \
                      typename itk::Image< T, 3 >::Pointer &gradY, \
                      typename itk::Image< T, 3 >::Pointer &gradZ, \
                      unsigned char zeroBoundary)
  {
    typedef itk::Image< T, 3 > InputImageType;

    typedef itk::DerivativeImageFilter<InputImageType, InputImageType> DerivativeImageFilterType;

    typename DerivativeImageFilterType::Pointer filterX = DerivativeImageFilterType::New();
    filterX->SetOrder( 1 );
    filterX->SetDirection( 0 ); // 0 is X dir
    filterX->SetInput( image );
    filterX->Update();
    gradX = filterX->GetOutput();


    typename DerivativeImageFilterType::Pointer filterY = DerivativeImageFilterType::New();
    filterY->SetOrder( 1 );
    filterY->SetDirection( 1 ); // 1 is Y dir
    filterY->SetInput( image );
    filterY->Update();
    gradY = filterY->GetOutput();


    typename DerivativeImageFilterType::Pointer filterZ = DerivativeImageFilterType::New();
    filterZ->SetOrder( 1 );
    filterZ->SetDirection( 2 ); // 2 is Z dir
    filterZ->SetInput( image );
    filterZ->Update();
    gradZ = filterZ->GetOutput();

    if ( zeroBoundary == 1 )
      {
        typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorWithIdxType;
        IteratorWithIdxType xIter( gradX, gradX->GetLargestPossibleRegion() );
        IteratorWithIdxType yIter( gradY, gradY->GetLargestPossibleRegion() );
        IteratorWithIdxType zIter( gradZ, gradZ->GetLargestPossibleRegion() );

        typename InputImageType::SizeType imSize = image->GetLargestPossibleRegion().GetSize();

        for (xIter.GoToBegin(), yIter.GoToBegin(), zIter.GoToBegin(); \
             !xIter.IsAtEnd(); ++xIter, ++yIter, ++zIter)
          {
            typename InputImageType::IndexType sub = xIter.GetIndex();
            long ix = static_cast< long >(sub[0]);
            long iy = static_cast< long >(sub[1]);
            long iz = static_cast< long >(sub[2]);

            if (ix == 0 || ix == imSize[0] - 1                    \
                || iy == 0 || iy == imSize[1] - 1                 \
                || iz == 0 || iz == imSize[2] - 1)
              {
                xIter.Set(0);
                yIter.Set(0);
                zIter.Set(0);
              }
          }
      }

  }


  /*************************************************************
   affine transform */
  template< typename T > typename itk::Image< T, 2 >::Pointer affineTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                                               vnl_matrix_fixed< double, 2, 2 > A, \
                                                                               vnl_vector_fixed< double, 2 > translation, \
                                                                               T fillValue, \
                                                                               unsigned char interpType, \
                                                                               long newNx, long newNy)
  {
    const int SpatialDimension_ = 2;


#ifndef NDEBUG
    if (vnl_det(A) < 0.1)
      {
        std::cout<<"Warning, determinent of A < 0.1 \n";
      }
#endif      

    typedef itk::Image< T, SpatialDimension_ > ImageType;

    /* allocate the new image */
    typename ImageType::Pointer newImage = ImageType::New();
    
    // Determine the size of the new image
    if ( -1 == newNx || -1 == newNy )
      {
        // not specific region info, use the original image's
        newImage->SetRegions( image->GetLargestPossibleRegion() );
      }
    else if ( -1 != newNx && -1 != newNy)
      {
        typename ImageType::IndexType start;
        start[0] =   0;  // first index on X
        start[1] =   0;  // first index on Y

        typename ImageType::SizeType size;
        size[0]  = newNx;  // size along X
        size[1]  = newNy;  // size along Y


        typename ImageType::RegionType region;
        region.SetSize( size );
        region.SetIndex( start );

        newImage->SetRegions( region );
      }
    else
      {
        std::cerr<<"Either specify all newNx, newNy and newNz, or assign none. No in-between. abort\n";
        raise(SIGABRT);
      }
    
    newImage->SetSpacing( image->GetSpacing() );
    newImage->Allocate();

    /* fill the new image using the default value */
    newImage->FillBuffer(fillValue);

    /*
      Go over all the pixels (x, y) in the new image, get the value: 
    
      f[ A*(x - cX, y - cY)^T + (transX, transY)^T + (cX, cY)^T ]
    */


    /* get the center of the image */
    typename ImageType::RegionType::SizeType size = newImage->GetLargestPossibleRegion().GetSize();
    vnl_vector_fixed< double, SpatialDimension_ > center;

    typename ImageType::RegionType::SizeType oriSize = image->GetLargestPossibleRegion().GetSize();
    vnl_vector_fixed< double, SpatialDimension_ > oriCenter;

    for (int i = 0; i < SpatialDimension_; ++i)
      {
        center[i] = (double)size[i]/2.0;
        oriCenter[i] = (double)oriSize[i]/2.0;
      }

    //    std::cout<<"center is "<<center<<std::endl;

    /* continuous idx is used to extract between-grid values*/
    typedef itk::ContinuousIndex<double, SpatialDimension_> ContinuousIndexType;

    /* set up the for the old image */
    typename itk::InterpolateImageFunction<ImageType, double>::Pointer interpolator;

    if (interpType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 2)
      {
        // bspline
        typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }


    interpolator->SetInputImage(image);

    /* Iterator goes thru new image: map every point BACK to where it
       comes from, then evaluate the value there and fill the new image.
    */
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorWithIdxType;
    IteratorWithIdxType newImageIter( newImage, newImage->GetLargestPossibleRegion() );
    for (newImageIter.GoToBegin(); !newImageIter.IsAtEnd(); ++newImageIter)
      {
        typename ImageType::IndexType newIdx = newImageIter.GetIndex();

        vnl_vector_fixed< double, SpatialDimension_ > contNewIdx;
        for (long id = 0; id < SpatialDimension_; ++id)
          {
            contNewIdx[id] = newIdx[id];
          }

        // de-center, affine, restore center
        contNewIdx = A*(contNewIdx - center) + translation + oriCenter;

        ContinuousIndexType contNewIdxMapBack;
        for (long id = 0; id < SpatialDimension_; ++id)
          {
            contNewIdxMapBack[id] = contNewIdx[id];
          }


        if (interpolator->IsInsideBuffer(contNewIdxMapBack))
          {
            newImageIter.Set( interpolator->EvaluateAtContinuousIndex( contNewIdxMapBack ) );
          }
      }

    return newImage;
  }

  /*============================================================*/
  template< typename T > typename itk::Image< T, 2 >::Pointer affineTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                                               vnl_matrix< double > A, \
                                                                               vnl_vector< double > translation, \
                                                                               T fillValue, \
                                                                               unsigned char interpType, \
                                                                               long newNx, long newNy)
  {
    if (translation.size() != 2 || A.rows() != 2 || A.columns() != 2)
      {
        std::cerr<<"affine matrix and translation vector are all 2D. abort. \n";
        raise(SIGABRT);
      }

    vnl_matrix_fixed< double, 2, 2 > Af;
    Af(0, 0) = A(0, 0);
    Af(0, 1) = A(0, 1);
    Af(1, 0) = A(1, 0);
    Af(1, 1) = A(1, 1);

    vnl_vector_fixed< double, 2 > Tf;
    Tf[0] = translation[0];
    Tf[1] = translation[1];


    return affineTransform2<T>(image, Af, Tf, fillValue, interpType, newNx, newNy);
  }

  /*============================================================*/
  template< typename T > typename itk::Image< T, 3 >::Pointer affineTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                                               vnl_matrix_fixed< double, 3, 3 > A, \
                                                                               vnl_vector_fixed< double, 3 > translation, \
                                                                               T fillValue, \
                                                                               unsigned char interpType, \
                                                                               long newNx, long newNy, long newNz)
  {
    const int SpatialDimension_ = 3;


#ifndef NDEBUG
    if (vnl_det(A) < 0.1)
      {
        std::cout<<"Warning, determinent of A < 0.1 \n";
      }
#endif      

    typedef itk::Image< T, SpatialDimension_ > ImageType;

    /* allocate the new image */
    typename ImageType::Pointer newImage = ImageType::New();
    // Determine the size of the new image
    if ( -1 == newNx || -1 == newNy || -1 == newNz )
      {
        // not specific region info, use the original image's
        newImage->SetRegions( image->GetLargestPossibleRegion() );
      }
    else if ( -1 != newNx && -1 != newNy && -1 != newNz )
      {
        typename ImageType::IndexType start;
        start[0] =   0;  // first index on X
        start[1] =   0;  // first index on Y
        start[2] =   0;  // first index on Z

        typename ImageType::SizeType size;
        size[0]  = newNx;  // size along X
        size[1]  = newNy;  // size along Y
        size[2]  = newNz;  // size along Z


        typename ImageType::RegionType region;
        region.SetSize( size );
        region.SetIndex( start );

        newImage->SetRegions( region );
      }
    else
      {
        std::cerr<<"Either specify all newNx, newNy and newNz, or assign none. No in-between. abort\n";
        raise(SIGABRT);
      }

    newImage->SetSpacing( image->GetSpacing() );
    newImage->Allocate();

    /* fill the new image using the default value */
    newImage->FillBuffer(fillValue);

    /*
      Go over all the pixels (x, y) in the new image, get the value: 
    
      f[ A*(x - cX, y - cY)^T + (transX, transY)^T + (cX, cY)^T ]
    */


    /* get the center of the image */
    typename ImageType::RegionType::SizeType size = newImage->GetLargestPossibleRegion().GetSize();
    vnl_vector_fixed< double, SpatialDimension_ > center;

    typename ImageType::RegionType::SizeType oriSize = image->GetLargestPossibleRegion().GetSize();
    vnl_vector_fixed< double, SpatialDimension_ > oriCenter;

    for (int i = 0; i < SpatialDimension_; ++i)
      {
        center[i] = (double)size[i]/2.0;
        oriCenter[i] = (double)oriSize[i]/2.0;
      }

    //Debug
    //    std::cout<<"center is "<<center<<std::endl;
    

    /* continuous idx is used to extract between-grid values*/
    typedef itk::ContinuousIndex<double, SpatialDimension_> ContinuousIndexType;

    /* set up the for the old image */
    typename itk::InterpolateImageFunction<ImageType, double>::Pointer interpolator;

    if (interpType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 2)
      {
        // bspline
        typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }


    interpolator->SetInputImage(image);

    /* Iterator goes thru new image: map every point BACK to where it
       comes from, then evaluate the value there and fill the new image.
    */


    //    std::cout<<"A = "<<A<<std::endl;

    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorWithIdxType;
    IteratorWithIdxType newImageIter( newImage, newImage->GetLargestPossibleRegion() );
    for (newImageIter.GoToBegin(); !newImageIter.IsAtEnd(); ++newImageIter)
      {
        typename ImageType::IndexType newIdx = newImageIter.GetIndex();

        vnl_vector_fixed< double, SpatialDimension_ > contNewIdx;
        //        vnl_vector< double > contNewIdx(SpatialDimension_);

        for (long id = 0; id < SpatialDimension_; ++id)
          {
            contNewIdx[id] = newIdx[id];
          }


        // de-center, affine, restore center
        contNewIdx = A*(contNewIdx - center) + translation + oriCenter;

        ContinuousIndexType contNewIdxMapBack;
        for (long id = 0; id < SpatialDimension_; ++id)
          {
            contNewIdxMapBack[id] = contNewIdx[id];
          }


        if (interpolator->IsInsideBuffer(contNewIdxMapBack))
          {
            newImageIter.Set( interpolator->EvaluateAtContinuousIndex( contNewIdxMapBack ) );
          }
      }

    return newImage;    
  }


  /*============================================================*/
  template< typename T > typename itk::Image< T, 3 >::Pointer affineTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                                               vnl_matrix< double > A, \
                                                                               vnl_vector< double > translation, \
                                                                               T fillValue, \
                                                                               unsigned char interpType = 1, \
                                                                               long newNx = -1, long newNy = -1, long newNz = -1)
  {
    if (translation.size() != 3 || A.rows() != 3 || A.columns() != 3)
      {
        std::cerr<<"affine matrix and translation vector are all 3D. abort. \n";
        raise(SIGABRT);
      }

    vnl_matrix_fixed< double, 3, 3 > Af;
    Af(0, 0) = A(0, 0);
    Af(0, 1) = A(0, 1);
    Af(0, 2) = A(0, 2);
    Af(1, 0) = A(1, 0);
    Af(1, 1) = A(1, 1);
    Af(1, 2) = A(1, 2);
    Af(2, 0) = A(2, 0);
    Af(2, 1) = A(2, 1);
    Af(2, 2) = A(2, 2);

    vnl_vector_fixed< double, 3 > Tf;
    Tf[0] = translation[0];
    Tf[1] = translation[1];
    Tf[2] = translation[2];


    return affineTransform3<T>(image, Af, Tf, fillValue, interpType, newNx, newNy, newNz);
  }




  /*============================================================
    SDF */
  template< typename Tin, typename Tout, unsigned int dim >
  typename itk::Image< Tout, dim >::Pointer sdfOfBin( typename itk::Image< Tin, dim >::Pointer binImg)
  {
    typedef itk::Image< Tin, dim > InputImageType;
    typedef itk::Image< Tout, dim > OutputImageType;

    typedef itk::SignedDanielssonDistanceMapImageFilter< InputImageType, OutputImageType >  FilterType;

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( binImg );
    filter->InsideIsPositiveOn();
    filter->Update();

    return filter->GetOutput();
  }

  /*============================================================*/
  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 2 >::Pointer sdfOfBin2D( typename itk::Image< Tin, 2 >::Pointer binImg)
  {
    return sdfOfBin< Tin, Tout, 2 >(binImg);
  }


  /*============================================================*/
  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 3 >::Pointer sdfOfBin3D( typename itk::Image< Tin, 3 >::Pointer binImg)
  {
    return sdfOfBin< Tin, Tout, 3 >(binImg);
  }



  /*============================================================
    tanh */
  template< typename Tin, unsigned int dim >
  typename itk::Image< double, dim >::Pointer tanhOfImage( typename itk::Image< Tin, dim >::Pointer img)
  {
    typedef itk::Image< Tin, dim > InputImageType;
    typedef itk::Image< double, dim > OutputImageType;

    typedef itk::SignedDanielssonDistanceMapImageFilter< InputImageType, OutputImageType >  FilterType;

    typename OutputImageType::Pointer tanhImg = OutputImageType::New();
    tanhImg->SetRegions(img->GetLargestPossibleRegion() );
    tanhImg->CopyInformation(img);
    tanhImg->Allocate();

    typedef itk::ImageRegionIterator< InputImageType > IteratorType;
    typedef itk::ImageRegionIterator< OutputImageType > OutputIteratorType;

    IteratorType imgIter(img, img->GetLargestPossibleRegion());
    OutputIteratorType tanhImgIter(tanhImg, tanhImg->GetLargestPossibleRegion());

    for (imgIter.GoToBegin(), tanhImgIter.GoToBegin();  \
         !imgIter.IsAtEnd();                             \
         ++imgIter, ++tanhImgIter)
      {
        Tin v = imgIter.Get();
        tanhImgIter.Set( tanh(v) );
      }

    return tanhImg;
  }

  /*============================================================*/
  template< typename Tin >
  typename itk::Image< double, 2 >::Pointer tanhOfImage2D( typename itk::Image< Tin, 2 >::Pointer img)
  {
    return tanhOfImage< Tin, 2 >(img);
  }

  /*============================================================*/
  template< typename Tin >
  typename itk::Image< double, 3 >::Pointer tanhOfImage3D( typename itk::Image< Tin, 3 >::Pointer img)
  {
    return tanhOfImage< Tin, 3 >(img);
  }




  /*--------------------------------------------------
    similarityTransform2
  */
  template< typename T > typename itk::Image< T, 2 >::Pointer similarityTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                                                   vnl_vector< double > parameters, \
                                                                                   T fillValue, \
                                                                                   unsigned char interpType)
  {
    if (parameters.size() < 4)
      {
        std::cerr<<"parameter size < 4. It should contain: tranlateX, tranlateY, rotation, scale. exit. \n";
        raise(SIGABRT);
      }

    return similarityTransform2<T>(image, parameters[0], parameters[1], parameters[2], parameters[3], fillValue, interpType);
  }

  /*--------------------------------------------------
    similarityTransform2
  */
  template< typename T > typename itk::Image< T, 2 >::Pointer similarityTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                                                   double translateX, double translateY, \
                                                                                   double rotate, \
                                                                                   double scale, \
                                                                                   T fillValue, \
                                                                                   unsigned char interpType)
  {
    /*
      We denote the original image as a function f(x, y), with origin at
      its center. So let's just say it's in the region of Omega=[-10:1:10]X[-10:1:10]  (21*21 samples)

      Then this function will return the image of 
      g(x, y) := f( s*R*(x, y)^T - (transX, transY)^T )
      where R = (cos, -sin; sin, cos)

      still in the region Omega, still 21*21 samples, sample rate does
      NOT change, it's just the signal itself changes it's frequence by
      s.

      Eg 1, if s > 1, then the result *image* of g looks *smaller*.

      Eg 2, the tree, if there is any, in the result image, would appear
      to be rotated by theta *CCW* (theta>0, and viewing from positive z-dir)

      Eg 3, if transX = 1, transY = 2; then the tree in the new image,
      appear to be moved *right* by 1 and *down* by 2.(assuming the
      postive dir for X is *right*, postive dir for Y is *down*, viewing
      from positive z-dir)


      So, just triverse all points in Omega, at (x, y)^T, what's the
      value of g(x, y)? easy: just f(s*R(x, y)^T - (tranX, transY)^T)

      No "mapping back" blabla which can only confuse me, so many times....
    */

    const long SpatialDimension = 2;

    if (scale <=0)
      {
#ifndef NDEBUG
        std::cerr<<"scale factor should > 0. set to 0.1 \n";
#endif      
        scale = 0.1;
      }



    typedef itk::Image< T, 2 > ImageType;

    /* allocate the new image */
    typename ImageType::Pointer newImage = ImageType::New();
    newImage->SetRegions( image->GetLargestPossibleRegion() );
    newImage->CopyInformation(image);
    newImage->Allocate();

    /* fill the new image using the default value */
    //   typename ImageType::PixelType fillValue;
    //   fillValue = 0.0;
    newImage->FillBuffer(fillValue);

    /*
      Go over all the pixels (x, y) in the new image, get the value: 
    
      f[ s*R*(x - cX, y - cY)^T - (transX, transY)^T + (cX, cY)^T ]
    */


    /* read parameters from passing in*/
    vnl_vector< double > translateCoeff(SpatialDimension);
    translateCoeff[0] = translateX;
    translateCoeff[1] = translateY;

    vnl_matrix< double > rotationMatrix(SpatialDimension, SpatialDimension);
    rotationMatrix(0, 0) = cos(rotate);
    rotationMatrix(1, 0) = sin(rotate);
    rotationMatrix(0, 1) = -rotationMatrix(1, 0);
    rotationMatrix(1, 1) = rotationMatrix(0, 0);
  

    /* get the center of the image */
    typename ImageType::RegionType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    vnl_vector< double > center(SpatialDimension);
    for (long id = 0; id < SpatialDimension; ++id)
      {
        center[id] = (double)size[id]/2.0;
      }

    /* continuous idx is used to extract between-grid values*/
    typedef itk::ContinuousIndex<double, SpatialDimension> ContinuousIndexType;

    /* set up the for the old image */
    typename itk::InterpolateImageFunction<ImageType, double>::Pointer interpolator;

    if (interpType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 2)
      {
        // bspline
        typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }


    interpolator->SetInputImage(image);

    /* Iterator goes thru new image: map every point BACK to where it
       comes from, then evaluate the value there and fill the new image.
    */
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorWithIdxType;
    IteratorWithIdxType newImageIter( newImage, newImage->GetLargestPossibleRegion() );
    for (newImageIter.GoToBegin(); !newImageIter.IsAtEnd(); ++newImageIter)
      {
        typename ImageType::IndexType newIdx = newImageIter.GetIndex();

        vnl_vector< double > contNewIdx(SpatialDimension);
        for (long id = 0; id < SpatialDimension; ++id)
          {
            contNewIdx[id] = newIdx[id];
          }

        // de-center
        contNewIdx -= center;

        // rotation
        contNewIdx.pre_multiply(rotationMatrix);

        // scaling
        contNewIdx *= scale;


        // trans
        contNewIdx -= translateCoeff;

        // restore center
        contNewIdx += center;


        ContinuousIndexType contNewIdxMapBack;
        for (long id = 0; id < SpatialDimension; ++id)
          {
            contNewIdxMapBack[id] = contNewIdx[id];
          }


        if (interpolator->IsInsideBuffer(contNewIdxMapBack))
          {
            newImageIter.Set( interpolator->EvaluateAtContinuousIndex( contNewIdxMapBack ) );
          }
      }

    return newImage;
  }




  /*--------------------------------------------------
    similarityTransform3
  */
  template< typename T > typename itk::Image< T, 3 >::Pointer similarityTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                                                   vnl_vector< double > parameters, \
                                                                                   T fillValue, \
                                                                                   unsigned char interpType)
  {
    if (parameters.size() < 7)
      {
        std::cerr<<"parameter size < 7. It should contain: tranlateX, tranlateY, tranlateZ, rotationAroundX, rotationAroundY, rotationAroundZ, scale. exit. \n";
        raise(SIGABRT);
      }

    return similarityTransform3<T>(image, parameters[0], parameters[1], parameters[2], \
                                   parameters[3], parameters[4], parameters[5], \
                                   parameters[6], \
                                   fillValue,	\
                                   interpType);
  }

  /*--------------------------------------------------
    similarityTransform3
  */
  template< typename T > typename itk::Image< T, 3 >::Pointer similarityTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                                                   double translateX, \
                                                                                   double translateY, \
                                                                                   double translateZ, \
                                                                                   double rotateAroundX, \
                                                                                   double rotateAroundY, \
                                                                                   double rotateAroundZ, \
                                                                                   double scale, \
                                                                                   T fillValue, \
                                                                                   unsigned char interpType)
  {
    /*
      We denote the original image as a function f(x, y, z), with origin at
      its center. So let's just say it's in the region of Omega=[-10:1:10]^3  (21^3 samples)

      Then this function will return the image of 
      g(x, y, z) := f( s*Rx*Ry*Rz*(x, y, z)^T - (transX, transY, transZ)^T )
      where 
      Rx = (1, 0, 0; 0, cosX, -sinX; 0, sinX, cosX)
      Ry = (cosY, 0, sinY; 0, 1, 0; -sinY, 0, cosY)
      Rz = (cosZ, -sinZ, 0; sinZ, cosZ, 0; 0, 0, 1)

      still in the region Omega, still 21*21 samples, sample rate does
      NOT change, it's just the signal itself changes it's frequence by
      s.

      Eg 1, if s > 1, then the result *image* of g looks *smaller*.

      Eg 2, the tree, if there is any, in the result image, would appear
      to be rotated by theta *CCW* around z-dir (if thetaZ>0, and
      viewing from positive z-dir), same for x and y.


      So, just triverse all points in Omega, at (x, y, z)^T, what's the
      value of g(x, y, z)? easy: just f(s*Rx*Ry*Rz(x, y, z)^T - (tranX, transY)^T)

      No "mapping back" blabla which can only confuse me, so many times....
    */

    const long SpatialDimension = 3;


    if (scale <=0)
      {
#ifndef NDEBUG
        std::cerr<<"scale factor should > 0. set to 0.1 \n";
#endif      
        scale = 0.1;
      }


    typedef itk::Image< T, 3 > ImageType;

    /* allocate the new image */
    typename ImageType::Pointer newImage = ImageType::New();
    newImage->SetRegions( image->GetLargestPossibleRegion() );
    newImage->CopyInformation(image);
    newImage->Allocate();

    /* fill the new image using the default value */
    //   typename ImageType::PixelType fillValue;
    //   fillValue = 0.0;
    newImage->FillBuffer(fillValue);

    /*
      Go over all the pixels (x, y, z) in the new image, get the value: 
    
      f[ s*Rx*Ry*Rz*(x - cX, y - cY, z - cZ)^T - (transX, transY, transZ)^T + (cX, cY, cZ)^T ]
    */



    /* read parameters from passing in*/
    vnl_vector< double > translateCoeff(SpatialDimension);
    translateCoeff[0] = translateX;
    translateCoeff[1] = translateY;
    translateCoeff[2] = translateZ;

    vnl_matrix< double >  RX(SpatialDimension, SpatialDimension, 0.0);
    RX(0, 0) = 1.0;
    RX(1, 1) = cos(rotateAroundX);
    RX(2, 1) = sin(rotateAroundX);
    RX(1, 2) = -RX(2, 1);
    RX(2, 2) = RX(1, 1);

    vnl_matrix< double > RY(SpatialDimension, SpatialDimension, 0.0);
    RY(1, 1) = 1.0;
    RY(0, 0) = cos(rotateAroundY);
    RY(0, 2) = sin(rotateAroundY);
    RY(2, 0) = -RY(0, 2);
    RY(2, 2) = RY(0, 0);

    vnl_matrix< double > RZ(SpatialDimension, SpatialDimension, 0.0);
    RZ(2, 2) = 1.0;
    RZ(0, 0) = cos(rotateAroundZ);
    RZ(1, 0) = sin(rotateAroundZ);
    RZ(0, 1) = -RZ(1, 0);
    RZ(1, 1) = RZ(0, 0);


    //    vnl_matrix< double > R = (RX*RY)*RZ;
    //    std::cout<<"R = "<<R<<std::endl;



    /* get the center of the image */
    typename ImageType::RegionType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    vnl_vector< double > center(SpatialDimension);
    for (long id = 0; id < SpatialDimension; ++id)
      {
        center[id] = (double)size[id]/2.0;
      }

    /* continuous idx is used to extract between-grid values*/
    typedef itk::ContinuousIndex<double, SpatialDimension> ContinuousIndexType;

    /* set up the linear interpolator for the old image */
    typename itk::InterpolateImageFunction<ImageType, double>::Pointer interpolator;

    if (interpType == 0)
      {
        // NN interpolation
#ifndef NDEBUG
        std::cerr<<"NN interpolation\n";
#endif
        typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 1)
      {
        // linear 
#ifndef NDEBUG
        std::cerr<<"linear interpolation\n";
#endif
        typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 2)
      {
#ifndef NDEBUG
        std::cerr<<"bspline interpolation \n";
#endif
        typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }

    interpolator->SetInputImage(image);

    /* Iterator goes thru new image: map every point BACK to where it
       comes from, then evaluate the value there and fill the new image.
    */
    typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorWithIdxType;
    IteratorWithIdxType newImageIter( newImage, newImage->GetLargestPossibleRegion() );
    for (newImageIter.GoToBegin(); !newImageIter.IsAtEnd(); ++newImageIter)
      {
        typename ImageType::IndexType newIdx = newImageIter.GetIndex();

        vnl_vector< double > contNewIdx(SpatialDimension);

        for (long id = 0; id < SpatialDimension; ++id)
          {
            contNewIdx[id] = newIdx[id];
          }

        // de-center
        contNewIdx -= center;

        // rotation
        contNewIdx.pre_multiply(RZ);
        contNewIdx.pre_multiply(RY);
        contNewIdx.pre_multiply(RX);

        // scaling
        contNewIdx *= scale;

        // trans
        contNewIdx -= translateCoeff;

        // restore center
        contNewIdx += center;


        ContinuousIndexType contNewIdxMapBack;
        for (long id = 0; id < SpatialDimension; ++id)
          {
            contNewIdxMapBack[id] = contNewIdx[id];
          }


        if (interpolator->IsInsideBuffer(contNewIdxMapBack))
          {
            newImageIter.Set( interpolator->EvaluateAtContinuousIndex( contNewIdxMapBack ) );
          }
      }

    return newImage;
  }



  /*==================================================
     threshold  */
  template< typename T, unsigned int dim > 
  typename itk::Image< unsigned char, dim >::Pointer thresholdImage(typename itk::Image< T, dim >::Pointer image, T thld)
  {
    typedef itk::Image< unsigned char, dim > UcharImageType;
    typedef itk::Image< T, dim > InputImageType;

    /* allocate the new image */
    typename UcharImageType::Pointer resultImg = UcharImageType::New();
    resultImg->SetRegions( image->GetLargestPossibleRegion() );
    resultImg->CopyInformation(image);
    resultImg->Allocate();


    typedef itk::ImageRegionIterator<InputImageType> InputIteratorType;
    InputIteratorType inputImIter(image, image->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<UcharImageType> UcharIteratorType;
    UcharIteratorType outputImIter(resultImg, resultImg->GetLargestPossibleRegion());

    for (inputImIter.GoToBegin(), outputImIter.GoToBegin(); \
         !inputImIter.IsAtEnd(); \
         ++inputImIter, ++outputImIter)
      {
        T v = inputImIter.Get();
        outputImIter.Set( (unsigned char)(v > thld));
      }

    return resultImg;
  }

  template< typename T > 
  typename itk::Image< unsigned char, 2 >::Pointer thresholdImage2(typename itk::Image< T, 2 >::Pointer image, T thld)
  {
    return thresholdImage<T, 2>(image, thld);
  }

  template< typename T > 
  typename itk::Image< unsigned char, 3 >::Pointer thresholdImage3(typename itk::Image< T, 3 >::Pointer image, T thld)
  {
    return thresholdImage<T, 3>(image, thld);
  }


  /* ============================================================
     off grid image value   */
  template< typename T, unsigned int dim > T offGridImageValue(typename itk::Image< T, dim >::Pointer image, \
                                                               itk::ContinuousIndex<double, dim> offGridIndex, \
                                                               unsigned char interpType)
  {
    typedef itk::Image< T, dim > ImageType;
    T offGridValue = 0;

    /* set up the for the old image */
    typename itk::InterpolateImageFunction<ImageType, double>::Pointer interpolator;

    if (interpType == 0)
      {
        // NN interpolation
        typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 1)
      {
        // linear 
        typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }
    else if(interpType == 2)
      {
#ifndef NDEBUG
        std::cerr<<"spline interpolation \n";
#endif
        typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> InterpolatorType;
        interpolator = InterpolatorType::New();
      }

    interpolator->SetInputImage(image);

    if (interpolator->IsInsideBuffer(offGridIndex))
      {
        offGridValue = interpolator->EvaluateAtContinuousIndex( offGridIndex ) ;
      }


    return offGridValue;
  }

  /*============================================================
    offGridImageValue 2D */
  template< typename T > T offGridImageValue(typename itk::Image< T, 2 >::Pointer image, \
                                             double x, double y, \
                                             unsigned char interpType)
  {
    typedef itk::ContinuousIndex<double, 2> ContinuousIndexType;

    ContinuousIndexType offGridIndex;
    offGridIndex[0] = x;
    offGridIndex[1] = y;
  

    return offGridImageValue<T, 2>(image, offGridIndex, interpType);
  }


  /*============================================================
    offGridImageValue 3D */
  template< typename T > T offGridImageValue(typename itk::Image< T, 3 >::Pointer image, \
                                             double x, double y, double z, \
                                             unsigned char interpType)
  {
    typedef itk::ContinuousIndex<double, 3> ContinuousIndexType;

    ContinuousIndexType offGridIndex;
    offGridIndex[0] = x;
    offGridIndex[1] = y;
    offGridIndex[2] = z;
  

    return offGridImageValue<T, 3>(image, offGridIndex, interpType);
  }



  /*============================================================
    Intensity linear transform */
  template< typename Tin, typename Tout, unsigned int dim >
  typename itk::Image< Tout, dim >::Pointer intensityLinearTransform( typename itk::Image< Tin, dim >::Pointer img, \
                                                                      double a, double b)
  {
    typedef itk::Image< Tin, dim > InputImageType;
    typedef itk::Image< double, dim > OutputImageType;

    typename OutputImageType::Pointer outputImg = OutputImageType::New();
    outputImg->SetRegions(img->GetLargestPossibleRegion() );
    outputImg->CopyInformation(img);
    outputImg->Allocate();

    typedef itk::ImageRegionIterator< InputImageType > IteratorType;
    typedef itk::ImageRegionIterator< OutputImageType > OutputIteratorType;

    IteratorType imgIter(img, img->GetLargestPossibleRegion());
    OutputIteratorType outputImgIter(outputImg, outputImg->GetLargestPossibleRegion());

    for (imgIter.GoToBegin(), outputImgIter.GoToBegin();  \
         !imgIter.IsAtEnd();                             \
         ++imgIter, ++outputImgIter)
      {
        Tin v = imgIter.Get();

        outputImgIter.Set( a*v + b ); 
      }

    return outputImg;    
  }


  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 2 >::Pointer intensityLinearTransform2D( typename itk::Image< Tin, 2 >::Pointer img, \
                                                                      double a, double b)
  {
    return intensityLinearTransform< Tin, Tout, 2 >(img, a, b);
  }

  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 3 >::Pointer intensityLinearTransform3D( typename itk::Image< Tin, 3 >::Pointer img, \
                                                                      double a, double b)
  {
    return intensityLinearTransform< Tin, Tout, 3 >(img, a, b);
  }




} // namespace

#endif
