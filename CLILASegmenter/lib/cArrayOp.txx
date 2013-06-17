#ifndef cArrayOp_txx_
#define cArrayOp_txx_

#include "cArrayOp.h"

#include "shuxue.h"

#include <iostream>

#include "itkVectorImage.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matlab_filewrite.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkContinuousIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

#include <cassert>
#include <csignal>


namespace douher
{
  /* ============================================================
     itkImageToVnlVector  */
  template< typename T, unsigned int dim > vnl_vector< T > itkImageToVnlVector(typename itk::Image< T, dim >::Pointer img)
  {
    typedef itk::Image< T, dim > ImageType;

    typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();

    long totalLen = 1;
    for (unsigned long i = 0; i < dim; ++i)
      {
        totalLen *= size[i];
      }

    //  std::cout<<"totalLen = "<<totalLen<<std::endl;

    typedef itk::ImageRegionIterator<ImageType> IteratorType;
    IteratorType imIter(img, img->GetLargestPossibleRegion());

    long longVectorIdx = 0;
    vnl_vector< T > vct(totalLen);

    for (imIter.GoToBegin(); !imIter.IsAtEnd(); ++imIter)
      {
        vct[longVectorIdx++] = imIter.Get();
      }

    return vct;
  }


  /* ============================================================
     itkImageToVnlVector, const pointer version 
  template< typename T, unsigned int dim > vnl_vector< T > itkImageToVnlVector(typename itk::Image< T, dim >::ConstPointer img)
  {
    typedef itk::Image< T, dim > ImageType;

    typename ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();

    long totalLen = 1;
    for (unsigned long i = 0; i < dim; ++i)
      {
        totalLen *= size[i];
      }

    //  std::cout<<"totalLen = "<<totalLen<<std::endl;

    typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType imIter(img, img->GetLargestPossibleRegion());

    long longVectorIdx = 0;
    vnl_vector< T > vct(totalLen);

    for (imIter.GoToBegin(); !imIter.IsAtEnd(); ++imIter)
      {
        vct[longVectorIdx++] = imIter.Get();
      }

    return vct;
  } */




  /*============================================================
     vnlVectorToItkImage  */
  template< typename T, unsigned int dim > typename itk::Image< T, dim >::Pointer vnlVectorToItkImage(const vnl_vector< T >& vct, \
                                                                                                      long imSize[])
  {
    typedef itk::Image< T, dim > ImageType;
    typename ImageType::Pointer image = ImageType::New();

    typename ImageType::IndexType start;
    start.Fill(0);

    typename ImageType::SizeType size;
    for (unsigned long i = 0; i < dim; ++i)
      {
        size[i] = imSize[i];
      }

    typename ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();

    typedef itk::ImageRegionIterator<ImageType> ImageIteratorType;
    ImageIteratorType Iter(image, image->GetLargestPossibleRegion());

    long longVectorIdx = 0;
    for (Iter.GoToBegin(); !Iter.IsAtEnd(); ++Iter)
      {
        Iter.Set( vct[longVectorIdx++] );
      }

    return image;
  }

  /* --------------------------------------------------
     readImage2
  */
  template< typename T > typename itk::Image< T, 2 >::Pointer readImage2(const char *fileName)
  {
    return readImage<T, 2>(fileName);
  }

  /* --------------------------------------------------
     readImage3
  */
  template< typename T > typename itk::Image< T, 3 >::Pointer readImage3(const char *fileName)
  {
    return readImage<T, 3>(fileName);
  }


  /* --------------------------------------------------
     readImage
  */
  template< typename T, unsigned int dim > typename itk::Image< T, dim >::Pointer readImage(const char *fileName)
  {
    typedef itk::Image< T, dim > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename ImageType::Pointer image;
    
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


  /*--------------------------------------------------
    writeImage2
  */
  template< typename T > void writeImage2(typename itk::Image< T, 2 >::Pointer img, const char *fileName)
  {
    writeImage<T, 2>(img, fileName);
  }


  /*--------------------------------------------------
    writeImage3
  */
  template< typename T > void writeImage3(typename itk::Image< T, 3 >::Pointer img, const char *fileName)
  {
    writeImage<T, 3>(img, fileName);
  }


  /*--------------------------------------------------
    writeImage
  */
  template< typename T, unsigned int dim > void writeImage(typename itk::Image< T, dim >::Pointer img, const char *fileName)
  {
    typedef itk::Image< T, dim > ImageType;
    typedef itk::ImageFileWriter< ImageType > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);

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


  /*============================================================
    readImageToArray2  */
  template< typename T > 
  typename cArray2D< T >::PointerType readImageToArray2(const char *fileName)
  {
    typedef itk::Image< T, 2 > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename ImageType::Pointer image;
    
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

    return itkImageToArray2< T >(image);
  }



  /*============================================================
    readImageToArray3  */
  template< typename T > 
  typename cArray3D< T >::PointerType readImageToArray3(const char *fileName)
  {
    typedef itk::Image< T, 3 > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename ImageType::Pointer image;
    
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

  
    return itkImageToArray3< T >(image);
  }



  /*============================================================
    readImageToArray4  */
  template< typename T > 
  typename cArray4D< T >::PointerType readImageToArray4(const char *fileName)
  {
    typedef itk::VectorImage< T, 3 > VectorImageType;
    typedef itk::ImageFileReader< VectorImageType > VectorImageReaderType;
    typename VectorImageReaderType::Pointer vectorImgReader = VectorImageReaderType::New();
    vectorImgReader->SetFileName(fileName);

    typename VectorImageType::Pointer vctImage;

    try
      {
        vectorImgReader->Update();
        vctImage = vectorImgReader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl; 
        std::cerr<< err << std::endl; 
        raise(SIGABRT);
      }

    //  std::cerr<<"dim is: "<<vctImage->GetImageDimension()<<std::endl;
    typename VectorImageType::RegionType::SizeType dataSz = (vctImage->GetLargestPossibleRegion()).GetSize();

    //   std::cerr<<"Size of each dim is: ";
    //   long imgDim = vctImage->GetImageDimension();
    //   for (long id = 0; id < imgDim; ++id)
    //     {
    //       std::cerr<<dataSz[id]<<"  ";
    //     }
    //   std::cerr<<std::endl;

    //  std::cerr<<"vector length is: "<<vctImage->GetVectorLength()<<std::endl;
    //  long vctLen = vctImage->GetVectorLength();

    long n4 = vctImage->GetVectorLength();
    typename cArray4D< T >::PointerType array4( new cArray4D< T >(dataSz[0], dataSz[1], dataSz[2], n4) );

    /* port the DWI image to a cArray4D structure */
    typedef itk::ImageRegionIterator< VectorImageType > VectorImageIteratorType;
    VectorImageIteratorType imIter(vctImage, vctImage->GetLargestPossibleRegion() );


    imIter.GoToBegin();
    for (long iz = 0; iz < dataSz[2]; ++iz)
      {
        for (long iy = 0; iy < dataSz[1]; ++iy)
          {
            for (long ix = 0; ix < dataSz[0]; ++ix)
              {
                /*
                  NOTE!!! the order: iz changes slowest while ix changes
                  fastest MATTERS! 

                  Since this is using together with iterator of
                  itkImage. There the iterator first goes thru x-dim,
                  then y-dim, then z-dim. So the fastest changing
                  dimension should be consistant.
                */

                typename VectorImageType::PixelType v = imIter.Get();
	      
                for (long i4 = 0; i4 < n4; ++i4)
                  {
                    (*array4)(ix, iy, iz, i4) = v[i4];
                  }

                ++imIter;
              }
          }
      }

    typename VectorImageType::SpacingType sp = vctImage->GetSpacing();
    array4->setSpacingX(sp[0]);
    array4->setSpacingY(sp[1]);
    array4->setSpacingZ(sp[2]);


    return array4;
  }


  /*============================================================
    saveAsImage2  */
  template< typename T > 
  void saveAsImage2(typename cArray2D< T >::PointerType array2, const char *fileName)
  {
    typedef itk::Image< T, 2 > ImageType; // 2 means 2D  

    typename ImageType::Pointer image = cArray2ToItkImage< T >(array2);

    typedef itk::ImageFileWriter< ImageType > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(image);

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

    return;
  }


  /*============================================================
    saveAsImage3  */
  template< typename T > 
  void saveAsImage3(typename cArray3D< T >::PointerType array3, const char *fileName)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D  

    typename ImageType::Pointer image = cArray3ToItkImage< T >(array3);

    typedef itk::ImageFileWriter< ImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(image);

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

    return;
  }


  /*============================================================
    saveAsImage4  */
  template< typename T > 
  void saveAsImage4(typename cArray4D< T >::PointerType array4, const char *fileName)
  {
    typedef itk::VectorImage< T, 3 > VectorImageType;
    typename VectorImageType::Pointer image = VectorImageType::New();

    long vctLength = array4->getSize4();
    image->SetVectorLength(vctLength);
  
    typename VectorImageType::IndexType imStart;
    imStart[0] = 0;  // first index on X
    imStart[1] = 0;  // first index on Y
    imStart[2] = 0;  // first index on Y

    typename VectorImageType::SizeType  imSize;
    imSize[0] = array4->getSizeX();  // size along X
    imSize[1] = array4->getSizeY();  // size along Y
    imSize[2] = array4->getSizeZ();  // size along Y

    typename VectorImageType::RegionType imRegion;
    imRegion.SetSize( imSize );
    imRegion.SetIndex( imStart );

    typename VectorImageType::SpacingType sp;
    sp[0] = array4->getSpacingX();
    sp[1] = array4->getSpacingY();
    sp[2] = array4->getSpacingZ();

    image->SetSpacing(sp);
    image->SetRegions( imRegion );
    image->Allocate();


    typedef itk::ImageRegionIterator< VectorImageType > VectorImageIteratorType;
    VectorImageIteratorType imIter(image, image->GetLargestPossibleRegion() );

    typename VectorImageType::RegionType::SizeType dataSz = (image->GetLargestPossibleRegion()).GetSize();

    imIter.GoToBegin();
    for (long iz = 0; iz < dataSz[2]; ++iz)
      {
        for (long iy = 0; iy < dataSz[1]; ++iy)
          {
            for (long ix = 0; ix < dataSz[0]; ++ix)
              {
                /*
                  NOTE!!! the order: iz changes slowest while ix changes
                  fastest MATTERS! 

                  Since this is using together with iterator of
                  itkImage. There the iterator first goes thru x-dim,
                  then y-dim, then z-dim. So the fastest changing
                  dimension should be consistant.
                */

                typename VectorImageType::PixelType v(vctLength);
	      
                for (long i4 = 0; i4 < vctLength; ++i4)
                  {
                    v[i4] = (*array4)(ix, iy, iz, i4);
                  }

                imIter.Set(v);

                ++imIter;
              }
          }
      }


    typedef itk::ImageFileWriter< VectorImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(image);

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


    return;
  }

  /*--------------------------------------------------
    saveAsMat2
  */
  template< typename T > 
  void saveAsMat2(typename cArray2D< T >::PointerType array2, const char *fileName)
  {
    /*
      Port to matlab .mat file:

      Then in matlab, when load in the .mat file, there will be two
      variables: data and size. Do data = reshape(data, size'); then
      imagesc(data) to show the data.
    */

    // The writer doesn't support int etc. so it's safe to hard code
    // here to avoid the problem when T is some un-supported type of the
    // writer.
    long nx = array2->getSizeX();
    long ny = array2->getSizeY();
    long n = nx*ny;

    vnl_vector<double> array_vnl( n );
    for (long i = 0; i < n; ++i)
      {
        array_vnl[i] = (double)(*array2)[i];
      }

    vnl_vector< double > size_vnl(2); 
    size_vnl[0] = nx;
    size_vnl[1] = ny;
  
    vnl_matlab_filewrite matlabFileWriter(fileName);
    matlabFileWriter.write(array_vnl, "data");
    matlabFileWriter.write(size_vnl, "sz"); // not "size" to avoid confliction with Matlab's size fn.
  }


  /*--------------------------------------------------
    saveAsMat3
  */
  template< typename T > 
  void saveAsMat3(typename cArray3D< T >::PointerType array3, const char *fileName)
  {
    /*
      Port to matlab .mat file:

      Then in matlab, when load in the .mat file, there will be two
      variables: data and size. Do data = reshape(data, size'); then
      imagesc(data) to show the data.
    */

    // The writer doesn't support int etc. so it's safe to hard code
    // here to avoid the problem when T is some un-supported type of the
    // writer.
    long nx = array3->getSizeX();
    long ny = array3->getSizeY();
    long nz = array3->getSizeZ();
    long n = nx*ny*nz;

    vnl_vector< double > array_vnl( n );
    for (long i = 0; i < n; ++i)
      {
        array_vnl[i] = (double)(*array3)[i];
      }

    vnl_vector<double> size_vnl(3); 
    size_vnl[0] = nx;
    size_vnl[1] = ny;
    size_vnl[2] = nz;
  
    vnl_matlab_filewrite matlabFileWriter(fileName);
    matlabFileWriter.write(array_vnl, "data");
    matlabFileWriter.write(size_vnl, "sz"); // not "size" to avoid confliction with Matlab's size fn.
  }


  /*--------------------------------------------------
    saveAsMat4  */
  template< typename T > 
  void saveAsMat4(typename cArray4D< T >::PointerType array4, const char *fileName)
  {
    /*
      Port to matlab .mat file:

      Then in matlab, when load in the .mat file, there will be two
      variables: data and size. Do data = reshape(data, size'); then
      imagesc(data) to show the data.
    */

    // The writer doesn't support int etc. so it's safe to hard code
    // here to avoid the problem when T is some un-supported type of the
    // writer.
    long nx = array4->getSizeX();
    long ny = array4->getSizeY();
    long nz = array4->getSizeZ();
    long n4 = array4->getSize4();

    long n = nx*ny*nz*n4;

    vnl_vector< double > array_vnl( n );
    for (long i = 0; i < n; ++i)
      {
        array_vnl[i] = (double)(*array4)[i];
      }

    vnl_vector<double> size_vnl(4); 
    size_vnl[0] = nx;
    size_vnl[1] = ny;
    size_vnl[2] = nz;
    size_vnl[3] = n4;
  
    vnl_matlab_filewrite matlabFileWriter(fileName);
    matlabFileWriter.write(array_vnl, "data");
    matlabFileWriter.write(size_vnl, "sz"); // not "size" to avoid confliction with Matlab's size fn.
  }



  /*--------------------------------------------------
    rescale array
  */
  template< typename T > cArray2D< T >* rescaleArray(cArray2D< T >* array2, T outputMin, T outputMax)
  {
    if (outputMin >= outputMax)
      {
        std::cerr<<"outputMin >= outputMax, exit. \n";
        raise(SIGABRT);
      }


    long nx = array2->getSizeX();
    long ny = array2->getSizeY();
    long n = nx*ny;
    if (n <= 1)
      {
        std::cerr<<"size of array <= 1, exit. \n";
        raise(SIGABRT);
      }


    T inputMin = (*array2)[0];
    T inputMax = (*array2)[0];

    for (long i = 0; i < n; ++i)
      {
        T v = (*array2)[i];
        inputMin = inputMin>v?v:inputMin;
        inputMax = inputMax>v?inputMax:v;
      }


    T outputDiff = outputMax - outputMin;
    T inputDiff = inputMax - inputMin;

    if (inputDiff == 0)
      {
        std::cerr<<"input data has inputMax == inputMin, exit.\n";
        raise(SIGABRT);
      }

  
    cArray2D< T >* newArray = new cArray2D< T >(*array2);
    for (long i = 0; i < n; ++i)
      {
        T v = (*newArray)[i];
        T newV = outputMin + (v - inputMin)*outputDiff/(inputDiff + vnl_math::eps);
        (*newArray)[i] = newV;
      }

    return newArray;
  }

  template< typename T > cArray3D< T >* rescaleArray(cArray3D< T >* array3, T outputMin, T outputMax)
  {
    if (outputMin >= outputMax)
      {
        std::cerr<<"outputMin >= outputMax, exit. \n";
        raise(SIGABRT);
      }


    long nx = array3->getSizeX();
    long ny = array3->getSizeY();
    long nz = array3->getSizeZ();
    long n = nx*ny*nz;

    if (n <= 1)
      {
        std::cerr<<"size of array <= 1, exit. \n";
        raise(SIGABRT);
      }

    T inputMin = (*array3)[0];
    T inputMax = (*array3)[0];

    for (long i = 0; i < n; ++i)
      {
        T v = (*array3)[i];
        inputMin = inputMin>v?v:inputMin;
        inputMax = inputMax>v?inputMax:v;
      }


    T outputDiff = outputMax - outputMin;
    T inputDiff = inputMax - inputMin;

    if (inputDiff == 0)
      {
        std::cerr<<"input data has inputMax == inputMin, exit.\n";
        raise(SIGABRT);
      }

  
    cArray3D<T>* newArray = new cArray3D< T >(*array3);
    for (long i = 0; i < n; ++i)
      {
        T v = (*newArray)[i];
        T newV = outputMin + (v - inputMin)*outputDiff/(inputDiff + vnl_math::eps);
        (*newArray)[i] = newV;
      }

    return newArray;
  }

  template< typename T > cArray4D< T >* rescaleArray(cArray4D< T >* array4, T outputMin, T outputMax)
  {
    if (outputMin >= outputMax)
      {
        std::cerr<<"outputMin >= outputMax, exit. \n";
        raise(SIGABRT);
      }


    long nx = array4->getSizeX();
    long ny = array4->getSizeY();
    long nz = array4->getSizeZ();
    long n4 = array4->getSize4();
    long n = nx*ny*nz*n4;

    if (n <= 1)
      {
        std::cerr<<"size of array <= 1, exit. \n";
        raise(SIGABRT);
      }

    T inputMin = (*array4)[0];
    T inputMax = (*array4)[0];

    for (long i = 0; i < n; ++i)
      {
        T v = (*array4)[i];
        inputMin = inputMin>v?v:inputMin;
        inputMax = inputMax>v?inputMax:v;
      }


    T outputDiff = outputMax - outputMin;
    T inputDiff = inputMax - inputMin;

    if (inputDiff == 0)
      {
        std::cerr<<"input data has inputMax == inputMin, exit.\n";
        raise(SIGABRT);
      }

  
    cArray4D<T>* newArray = new cArray4D< T >(*array4);
    for (long i = 0; i < n; ++i)
      {
        T v = (*newArray)[i];
        T newV = outputMin + (v - inputMin)*outputDiff/(inputDiff + vnl_math::eps);
        (*newArray)[i] = newV;
      }

    return newArray;
  }


  /* ==================================================
     clamp arrays 2D */
  template< typename T > 
  typename cArray2D< T >::PointerType clampBelow2D( const typename cArray2D< T >::PointerType array, T clampValue, T setToValue = 0.0)
  {
    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long n = nx*ny;

    for (long i = 0; i < n; ++i)
      {
        T v = (*array)[i];

        (*array)[i] = v>=clampValue?v:setToValue;
      }
  }

  /* ==================================================
     clamp arrays 3D */
  template< typename T > 
  typename cArray3D< T >::PointerType clampBelow3D( const typename cArray3D< T >::PointerType array, T clampValue, T setToValue = 0.0)
  {
    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long nz = array->getSizeZ();
    long n = nx*ny*nz;

    for (long i = 0; i < n; ++i)
      {
        T v = (*array)[i];

        (*array)[i] = v>=clampValue?v:setToValue;
      }
  }


  /* ==================================================
     clamp arrays 4D */
  template< typename T > 
  typename cArray4D< T >::PointerType clampBelow4D( const typename cArray4D< T >::PointerType array, T clampValue, T setToValue = 0.0)
  {
    long nx = array->getSizeX();
    long ny = array->getSizeY();
    long nz = array->getSizeZ();
    long n4 = array->getSize4();
    long n = nx*ny*nz*n4;

    for (long i = 0; i < n; ++i)
      {
        T v = (*array)[i];

        (*array)[i] = v>=clampValue?v:setToValue;
      }
  }





  /* ============================================================
     itkImageToArray2  */
  template< typename T > 
  typename cArray2D< T >::PointerType itkImageToArray2(typename itk::Image< T, 2 >::Pointer image)
  {
    typedef itk::Image< T, 2 > ImageType; // 2 means 2D  
  
    typename ImageType::RegionType::SizeType dataSz = (image->GetLargestPossibleRegion()).GetSize();

    typename cArray2D<T>::PointerType array2( new cArray2D< T >(dataSz[0], dataSz[1]) );

    typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;
    ImageIteratorType imIter(image, image->GetLargestPossibleRegion() );

    imIter.GoToBegin();

    for (unsigned long iy = 0; iy < dataSz[1]; ++iy)
      {
        for (unsigned long ix = 0; ix < dataSz[0]; ++ix)
          {
            /*
              NOTE!!! the order: iz changes slowest while ix changes
              fastest MATTERS! 

              Since this is using together with iterator of
              itkImage. There the iterator first goes thru x-dim,
              then y-dim, then z-dim. So the fastest changing
              dimension should be consistant.
            */
            (*array2)(ix, iy) = imIter.Get();

            ++imIter;
          }
      }

    typename ImageType::SpacingType sp = image->GetSpacing();
    array2->setSpacingX(sp[0]);
    array2->setSpacingY(sp[1]);

    return array2;
  }

  /*============================================================
    itkImageToArray3 */
  template< typename T > 
  typename cArray3D< T >::PointerType itkImageToArray3(typename itk::Image< T, 3 >::Pointer image)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D  
    typename ImageType::RegionType::SizeType dataSz = (image->GetLargestPossibleRegion()).GetSize();

    typename cArray3D<T>::PointerType array3( new cArray3D< T >(dataSz[0], dataSz[1], dataSz[2]) );

    typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;
    ImageIteratorType imIter(image, image->GetLargestPossibleRegion() );

    imIter.GoToBegin();

    for (unsigned long iz = 0; iz < dataSz[2]; ++iz)
      {
        for (unsigned long iy = 0; iy < dataSz[1]; ++iy)
          {
            for (unsigned long ix = 0; ix < dataSz[0]; ++ix)
              {
                /*
                  NOTE!!! the order: iz changes slowest while ix changes
                  fastest MATTERS! 

                  Since this is using together with iterator of
                  itkImage. There the iterator first goes thru x-dim,
                  then y-dim, then z-dim. So the fastest changing
                  dimension should be consistant.
                */
                (*array3)(ix, iy, iz) = imIter.Get();

                ++imIter;
              }
          }
      }

    typename ImageType::SpacingType sp = image->GetSpacing();
    array3->setSpacingX(sp[0]);
    array3->setSpacingY(sp[1]);
    array3->setSpacingZ(sp[2]);


    return array3;
  }


  /*============================================================
    cArray2ToItkImage  */
  template< typename T > 
  typename itk::Image< T, 2 >::Pointer cArray2ToItkImage(typename cArray2D< T >::PointerType array2)
  {
    typedef itk::Image< T, 2 > ImageType; // 2 means 2D  

    typename ImageType::Pointer image = ImageType::New();
  
    typename ImageType::IndexType imStart;
    imStart[0] = 0;  // first index on X
    imStart[1] = 0;  // first index on Y

    typename ImageType::SizeType  imSize;
    imSize[0] = array2->getSizeX();  // size along X
    imSize[1] = array2->getSizeY();  // size along Y

    typename ImageType::RegionType imRegion;
    imRegion.SetSize( imSize );
    imRegion.SetIndex( imStart );


    typename ImageType::SpacingType sp;
    sp[0] = array2->getSpacingX();
    sp[1] = array2->getSpacingY();

    image->SetSpacing(sp);
    image->SetRegions( imRegion );
    image->Allocate();

    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType imIt(image, image->GetLargestPossibleRegion() );
    imIt.GoToBegin();

    for (long i = 0; !imIt.IsAtEnd(); ++imIt)
      {
        imIt.Set((*array2)[i]);
        ++i;
      }


    return image;
  }


  /*============================================================
    cArray3ToItkImage  */
  template< typename T > 
  typename itk::Image< T, 3 >::Pointer cArray3ToItkImage(typename cArray3D< T >::PointerType array3)
  {
    typedef itk::Image< T, 3 > ImageType; // 3 means 3D  

    typename ImageType::Pointer image = ImageType::New();
  
    typename ImageType::IndexType imStart;
    imStart[0] = 0;  // first index on X
    imStart[1] = 0;  // first index on Y
    imStart[2] = 0;  // first index on Y

    typename ImageType::SizeType  imSize;
    imSize[0] = array3->getSizeX();  // size along X
    imSize[1] = array3->getSizeY();  // size along Y
    imSize[2] = array3->getSizeZ();  // size along Y

    typename ImageType::RegionType imRegion;
    imRegion.SetSize( imSize );
    imRegion.SetIndex( imStart );

    typename ImageType::SpacingType sp;
    sp[0] = array3->getSpacingX();
    sp[1] = array3->getSpacingY();
    sp[2] = array3->getSpacingZ();

    image->SetSpacing(sp);
    image->SetRegions( imRegion );
    image->Allocate();

    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType imIt(image, image->GetLargestPossibleRegion() );
    imIt.GoToBegin();

    for (long i = 0; !imIt.IsAtEnd(); ++imIt)
      {
        imIt.Set((*array3)[i]);
        ++i;
      }


    return image;
  }


  /*============================================================
    crop array  */
  template< typename T > 
  typename cArray2D< T >::PointerType cropArray(typename cArray2D< T >::PointerType array2, \
                                                long xMin, long xMax, long yMin, long yMax)
  {
    /*
      xMin is the min index, in the original array, to be *included* in
      the cropped array. While xMax is the max, in the original array,
      to be *included* in the cropped array.

      So the limit of xMin is 0, indicating: no crop at head. And the
      limit of xMax is nx-1, indicating: no crop at end.

      Same for y-dim.
    */
    long nx = array2->getSizeX();
    long ny = array2->getSizeY();

    if (xMin < 0 || xMax > nx-1 || yMin < 0 || yMax > ny-1 )
      {
        std::cerr<<"cropping out of bounds, exit.\n";
        raise(SIGABRT);
      }

    if (xMin > xMax || yMin > yMax)
      {
        std::cerr<<"xMin > xMax || yMin > yMax, exit.\n";
        raise(SIGABRT);
      }

    long szX = xMax - xMin + 1;
    long szY = yMax - yMin + 1;

    cArray2D< T >* cropped = new cArray2D< T >(szX, szY);
  
    for (long ix = xMin; ix <= xMax; ++ix)
      {
        for (long iy = yMin; iy <= yMax; ++iy)
          {
            // Notice the "=" in the upper limit, since we want to
            // include the xMax and yMax
            long iix = ix - xMin;
            long iiy = iy - yMin;

            (*cropped)(iix, iiy) = (*array2)(ix, iy);
          }
      }


    return cropped;
  }


  
  /* ============================================================
     cropArray 3D */
  template< typename T > 
  typename cArray3D< T >::PointerType cropArray(typename cArray3D< T >::PointerType array3, \
                                                long xMin, long xMax, long yMin, long yMax, long zMin, long zMax)
  {
    /*
      xMin is the min index, in the original array, to be *included* in
      the cropped array. While xMax is the max, in the original array,
      to be *included* in the cropped array.

      So the limit of xMin is 0, indicating: no crop at head. And the
      limit of xMax is nx-1, indicating: no crop at end.

      Same for y-dim.
    */
    long nx = array3->getSizeX();
    long ny = array3->getSizeY();
    long nz = array3->getSizeZ();

    if (xMin < 0 || xMax > nx-1 || yMin < 0 || yMax > ny-1 || zMin < 0 || zMax > nz-1)
      {
        std::cerr<<"cropping out of bounds, exit.\n";
        raise(SIGABRT);
      }

    if (xMin > xMax || yMin > yMax || zMin > zMax)
      {
        std::cerr<<"xMin > xMax || yMin > yMax || zMin > zMax, exit.\n";
        raise(SIGABRT);
      }

    long szX = xMax - xMin + 1;
    long szY = yMax - yMin + 1;
    long szZ = zMax - zMin + 1;

    cArray3D< T >* cropped = new cArray3D< T >(szX, szY, szZ);
  
    for (long ix = xMin; ix <= xMax; ++ix)
      {
        for (long iy = yMin; iy <= yMax; ++iy)
          {
            for (long iz = zMin; iz <= zMax; ++iz)
              {
                // Notice the "=" in the upper limit, since we want to
                // include the xMax and yMax and zMax
                long iix = ix - xMin;
                long iiy = iy - yMin;
                long iiz = iz - zMin;

                (*cropped)(iix, iiy, iiz) = (*array3)(ix, iy, iz);
              }
          }
      }


    return cropped;
  }

  /* ============================================================
     cropArray 4D */
  template< typename T > 
  typename cArray4D< T >::PointerType cropArray(typename cArray4D< T >::PointerType array4, \
                                                long xMin, long xMax, long yMin, long yMax, long zMin, long zMax, \
                                                long fourMin, long fourMax)
  {
    /*
      xMin is the min index, in the original array, to be *included* in
      the cropped array. While xMax is the max, in the original array,
      to be *included* in the cropped array.

      So the limit of xMin is 0, indicating: no crop at head. And the
      limit of xMax is nx-1, indicating: no crop at end.

      Same for y-dim.
    */
    long nx = array4->getSizeX();
    long ny = array4->getSizeY();
    long nz = array4->getSizeZ();
    long n4 = array4->getSize4();


    if (xMin < 0 || xMax > nx-1 || yMin < 0 || yMax > ny-1 || zMin < 0 || zMax > nz-1 || fourMin < 0 || fourMax > n4-1)
      {
        std::cerr<<"cropping out of bounds, exit.\n";
        raise(SIGABRT);
      }

    if (xMin > xMax || yMin > yMax || zMin > zMax || fourMin > fourMax)
      {
        std::cerr<<"xMin > xMax || yMin > yMax || zMin > zMax || fourMin > fourMax, exit.\n";
        raise(SIGABRT);
      }

    long szX = xMax - xMin + 1;
    long szY = yMax - yMin + 1;
    long szZ = zMax - zMin + 1;
    long sz4 = fourMax - fourMin + 1;

    cArray4D< T >* cropped = new cArray4D< T >(szX, szY, szZ, sz4);
  
    for (long ix = xMin; ix <= xMax; ++ix)
      {
        for (long iy = yMin; iy <= yMax; ++iy)
          {
            for (long iz = zMin; iz <= zMax; ++iz)
              {
                for (long i4 = fourMin; i4 <= fourMax; ++i4)
                  {
                    // Notice the "=" in the upper limit, since we want to
                    // include the xMax, yMax, zMax and fourMax
                    long iix = ix - xMin;
                    long iiy = iy - yMin;
                    long iiz = iz - zMin;
                    long ii4 = i4 - fourMin;

                    (*cropped)(iix, iiy, iiz, ii4) = (*array4)(ix, iy, iz, i4);
                  }
              }
          }
      }


    return cropped;
  }


  /*============================================================
    const pad array 2D */
  template< typename T > 
  typename cArray2D< T >::PointerType constPadArray(typename cArray2D< T >::PointerType array2, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum)
  {
    long nx = array2->getSizeX();
    long ny = array2->getSizeY();

    long newNx = nx + xHeadNum + xTailNum;
    long newNy = ny + yHeadNum + yTailNum;

    typename cArray2D< T >::PointerType padded( new cArray2D< T >(newNx, newNy, padValue) );

    for (long ix = 0; ix < nx; ++ix)
      {
        for (long iy = 0; iy < ny; ++iy)
          {
            long newIx = ix + xHeadNum;
            long newIy = iy + yHeadNum;

            (*padded)(newIx, newIy) = (*array2)(ix, iy);
          }
      }

    return padded;
  }

  /*============================================================
    const pad array 3D  */
  template< typename T > 
  typename cArray3D< T >::PointerType constPadArray(typename cArray3D< T >::PointerType array3, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum, \
                                                    long zHeadNum, long zTailNum)
  {
    long nx = array3->getSizeX();
    long ny = array3->getSizeY();
    long nz = array3->getSizeZ();

    long newNx = nx + xHeadNum + xTailNum;
    long newNy = ny + yHeadNum + yTailNum;
    long newNz = nz + zHeadNum + zTailNum;

    typename cArray3D< T >::PointerType padded( new cArray3D< T >(newNx, newNy, newNz, padValue) );

    for (long ix = 0; ix < nx; ++ix)
      {
        for (long iy = 0; iy < ny; ++iy)
          {
            for (long iz = 0; iz < nz; ++iz)
              {
                long newIx = ix + xHeadNum;
                long newIy = iy + yHeadNum;
                long newIz = iz + zHeadNum;

                (*padded)(newIx, newIy, newIz) = (*array3)(ix, iy, iz);
              }
          }
      }

    return padded;
  }

  /*============================================================
    const pad array 4D  */
  template< typename T > 
  typename cArray4D< T >::PointerType constPadArray(typename cArray4D< T >::PointerType array4, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum, \
                                                    long zHeadNum, long zTailNum, \
                                                    long fourHeadNum, long fourTailNum)
  {
    long nx = array4->getSizeX();
    long ny = array4->getSizeY();
    long nz = array4->getSizeZ();
    long n4 = array4->getSize4();

    long newNx = nx + xHeadNum + xTailNum;
    long newNy = ny + yHeadNum + yTailNum;
    long newNz = nz + zHeadNum + zTailNum;
    long newN4 = n4 + fourHeadNum + fourTailNum;

    typename cArray4D< T >::PointerType padded( new cArray4D< T >(newNx, newNy, newNz, newN4, padValue) );

    for (long ix = 0; ix < nx; ++ix)
      {
        for (long iy = 0; iy < ny; ++iy)
          {
            for (long iz = 0; iz < nz; ++iz)
              {
                for (long i4 = 0; i4 < n4; ++i4)
                  {
                    long newIx = ix + xHeadNum;
                    long newIy = iy + yHeadNum;
                    long newIz = iz + zHeadNum;
                    long newI4 = i4 + fourHeadNum;

                    (*padded)(newIx, newIy, newIz, newI4) = (*array4)(ix, iy, iz, i4);
                  }
              }
          }
      }

    return padded;
  }

  /* ============================================================
     histogram */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  hist1D(const typename cArray2D< DataType >::PointerType data, DataType min, DataType max)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY());
    return hist1D(data->getDataPointer(), n, min, max);
  }

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  hist1D(const typename cArray2D< DataType >::PointerType data)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY());
    return hist1D(data->getDataPointer(), n);
  }


  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  hist1D(const typename cArray3D< DataType >::PointerType data, DataType min, DataType max)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ());
    return hist1D(data->getDataPointer(), n, min, max);
  }

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  hist1D(const typename cArray3D< DataType >::PointerType data)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ());
    return hist1D(data->getDataPointer(), n);
  }

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  hist1D(const typename cArray4D< DataType >::PointerType data, DataType min, DataType max)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ())*(data->getSize4());
    return hist1D(data->getDataPointer(), n, min, max);
  }


  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  hist1D(const typename cArray4D< DataType >::PointerType data)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ())*(data->getSize4());
    return hist1D(data->getDataPointer(), n);
  }


  /*================================================================================
     Gaussian kernel density estimation of 1D data  */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray2D< DataType >::PointerType data, double stdDev, DataType min, DataType max)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY());
    return kde1D(data->getDataPointer(), n, stdDev, min, max);
  }

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray2D< DataType >::PointerType data, double stdDev)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY());
    return kde1D(data->getDataPointer(), n, stdDev);
  }


  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray3D< DataType >::PointerType data, double stdDev, DataType min, DataType max)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ());
    return kde1D(data->getDataPointer(), n, stdDev, min, max);
  }

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray3D< DataType >::PointerType data, double stdDev)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ());
    return kde1D(data->getDataPointer(), n, stdDev);
  }

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray4D< DataType >::PointerType data, double stdDev, DataType min, DataType max)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ())*(data->getSize4());
    return kde1D(data->getDataPointer(), n, stdDev, min, max);
  }


  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray4D< DataType >::PointerType data, double stdDev)
  {
    unsigned long n = (data->getSizeX())*(data->getSizeY())*(data->getSizeZ())*(data->getSize4());
    return kde1D(data->getDataPointer(), n, stdDev);
  }



} // douher

////////////////////////////////////////////////////////////
#endif


