/* ==================================================
   Operations of cArray classes.

   I'm not sure how much including those operations into cArray
   objects would harm their speed, so I just separate them here.

   20081031 Yi Gao

   ================================================== */

#ifndef cArrayOp_h_
#define cArrayOp_h_

#include "itkImage.h"

#include <vnl/vnl_vector.h>


#include "basicHeaders.h"


#include "cArray2D.h"
#include "cArray3D.h"
#include "cArray4D.h"


namespace douher
{
  /* IO with images   */
  template< typename T > 
  typename cArray2D< T >::PointerType readImageToArray2(const char *fileName);

  template< typename T > 
  typename cArray3D< T >::PointerType readImageToArray3(const char *fileName);

  template< typename T > 
  typename cArray4D< T >::PointerType readImageToArray4(const char *fileName);

  template< typename T > 
  void saveAsImage2(typename cArray2D< T >::PointerType array2, const char *fileName);

  template< typename T > 
  void saveAsImage3(typename cArray3D< T >::PointerType array3, const char *fileName);

  template< typename T > 
  void saveAsImage4(typename cArray4D< T >::PointerType array4, const char *fileName);


  /* IO: image file <===> itk::Image  */
  template< typename T > typename itk::Image< T, 2 >::Pointer readImage2(const char *fileName);
  template< typename T > typename itk::Image< T, 3 >::Pointer readImage3(const char *fileName);
  template< typename T, unsigned int dim > typename itk::Image< T, dim >::Pointer readImage(const char *fileName);

  template< typename T > void writeImage2(typename itk::Image< T, 2 >::Pointer img, const char *fileName);
  template< typename T > void writeImage3(typename itk::Image< T, 3 >::Pointer img, const char *fileName);
  template< typename T, unsigned int dim > void writeImage(typename itk::Image< T, dim >::Pointer img, const char *fileName);


  /* IO: vnl_vector <===> itk::Image  */
  template< typename T, unsigned int dim > 
  vnl_vector< T > itkImageToVnlVector(typename itk::Image< T, dim >::Pointer img);

  template< typename T, unsigned int dim > 
  vnl_vector< T > itkImageToVnlVector(typename itk::Image< T, dim >::ConstPointer img);


  template< typename T, unsigned int dim > 
  typename itk::Image< T, dim >::Pointer vnlVectorToItkImage(const vnl_vector< T >& vct, long imSize[]);


  /*  IO: cArray <===> itk::Image  */
  template< typename T > 
  typename cArray2D< T >::PointerType itkImageToArray2(typename itk::Image< T, 2 >::Pointer image);

  template< typename T > 
  typename cArray3D< T >::PointerType itkImageToArray3(typename itk::Image< T, 3 >::Pointer image);

  template< typename T > 
  typename itk::Image< T, 2 >::Pointer cArray2ToItkImage(typename cArray2D< T >::PointerType array2);

  template< typename T > 
  typename itk::Image< T, 3 >::Pointer cArray3ToItkImage(typename cArray3D< T >::PointerType array3);


  /* IO with Matlab mat files  */
  template< typename T > 
  void saveAsMat2(typename cArray2D< T >::PointerType array2, const char *fileName);

  template< typename T > 
  void saveAsMat3(typename cArray3D< T >::PointerType array3, const char *fileName);

  template< typename T > 
  void saveAsMat4(typename cArray4D< T >::PointerType array4, const char *fileName);


  /* rescale array   */
  template< typename T > 
  typename cArray2D< T >::PointerType rescaleArray(typename cArray2D< T >::PointerType array2, T outputMin, T outputMax);

  template< typename T > 
  typename cArray3D< T >::PointerType rescaleArray(typename cArray3D< T >::PointerType array3, T outputMin, T outputMax);

  template< typename T > 
  typename cArray4D< T >::PointerType rescaleArray(typename cArray4D< T >::PointerType array4, T outputMin, T outputMax);


  /* crop array  */
  template< typename T > 
  typename cArray2D< T >::PointerType cropArray(typename cArray2D< T >::PointerType array2, \
                                                long xMin, long xMax, long yMin, long yMax);

  template< typename T > 
  typename cArray3D< T >::PointerType cropArray(typename cArray3D< T >::PointerType array3, \
                                                long xMin, long xMax, long yMin, long yMax, long zMin, long zMax);
  
  template< typename T > 
  typename cArray4D< T >::PointerType cropArray(typename cArray4D< T >::PointerType array4, \
                                                long xMin, long xMax, long yMin, long yMax, long zMin, long zMax, \
                                                long fourMin, long fourMax);


  /* const pad array  */
  template< typename T > 
  typename cArray2D< T >::PointerType constPadArray(typename cArray2D< T >::PointerType array2, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum);

  template< typename T > 
  typename cArray3D< T >::PointerType constPadArray(typename cArray3D< T >::PointerType array3, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum, \
                                                    long zHeadNum, long zTailNum);

  template< typename T > 
  typename cArray4D< T >::PointerType constPadArray(typename cArray4D< T >::PointerType array4, T padValue, \
                                                    long xHeadNum, long xTailNum, \
                                                    long yHeadNum, long yTailNum, \
                                                    long zHeadNum, long zTailNum, \
                                                    long fourHeadNum, long fourTailNum);


  /* ============================================================
     clamp arrays:
     All array value below clampValue will be set to setToValue */
  template< typename T > 
  typename cArray2D< T >::PointerType clampBelow2D( const typename cArray2D< T >::PointerType array, T clampValue, T setToValue = 0.0);

  template< typename T > 
  typename cArray3D< T >::PointerType clampBelow3D( const typename cArray3D< T >::PointerType array, T clampValue, T setToValue = 0.0);

  template< typename T > 
  typename cArray4D< T >::PointerType clampBelow4D( const typename cArray4D< T >::PointerType array, T clampValue, T setToValue = 0.0);

  /* ============================================================
     histogram */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const typename cArray2D< DataType >::PointerType data, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const typename cArray2D< DataType >::PointerType data);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const typename cArray3D< DataType >::PointerType data, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const typename cArray3D< DataType >::PointerType data);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const typename cArray4D< DataType >::PointerType data, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const typename cArray4D< DataType >::PointerType data);


  /* Gaussian kernel density estimation of scalar data */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray2D< DataType >::PointerType data, double stdDev, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray2D< DataType >::PointerType data, double stdDev);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray3D< DataType >::PointerType data, double stdDev, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray3D< DataType >::PointerType data, double stdDev);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray4D< DataType >::PointerType data, double stdDev, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > 
  kde1D(const typename cArray4D< DataType >::PointerType data, double stdDev);
}

#include "cArrayOp.txx"

#endif
/* definition in .txx*/
