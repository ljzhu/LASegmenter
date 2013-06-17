#ifndef imageProcessing_h_
#define imageProcessing_h_


#include "cArrayOp.h"


#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_matrix.h"

#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_vector.h"


namespace douher
{
  /*============================================================
    gradient2 */
  template< typename T > 
  void gradientImage2(typename itk::Image< T, 2 >::Pointer image, \
                      typename itk::Image< T, 2 >::Pointer &gradX, \
                      typename itk::Image< T, 2 >::Pointer &gradY, \
                      unsigned char zeroBoundary = 0);
  // When returning newly inner generate function by pointer in input
  // parameters, remember to use &


  /*==================================================
    gradient3 */
  template< typename T > 
  void gradientImage3(typename itk::Image< T, 3 >::Pointer image, \
                      typename itk::Image< T, 3 >::Pointer &gradX, \
                      typename itk::Image< T, 3 >::Pointer &gradY, \
                      typename itk::Image< T, 3 >::Pointer &gradZ, \
                      unsigned char zeroBoundary = 0);


  /* ============================================================
     similarity Transformation  */
  template< typename T > typename itk::Image< T, 2 >::Pointer similarityTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                                                   vnl_vector< double > parameters, \
                                                                                   T fillValue = 0, \
                                                                                   unsigned char interpType = 1);

  template< typename T > typename itk::Image< T, 2 >::Pointer similarityTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                                                   double translateX, double translateY, \
                                                                                   double rotate, \
                                                                                   double scale, \
                                                                                   T fillValue = 0, \
                                                                                   unsigned char interpType = 1);


  template< typename T > typename itk::Image< T, 3 >::Pointer similarityTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                                                   vnl_vector< double > parameters, \
                                                                                   T fillValue = 0, \
                                                                                   unsigned char interpType = 1);

  template< typename T > typename itk::Image< T, 3 >::Pointer similarityTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                                                   double translateX, \
                                                                                   double translateY, \
                                                                                   double translateZ, \
                                                                                   double rotateAroundX, \
                                                                                   double rotateAroundY, \
                                                                                   double rotateAroundZ, \
                                                                                   double scale, \
                                                                                   T fillValue = 0, \
                                                                                   unsigned char interpType = 1);
  /* threshold Images*/
  template< typename T, unsigned int dim > 
  typename itk::Image< unsigned char, dim >::Pointer thresholdImage(typename itk::Image< T, dim >::Pointer image, T thld);

  template< typename T > 
  typename itk::Image< unsigned char, 2 >::Pointer thresholdImage2(typename itk::Image< T, 2 >::Pointer image, T thld);

  template< typename T > 
  typename itk::Image< unsigned char, 3 >::Pointer thresholdImage3(typename itk::Image< T, 3 >::Pointer image, T thld);




  /*************************************************************
   affine transform */
  template< typename T > 
  typename itk::Image< T, 2 >::Pointer affineTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                        vnl_matrix_fixed< double, 2, 2 > A, \
                                                        vnl_vector_fixed< double, 2 > translation, \
                                                        T fillValue, \
                                                        unsigned char interpType = 1, \
                                                        long newNx = -1, long newNy = -1);

  template< typename T > 
  typename itk::Image< T, 2 >::Pointer affineTransform2(typename itk::Image< T, 2 >::Pointer image, \
                                                        vnl_matrix< double > A, \
                                                        vnl_vector< double > translation, \
                                                        T fillValue, \
                                                        unsigned char interpType = 1, \
                                                        long newNx = -1, long newNy = -1);

  template< typename T > 
  typename itk::Image< T, 3 >::Pointer affineTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                        vnl_matrix_fixed< double, 3, 3 > A, \
                                                        vnl_vector_fixed< double, 3 > translation, \
                                                        T fillValue, \
                                                        unsigned char interpType = 1,\
                                                        long newNx = -1, long newNy = -1, long newNz = -1);

  template< typename T > 
  typename itk::Image< T, 3 >::Pointer affineTransform3(typename itk::Image< T, 3 >::Pointer image, \
                                                        vnl_matrix< double > A, \
                                                        vnl_vector< double > translation, \
                                                        T fillValue, \
                                                        unsigned char interpType = 1, \
                                                        long newNx = -1, long newNy = -1, long newNz = -1);


  /*============================================================
    Intensity linear transform */
  template< typename Tin, typename Tout, unsigned int dim >
  typename itk::Image< Tout, dim >::Pointer intensityLinearTransform( typename itk::Image< Tin, dim >::Pointer img, \
                                                                      double a = 1, double b = 0);

  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 2 >::Pointer intensityLinearTransform2D( typename itk::Image< Tin, 2 >::Pointer img, \
                                                                      double a = 1, double b = 0);

  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 3 >::Pointer intensityLinearTransform3D( typename itk::Image< Tin, 3 >::Pointer img, \
                                                                      double a = 1, double b = 0);



  /*============================================================
    SDF */
  template< typename Tin, typename Tout, unsigned int dim >
  typename itk::Image< Tout, dim >::Pointer sdfOfBin( typename itk::Image< Tin, dim >::Pointer binImg);

  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 2 >::Pointer sdfOfBin2D( typename itk::Image< Tin, 2 >::Pointer binImg);

  template< typename Tin, typename Tout >
  typename itk::Image< Tout, 3 >::Pointer sdfOfBin3D( typename itk::Image< Tin, 3 >::Pointer binImg);


  /*============================================================
    tanh */
  template< typename Tin, unsigned int dim >
  typename itk::Image< double, dim >::Pointer tanhOfImage( typename itk::Image< Tin, dim >::Pointer img);

  template< typename Tin >
  typename itk::Image< double, 2 >::Pointer tanhOfImage2D( typename itk::Image< Tin, 2 >::Pointer img);

  template< typename Tin >
  typename itk::Image< double, 3 >::Pointer tanhOfImage3D( typename itk::Image< Tin, 3 >::Pointer img);



  
  /* off grid image value */
  template< typename T, unsigned int dim > T offGridImageValue(typename itk::Image< T, dim >::Pointer image, \
                                                               itk::ContinuousIndex<double, dim> offGridIndex, \
                                                               unsigned char interpType = 1);

  template< typename T > T offGridImageValue(typename itk::Image< T, 2 >::Pointer image, \
                                             double x, double y, \
                                             unsigned char interpType = 1);

  template< typename T > T offGridImageValue(typename itk::Image< T, 3 >::Pointer image, \
                                             double x, double y, double z, \
                                             unsigned char interpType = 1);


} // namespace



#include "imageProcessing.txx"

#endif
