/* ==================================================
   Operations of cArray classes.

   I'm not sure how much including those operations into cArray
   objects would harm their speed, so I just separate them here.

   20081031 Yi Gao

   ================================================== */

#ifndef shuxue_H_
#define shuxue_H_

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>

#include "basicHeaders.h"


namespace douher
{
  /* padding and crop */
  template< typename T > vnl_vector< T > *pad_zero(const vnl_vector< T > &f, long newSize);
  template< typename T > vnl_vector< T > *pad_const(const vnl_vector< T > &f, long newSize);
  //  template< typename T > vnl_vector< T > *pad_period(const vnl_vector< T > &f, long padSize);


  template< typename T > vnl_vector< T > *padToPowerOf2_zero(const vnl_vector< T > &f);
  template< typename T > vnl_vector< T > *padToPowerOf2_const(const vnl_vector< T > &f);

  /*  some simple special functions */
  template< typename T > inline T heaviside(T t, double epsilon = 1.0);
  template< typename T > inline T delta(T t, double epsilon = 1.0);



  /* DWT with CDF97 wavelet  */
  template< typename T > vnl_vector< T > *dwt1_cdf97(const vnl_vector< T > &f);
  template< typename T > vnl_vector< T > *idwt1_cdf97(const vnl_vector< T > &f);

  /* linspace as in Matlab  */
  template< typename T > vnl_vector< double > *linspace(double firstVal, double lastVal, long numOfVal);

  /* histogram of 1D data*/
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const DataType* data, unsigned long numData, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const DataType* data, unsigned long numData);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const boost::shared_ptr< vnl_vector< DataType > > data, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const boost::shared_ptr< vnl_vector< DataType > > data);


  /* Gaussian kernel density estimation of 1D data */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const DataType* data, unsigned long numData, double stdDev, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const DataType* data, unsigned long numData, double stdDev);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const boost::shared_ptr< vnl_vector< DataType > > data, double stdDev, DataType min, DataType max);

  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const boost::shared_ptr< vnl_vector< DataType > > data, double stdDev);


} // end namespace 

#include "shuxue.txx"

#endif

