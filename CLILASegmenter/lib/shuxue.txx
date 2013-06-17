#ifndef _shuxue_txx_
#define _shuxue_txx_

#include "shuxue.h"

#include <iostream>
#include <cmath>

#include <vnl/vnl_matlab_filewrite.h>

namespace douher
{
  template< typename T > 
  inline T heaviside(T t, double epsilon)
  {
    return 0.5 + atan(t/epsilon)/(vnl_math::pi);
  }

  template< typename T > 
  inline T delta(T t, double epsilon)
  {
    return epsilon*epsilon/(vnl_math::pi)/(epsilon*epsilon + t*t);
  }

  template< typename T > 
  vnl_vector< T > *dwt1_cdf97(const vnl_vector< T > &f)
  {
    long n = f.size();
 
    if ( !(n&&!(n&(n-1))) ) // test if n is power of 2, see workLog for detail
      {
        std::cerr<<"size has to be power of 2, exit.\n";
        exit(-1);
      }

    vnl_vector< T > *x = new vnl_vector< T >(f);

    // Predict 1
    double a = 0;
    long i = 0;

    a = -1.586134342;
    for(i=1; i<n-2; i+=2) 
      {
        (*x)[i] += a*((*x)[i-1] + (*x)[i+1]);
      } 
    (*x)[n-1] += 2*a*(*x)[n-2];

    // Update 1
    a = -0.05298011854;
    for(i=2; i<n; i+=2) 
      {
        (*x)[i] += a*((*x)[i-1] + (*x)[i+1]);
      }
    (*x)[0] += 2*a*(*x)[1];

    // Predict 2
    a = 0.8829110762;
    for(i=1; i<n-2; i+=2) 
      {
        (*x)[i] += a*((*x)[i-1] + (*x)[i+1]);
      }
    (*x)[n-1] += 2*a*(*x)[n-2];

    // Update 2
    a = 0.4435068522;
    for(i=2; i<n; i+=2) 
      {
        (*x)[i] += a*((*x)[i-1] + (*x)[i+1]);
      }
    (*x)[0] += 2*a*(*x)[1];

    // Scale
    a = 1/1.149604398;
    for(i=0; i<n; i++) 
      {
        if(i%2) 
          {
            (*x)[i] *= a;
          }
        else 
          {
            (*x)[i] /= a;
          }
      }

    // Pack
    T *tempbank = new T[n];
  
    for(i=0;i<n;i++) 
      {
        if (i%2==0) 
          {
            tempbank[i/2]=(*x)[i];
          }
        else 
          {
            tempbank[n/2+i/2]=(*x)[i];
          }
      }
    for(i=0;i<n;i++) 
      {
        (*x)[i]=tempbank[i];
      }

    delete tempbank;

    return x;


  }

  /*==================================================*/
  template< typename T > 
  vnl_vector< T > *idwt1_cdf97(const vnl_vector< T > &coef)
  {
    long n = coef.size();
    
    if (!(n&&!(n&(n-1)))) // test if n is power of 2, see workLog for detail
      {
        std::cerr<<"size has to be power of 2. exit.\n";
        exit(-1);
      }

    vnl_vector< T > *f = new vnl_vector< T >(coef);


    // Unpack
    T *tempbank = new T[n];

    double a = 0;
    long i = 0;

    for (i=0; i<n/2; i++) 
      {
        tempbank[i*2]=(*f)[i];
        tempbank[i*2+1]=(*f)[i+n/2];
      }

    for (i=0; i<n; i++)
      {
        (*f)[i]=tempbank[i];
      }

    delete tempbank;


    // Undo scale
    a=1.149604398;
    for (i=0; i<n; i++) 
      {
        if (i%2) 
          {
            (*f)[i] *= a;
          }
        else
          {
            (*f)[i] /= a;
          }
      }

    // Undo update 2
    a=-0.4435068522;
    for (i=2; i<n; i+=2) 
      {
        (*f)[i] += a*((*f)[i-1] + (*f)[i+1]);
      }
    (*f)[0] += 2*a*(*f)[1];

    // Undo predict 2
    a=-0.8829110762;
    for (i=1; i<n-2; i+=2) 
      {
        (*f)[i] += a*((*f)[i-1] + (*f)[i+1]);
      }
    (*f)[n-1] += 2*a*(*f)[n-2];

    // Undo update 1
    a=0.05298011854;
    for (i=2; i<n; i+=2) 
      {
        (*f)[i] += a*((*f)[i-1] + (*f)[i+1]);
      }
    (*f)[0] += 2*a*(*f)[1];

    // Undo predict 1
    a=1.586134342;
    for (i=1; i<n-2; i+=2) 
      {
        (*f)[i] += a*((*f)[i-1] + (*f)[i+1]);
      } 
    (*f)[n-1] += 2*a*(*f)[n-2];

    return f;
  }
                
  template< typename T > 
  vnl_vector< T > *pad_zero(const vnl_vector< T > &f, long newSize)
  {
    long oldSize = f.size();
    if (newSize < oldSize)
      {
        std::cerr<<"newSize should be >= oldSize. exit. \n";
        exit(-1);
      }
    
    vnl_vector< T > *newVct = new vnl_vector< T >(newSize, 0.0);

    for (long i = 0; i < oldSize; ++i)
      {
        (*newVct)[i] = f[i];
      }

    return newVct;
  }

  template< typename T > 
  vnl_vector< T > *pad_const(const vnl_vector< T > &f, long newSize)
  {
    long oldSize = f.size();
    if (newSize < oldSize)
      {
        std::cerr<<"newSize should be >= oldSize. exit. \n";
        exit(-1);
      }

    T lastVal = f[oldSize-1];    
    vnl_vector< T > *newVct = new vnl_vector< T >(newSize, lastVal);

    for (long i = 0; i < oldSize; ++i)
      {
        (*newVct)[i] = f[i];
      }

    return newVct;
  }


  template< typename T > 
  vnl_vector< T > *padToPowerOf2_zero(const vnl_vector< T > &f)
  {
    long n = f.size();
    if ( n&&!(n&(n-1)) ) // if n is already power of 2, just return a copy of it
      {
        vnl_vector< T > *padded = new vnl_vector< T >(f);
        return padded;
      }

    double b = log( (double)n )/log(2.0); // log base 2
    double intPartOfb = 0;
    modf(b, &intPartOfb);

    long newSize = (long)pow(2, intPartOfb+1);
    return pad_zero(f, newSize);
  }

  template< typename T > 
  vnl_vector< T > *padToPowerOf2_const(const vnl_vector< T > &f)
  {
    long n = f.size();
    if ( n&&!(n&(n-1)) ) // if n is already power of 2, just return a copy of it
      {
        vnl_vector< T > *padded = new vnl_vector< T >(f);
        return padded;
      }

    double b = log( (double)n )/log(2.0); // log base 2
    double intPartOfb = 0;
    modf(b, &intPartOfb);

    long newSize = (long)pow(2, intPartOfb+1);
    return pad_const(f, newSize);
  }


  /*----------------------------------------
    linspace   */
  //  VnlDoubleVectorType *linspace(double firstVal, double lastVal, long numOfVal)
  template< typename T > 
  vnl_vector< double > *linspace(double firstVal, double lastVal, long numOfVal)
  {
    // We exclude numOfVal = 1 case.
    if (numOfVal <= 1)
      {
        std::cerr<<"numOfVal <= 1, exit. \n";
        exit(-1);
      }

    vnl_vector< double > *f = new vnl_vector< double >(numOfVal);
    f->fill(firstVal);

    double stepSize = (lastVal - firstVal)/(double)(numOfVal - 1);
    for (long ii = 1; ii < numOfVal; ++ii)
      {
        (*f)[ii] = (*f)[ii - 1] + stepSize;
      }
    

    return f;
  }

  /* histogram of 1D data*/
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const DataType* data, unsigned long numData, DataType min, DataType max)
  {
    if (min >= max)
      {
	std::cerr<<"min >= max!!!??? abort. \n";
	raise(SIGABRT);
      }

    unsigned long lb = static_cast<unsigned long>(floor(min));
    unsigned long ub = static_cast<unsigned long>(ceil(max));
    unsigned long histoRange = ub - lb + 1;

    typedef vnl_vector< double > HistoFunctionType;
    typedef boost::shared_ptr< HistoFunctionType > HistoFunctionPointerType;

    double initH = 0.0;
    HistoFunctionPointerType h(new HistoFunctionType(histoRange, initH) );

    for (unsigned long i = 0; i < numData; ++i)
      {
	unsigned long pos = static_cast<unsigned long>( round(data[i]) );

	if (lb <= pos && pos <= ub)
	  {
	    h->put( pos-lb, h->get(pos-lb) + 1 );
	  }
      }

    return h;
  }

  /* histogram of 1D data, using its min and max as the bounds */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const DataType* data, unsigned long numData)
  {
    DataType minData = data[0];
    DataType maxData = minData;

    for (unsigned long i = 0; i < numData; ++i)
      {
	minData = minData<=data[i]?minData:data[i];
	maxData = maxData>=data[i]?maxData:data[i];
      }

    return hist1D(data, numData, minData, maxData);
  }

  /* histogram of vnl vector */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const boost::shared_ptr< vnl_vector< DataType > > data, DataType min, DataType max)
  {
    unsigned long n = data->size();

    return hist1D(data->data_block(), n, min, max);	
  }


  /* histogram of vnl vector, using its min and max as the bounds */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > hist1D(const boost::shared_ptr< vnl_vector< DataType > > data)
  {
    return hist1D(data, data->min_value(), data->max_value());
  }


  /* Gaussian kernel density estimation of 1D data */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const DataType* data, unsigned long numData, double stdDev, DataType min, DataType max)
  {
    if (min >= max)
      {
	std::cerr<<"min >= max!!!??? abort. \n";
	raise(SIGABRT);
      }

    if (stdDev <= 1e-6)
      {
	std::cerr<<"stdDev too small <= 1e-6. abort.\n";
	raise(SIGABRT);
      }

    /*----------------------------------------------------------------------
      Compute the neighbor size from the standard deviation. Assume a
      value can affect the histgram 3*stdDev away. Althogh the pdf may
      extend to region outside [min, max], we ignore those. Just to
      make sure that the size of the region of the returned vector is
      easier to handle. */
    unsigned long lb = static_cast<unsigned long>(floor(min));
    unsigned long ub = static_cast<unsigned long>(ceil(max));
    unsigned long histoRange = ub - lb + 1;

    unsigned int nbhdSize = static_cast< unsigned int >(3*stdDev);

    typedef vnl_vector< double > HistoFunctionType;
    typedef boost::shared_ptr< HistoFunctionType > HistoFunctionPointerType;

    double initH = 0.0;
    HistoFunctionPointerType h(new HistoFunctionType(histoRange, initH) );


    double coef = 1.0/(stdDev*sqrt(2.0*vnl_math::pi));
    double var2rec = -1.0/(2.0*stdDev*stdDev);

    for (unsigned long i = 0; i < numData; ++i)
      {
	DataType d = data[i];
	unsigned long floorD = static_cast<unsigned long>(floor(d));
	unsigned long ceilD = static_cast<unsigned long>(ceil(d));

	for (unsigned long ii = floorD - nbhdSize; ii <= ceilD + nbhdSize; ++ii)
	  {
	    if (lb <= ii && ii <= ub)
	      {
		double v = exp((d - ii)*(d - ii)*var2rec)*coef;
		h->put( ii, h->get(ii) + v );
	      }
	  }
      }

    (*h) /= numData;

    return h;    
  }

  /* Gaussian kernel density estimation of 1D data, using its min and max as the bounds */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const DataType* data, unsigned long numData, double stdDev)
  {
    DataType minData = data[0];
    DataType maxData = minData;

    for (unsigned long i = 0; i < numData; ++i)
      {
	minData = minData<=data[i]?minData:data[i];
	maxData = maxData>=data[i]?maxData:data[i];
      }

    return hist1D(data, numData, stdDev, minData, maxData);
  }


  /* Gaussian kernel density estimation of 1D data, for vnl_vector */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const boost::shared_ptr< vnl_vector< DataType > > data, double stdDev, DataType min, DataType max)
  {
    unsigned long n = data->size();

    return kde1D(data->data_block(), n, stdDev, min, max);	
  }


  /* Gaussian kernel density estimation of 1D data, using its min and max as the bounds */
  template< typename DataType > 
  boost::shared_ptr< vnl_vector< double > > kde1D(const boost::shared_ptr< vnl_vector< DataType > > data, double stdDev)
  {
    return kde1D(data, stdDev, data->min_value(), data->max_value());
  }


} // end of namespace

#endif
