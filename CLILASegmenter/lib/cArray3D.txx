#ifndef cArray3D_txx_
#define cArray3D_txx_

#include "cArray3D.h"

#include <iostream>

namespace douher
{
  template<typename T> 
  cArray3D<T>::cArray3D()
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;

    m_data = 0;

    m_nx = 0;
    m_ny = 0;
    m_nz = 0;
  }


  template<typename T> 
  cArray3D<T>::~cArray3D()
  {
    clear();
  }


  template<typename T> 
  cArray3D<T>::cArray3D(long nx, long ny, long nz)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;

    m_data = new T[nx*ny*nz];
  }

  template<typename T>
  cArray3D<T>::cArray3D(long nx, long ny, long nz, T initVal)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;

    m_data = new T[nx*ny*nz];

    fill(initVal);
  }


  template<typename T>
  cArray3D<T>::cArray3D(const cArray3D &a)
  {
    m_spacingX = a.m_spacingX;
    m_spacingY = a.m_spacingY;
    m_spacingZ = a.m_spacingZ;

    m_nx = a.m_nx;
    m_ny = a.m_ny;
    m_nz = a.m_nz;

    long n = m_nx*m_ny*m_nz;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = a.m_data[i];
      }
  }



  template<typename T>
  cArray3D<T>::cArray3D(const PointerType p)
  {
    m_spacingX = p->m_spacingX;
    m_spacingY = p->m_spacingY;
    m_spacingZ = p->m_spacingZ;

    m_nx = p->m_nx;
    m_ny = p->m_ny;
    m_nz = p->m_nz;

    long n = m_nx*m_ny*m_nz;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = p->m_data[i];
      }
  }


  template<typename T> 
  cArray3D< T > &cArray3D< T >::operator=(cArray3D< T > const &rhs)
  {
    if (this == &rhs)
      {
        return *this;
      }

    // this != &rhs
    if (rhs.m_data) 
      {
        // copy nx, ny, nz
        if (this->m_nx != rhs.m_nx || this->m_ny != rhs.m_ny || this->m_nz != rhs.m_nz)
          {
            this->setSize(rhs.m_nx, rhs.m_ny, rhs.m_nz);
          }

        // copy data
        long n = m_nx*m_ny*m_nz;
        for (long i = 0; i <= n-1; i++)
          {
            this->m_data[i] = rhs.m_data[i];
          }

        // copy spacing
        m_spacingX = rhs.m_spacingX;
        m_spacingY = rhs.m_spacingY;
        m_spacingZ = rhs.m_spacingZ;
      }

    return *this;
  }


  template<typename T> 
  void cArray3D<T>::clear()
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        delete[] m_data;
        m_data = 0;
      }

    return;
  }

  template<typename T> 
  bool cArray3D<T>::setSize(long nx, long ny, long nz)
  {
     if (this->m_data) 
       {
         // m_data already allocated
         if ( (this->m_nx == nx) && (this->m_ny == ny) && (this->m_nz == nz) )
           {
             // if no change in size, do not reallocate.
             return false;
           }
 
         this->clear();
         
         this->m_nx = nx;
         this->m_ny = ny;
         this->m_nz = nz;

         this->m_data = new T[nx*ny*nz];
       }
     else 
       {
         this->m_nx = nx;
         this->m_ny = ny;
         this->m_nz = nz;

         this->m_data = new T[nx*ny*nz];
       }
  
     return true;
  }


  template<typename T> 
  void cArray3D<T>::fill(T val)
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        long n = m_nx*m_ny*m_nz;
        for (long i = 0; i <= n-1; ++i)
          {
            m_data[i] = val;
          }
      }

    return;
  }

  template<typename T> 
  T cArray3D< T >::max_value() const
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }

    T maxValue = m_data[0];
    long n = m_nx*m_ny*m_nz;
    for (long i = 0; i <= n-1; ++i)
      {
        T currentValue = m_data[i];
        maxValue = maxValue>currentValue?maxValue:currentValue;
      }

    return maxValue;
  }

  template<typename T> 
  T cArray3D< T >::min_value() const
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }

    T minValue = m_data[0];
    long n = m_nx*m_ny*m_nz;
    for (long i = 0; i <= n-1; ++i)
      {
        T currentValue = m_data[i];
        minValue = minValue>currentValue?currentValue:minValue;
      }

    return minValue;
  }
  
  template<typename T>
  void cArray3D< T >::extreme_values(T &minValue, T &maxValue) const
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }

    minValue = m_data[0];
    maxValue = minValue;

    long n = m_nx*m_ny*m_nz;
    for (long i = 0; i <= n-1; ++i)
      {
        T currentValue = m_data[i];        
        maxValue = maxValue>currentValue?maxValue:currentValue;
        minValue = minValue>currentValue?currentValue:minValue;
      }

    return;
  }

  // template<typename T>
  // void cArray3D<T>::saveAsMat(const char* fileName)
  // {
  //     /*
  //     Port to vnl for visualization by matlab:

  //     Then in matlab, when load in the .mat file, there will be two
  //     variables: data and size. Do data = reshape(data, size'); then
  //     imagesc(data) to show the data.

  //   */
  //   vnl_vector<T> array_vnl(m_data, m_nx*m_ny*m_nz);
  //   vnl_vector<double> size_vnl(3); // the writer doesn't support int etc. so it's safe to hard code here to avoid the problem when T = int.
  //   size_vnl[0] = m_nx;
  //   size_vnl[1] = m_ny;
  //   size_vnl[2] = m_nz;
  
  //   vnl_matlab_filewrite matlabFileWriter(fileName);
  //   matlabFileWriter.write(array_vnl, "data");
  //   matlabFileWriter.write(size_vnl, "size");
  // }

  // template<typename T>
  // void cArray3D<T>::saveAsNrrd(const char* fileName)
  // {
  //   typedef itk::Image< T , 3 > Function3Type; // 3 means 3D   // I don't know why this doesn't work
  //   typename Function3Type::Pointer finalVolume = Function3Type::New();
  
  //   typename Function3Type::IndexType volStart;
  //   volStart[0] =   0;  // first index on X
  //   volStart[1] =   0;  // first index on Y
  //   volStart[2] =   0;  // first index on Z

  //   typename Function3Type::SizeType  volSize;
  //   volSize[0]  = m_nx;  // size along X
  //   volSize[1]  = m_ny;  // size along Y
  //   volSize[2]  = m_nz;  // size along Z

  //   typename Function3Type::RegionType volRegion;
  //   volRegion.SetSize( volSize );
  //   volRegion.SetIndex( volStart );

  //   finalVolume->SetRegions( volRegion );
  //   finalVolume->Allocate();

  //   typedef itk::ImageRegionIterator< Function3Type > Function3IteratorType;
  //   Function3IteratorType volIt(finalVolume, finalVolume->GetLargestPossibleRegion() );
  //   volIt.GoToBegin();

  //   for (long i = 0; !volIt.IsAtEnd(); ++volIt)
  //     {
  //       volIt.Set(m_data[i]);
  //       ++i;
  //     }


  //   typedef itk::ImageFileWriter< Function3Type > WriterType;
  //   typename WriterType::Pointer writer = WriterType::New();
  //   writer->SetFileName( fileName );
  //   writer->SetInput(finalVolume);

  //   try
  //     {
  //       writer->Update();
  //     }
  //   catch ( itk::ExceptionObject &err)
  //     {
  //       std::cout << "ExceptionObject caught !" << std::endl; 
  //       std::cout << err << std::endl; 
  //       raise(SIGABRT);   
  //     }

  // }



} //namespace 

#endif
