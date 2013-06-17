#ifndef cArray4D_txx_
#define cArray4D_txx_

#include "cArray4D.h"

#include <iostream>

namespace douher
{
  template<typename T> 
  cArray4D<T>::cArray4D()
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;
    m_spacing4 = 1.0;

    m_data = 0;

    m_nx = 0;
    m_ny = 0;
    m_nz = 0;
    m_n4 = 0;
  }


  template<typename T> 
  cArray4D<T>::~cArray4D()
  {
    clear();
  }


  template<typename T> 
  cArray4D<T>::cArray4D(long nx, long ny, long nz, long n4)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;
    m_spacing4 = 1.0;

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    m_n4 = n4;

    m_data = new T[nx*ny*nz*n4];
  }

  template<typename T>
  cArray4D<T>::cArray4D(long nx, long ny, long nz, long n4, T initVal)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;
    m_spacingZ = 1.0;
    m_spacing4 = 1.0;

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    m_n4 = n4;

    m_data = new T[nx*ny*nz*n4];

    fill(initVal);
  }


  template<typename T>
  cArray4D<T>::cArray4D(const cArray4D &a)
  {
    m_spacingX = a.m_spacingX;
    m_spacingY = a.m_spacingY;
    m_spacingZ = a.m_spacingZ;
    m_spacing4 = a.m_spacing4;

    m_nx = a.m_nx;
    m_ny = a.m_ny;
    m_nz = a.m_nz;
    m_n4 = a.m_n4;

    long n = m_nx*m_ny*m_nz*m_n4;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = a.m_data[i];
      }
  }


  template<typename T> 
  cArray4D< T > &cArray4D< T >::operator=(cArray4D< T > const &rhs)
  {
    if (this == &rhs)
      {
        return *this;
      }

    // this != &rhs
    if (rhs.m_data) 
      {
        // copy nx, ny, nz, n4
        if (this->m_nx != rhs.m_nx || this->m_ny != rhs.m_ny || this->m_nz != rhs.m_nz || this->m_n4 != rhs.m_n4)
          {
            this->setSize(rhs.m_nx, rhs.m_ny, rhs.m_nz, rhs.m_n4);
          }

        // copy data
        long n = m_nx*m_ny*m_nz*m_n4;
        for (long i = 0; i <= n-1; i++)
          {
            this->m_data[i] = rhs.m_data[i];
          }

        // copy spacing
        m_spacingX = rhs.m_spacingX;
        m_spacingY = rhs.m_spacingY;
        m_spacingZ = rhs.m_spacingZ;
        m_spacing4 = rhs.m_spacing4;
      }

    return *this;
  }



  template<typename T> 
  void cArray4D<T>::clear()
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
  bool cArray4D<T>::setSize(long nx, long ny, long nz, long n4)
  {
     if (this->m_data) 
       {
         // m_data already allocated
         if ( (this->m_nx == nx) && (this->m_ny == ny) && (this->m_nz == nz) && (this->m_n4 == n4) )
           {
             // if no change in size, do not reallocate.
             return false;
           }
 
         this->clear();
         
         this->m_nx = nx;
         this->m_ny = ny;
         this->m_nz = nz;
         this->m_n4 = n4;

         this->m_data = new T[nx*ny*nz*n4];
       }
     else 
       {
         this->m_nx = nx;
         this->m_ny = ny;
         this->m_nz = nz;
         this->m_n4 = n4;

         this->m_data = new T[nx*ny*nz*n4];
       }
  
     return true;
  }


  template<typename T> 
  void cArray4D<T>::fill(T val)
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        long n = m_nx*m_ny*m_nz*m_n4;
        for (long i = 0; i <= n-1; ++i)
          {
            m_data[i] = val;
          }
      }

    return;
  }



  // template<typename T>
  // void cArray4D<T>::saveAsMat(const char* fileName)
  // {
  //     /*
  //     Port to vnl for visualization by matlab:

  //     Then in matlab, when load in the .mat file, there will be two
  //     variables: data and size. Do data = reshape(data, size'); then
  //     imagesc(data) to show the data.

  //   */
  //   vnl_vector<T> array_vnl(m_data, m_nx*m_ny*m_nz*m_n4);
  //   vnl_vector<double> size_vnl(4); // the writer doesn't support int etc. so it's safe to hard code here to avoid the problem when T = int.
  //   size_vnl[0] = m_nx;
  //   size_vnl[1] = m_ny;
  //   size_vnl[2] = m_nz;
  //   size_vnl[3] = m_n4;
  
  //   vnl_matlab_filewrite matlabFileWriter(fileName);
  //   matlabFileWriter.write(array_vnl, "data");
  //   matlabFileWriter.write(size_vnl, "size");
  // }

  // template<typename T>
  // void cArray4D<T>::saveAsNrrd(const char* fileName)
  // {
  //   //  typedef itk::Image< T , 4 > Function4Type; // 4 means 4D   // I don't know why this doesn't work
  //   typedef itk::Image< double , 4 > Function4Type; // 4 means 4D  
  //   Function4Type::Pointer finalVolume = Function4Type::New();
  
  //   Function4Type::IndexType volStart;
  //   volStart[0] =   0;  // first index on X
  //   volStart[1] =   0;  // first index on Y
  //   volStart[2] =   0;  // first index on Z

  //   Function4Type::SizeType  volSize;
  //   volSize[0]  = m_nx;  // size along X
  //   volSize[1]  = m_ny;  // size along Y
  //   volSize[2]  = m_nz;  // size along Z

  //   Function4Type::RegionType volRegion;
  //   volRegion.SetSize( volSize );
  //   volRegion.SetIndex( volStart );

  //   finalVolume->SetRegions( volRegion );
  //   finalVolume->Allocate();

  //   typedef itk::ImageRegionIterator< Function4Type > Function4IteratorType;
  //   Function4IteratorType volIt(finalVolume, finalVolume->GetLargestPossibleRegion() );
  //   volIt.GoToBegin();

  //   for (long i = 0; !volIt.IsAtEnd(); ++volIt)
  //     {
  //       volIt.Set(m_data[i]);
  //       ++i;
  //     }


  //   typedef itk::ImageFileWriter< Function4Type > WriterType;
  //   WriterType::Pointer writer = WriterType::New();
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

} // douher

#endif
