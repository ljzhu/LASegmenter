#ifndef cArray2D_txx_
#define cArray2D_txx_

#include "cArray2D.h"

#include <iostream>

namespace douher
{
  template<typename T> 
  cArray2D<T>::cArray2D()
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;

    m_data = 0;

    m_nx = 0;
    m_ny = 0;
  }


  template<typename T> 
  cArray2D<T>::~cArray2D()
  {
    clear();
  }


  template<typename T> 
  cArray2D<T>::cArray2D(long nx, long ny)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;

    m_nx = nx;
    m_ny = ny;

    m_data = new T[nx*ny];
  }

  template<typename T>
  cArray2D<T>::cArray2D(long nx, long ny, T initVal)
  {
    m_spacingX = 1.0;
    m_spacingY = 1.0;

    m_nx = nx;
    m_ny = ny;

    m_data = new T[nx*ny];

    fill(initVal);
  }


  template<typename T>
  cArray2D<T>::cArray2D(const cArray2D& a)
  {
    m_spacingX = a.m_spacingX;
    m_spacingY = a.m_spacingY;

    m_nx = a.m_nx;
    m_ny = a.m_ny;

    long n = m_nx*m_ny;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = a.m_data[i];
      }
  }

  template<typename T>
  cArray2D<T>::cArray2D(const PointerType p)
  {
    m_spacingX = p->m_spacingX;
    m_spacingY = p->m_spacingY;

    m_nx = p->m_nx;
    m_ny = p->m_ny;

    long n = m_nx*m_ny;

    m_data = new T[n];

    for (long i = 0; i < n; ++i)
      {
        this->m_data[i] = p->m_data[i];
      }
  }

  template<typename T> 
  cArray2D< T >& cArray2D< T >::operator=(cArray2D< T > const& rhs)
  {
    if (this == &rhs)
      {
        return *this;
      }

    // this != &rhs
    if (rhs.m_data) 
      {
        // copy nx, ny
        if (this->m_nx != rhs.m_nx || this->m_ny != rhs.m_ny)
          {
            this->setSize(rhs.m_nx, rhs.m_ny);
          }

        // copy data
        long n = m_nx*m_ny;
        for (long i = 0; i <= n-1; i++)
          {
            this->m_data[i] = rhs.m_data[i];
          }

        // copy spacing
        m_spacingX = rhs.m_spacingX;
        m_spacingY = rhs.m_spacingY;
      }

    return *this;
  }



  template<typename T> 
  void cArray2D<T>::clear()
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

        m_nx = 0;
        m_ny = 0;

        m_spacingX = 0;
        m_spacingY = 0;
      }

    return;
  }


  template<typename T> 
  bool cArray2D<T>::setSize(long nx, long ny)
  {
     if (this->m_data) 
       {
         // m_data already allocated
         if ( (this->m_nx == nx) && (this->m_ny == ny) )
           {
             // if no change in size, do not reallocate.
             return false;
           }
 
         this->clear();
         
         this->m_nx = nx;
         this->m_ny = ny;

         this->m_data = new T[nx*ny];
       }
     else 
       {
         this->m_nx = nx;
         this->m_ny = ny;

         this->m_data = new T[nx*ny];
       }
  
     return true;
  }


  template<typename T> 
  void cArray2D<T>::fill(T val)
  {
    if (!m_data)
      {
        std::cerr<<"memory not allocated. exit.\n";
        raise(SIGABRT);
      }
    else
      {
        long n = m_nx*m_ny;
        for (long i = 0; i <= n-1; ++i)
          {
            m_data[i] = val;
          }
      }

    return;
  }

  // template<typename T>
  // void cArray2D<T>::saveAsMat(const char* fileName)
  // {
  //     /*
  //     Port to vnl for visualization by matlab:

  //     Then in matlab, when load in the .mat file, there will be two
  //     variables: data and size. Do data = reshape(data, size'); then
  //     imagesc(data) to show the data.

  //   */
  //   vnl_vector<T> array_vnl(m_data, m_nx*m_ny);
  //   vnl_vector<double> size_vnl(2); // the writer doesn't support int etc. so it's safe to hard code here to avoid the problem when T = int.
  //   size_vnl[0] = m_nx;
  //   size_vnl[1] = m_ny;
  
  //   vnl_matlab_filewrite matlabFileWriter(fileName);
  //   matlabFileWriter.write(array_vnl, "data");
  //   matlabFileWriter.write(size_vnl, "size");
  // }

  // template<typename T>
  // void cArray2D<T>::saveAsNrrd(const char* fileName)
  // {
  //   typedef itk::Image< T , 2 > Function2Type; // 2 means 2D  

  //   typename Function2Type::Pointer finalVolume = Function2Type::New();
  
  //   typename Function2Type::IndexType volStart;
  //   volStart[0] =   0;  // first index on X
  //   volStart[1] =   0;  // first index on Y

  //   typename Function2Type::SizeType  volSize;
  //   volSize[0]  = m_nx;  // size along X
  //   volSize[1]  = m_ny;  // size along Y

  //   typename Function2Type::RegionType volRegion;
  //   volRegion.SetSize( volSize );
  //   volRegion.SetIndex( volStart );

  //   finalVolume->SetRegions( volRegion );
  //   finalVolume->Allocate();

  //   typedef itk::ImageRegionIterator< Function2Type > Function2IteratorType;
  //   Function2IteratorType volIt(finalVolume, finalVolume->GetLargestPossibleRegion() );
  //   volIt.GoToBegin();

  //   for (long i = 0; !volIt.IsAtEnd(); ++volIt)
  //     {
  //       volIt.Set(m_data[i]);
  //       ++i;
  //     }


  //   typedef itk::ImageFileWriter< Function2Type > WriterType;
  // //   typedef typename WriterType::Pointer WriterPointerType;

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

} //douher

#endif
