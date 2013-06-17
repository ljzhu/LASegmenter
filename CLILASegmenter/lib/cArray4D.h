#ifndef cArray4D_h_
#define cArray4D_h_

#include "basicHeaders.h"

namespace douher
{
  template<typename T> 
  class cArray4D
  {
  public:
    typedef cArray4D< T > Self;
    typedef boost::shared_ptr< Self > Pointer;
    typedef Pointer PointerType; // just for back compatability


    /* functions */
    cArray4D();
    cArray4D(long nx, long ny, long nz, long n4);
    cArray4D(long nx, long ny, long nz, long n4, T initVal);
    cArray4D(const cArray4D &a);
    ~cArray4D();


    cArray4D< T > &operator=(cArray4D< T > const &rhs);

    inline T const& operator()(long x, long y, long z, long i4) const
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );
      assert(m_n4-1 >= i4 && 0 <= i4 );

      return m_data[i4*m_nz*m_ny*m_nx + z*m_ny*m_nx + y*m_nx + x];
    }

    inline T& operator()(long x, long y, long z, long i4)
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );
      assert(m_n4-1 >= i4 && 0 <= i4 );

      return m_data[i4*m_nz*m_ny*m_nx + z*m_ny*m_nx + y*m_nx + x];
    }


    inline T const& operator[](long idx) const
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz*m_n4-1);

      return m_data[idx];
    }

    inline T& operator[](long idx)
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz*m_n4-1);

      return m_data[idx];
    }


    //////////////////////////////////////////////////
    inline T get(long x, long y, long z, long i4) const
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );
      assert(m_n4-1 >= i4 && 0 <= i4 );

      return m_data[i4*m_nz*m_ny*m_nx + z*m_ny*m_nx + y*m_nx + x];
    }

    inline T get(long idx) const
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz*m_n4-1);

      return m_data[idx];
    }


    void set(long x, long y, long z, long i4, T v)
    {
      assert(m_nx-1 >= x && 0 <= x );
      assert(m_ny-1 >= y && 0 <= y );
      assert(m_nz-1 >= z && 0 <= z );
      assert(m_n4-1 >= i4 && 0 <= i4 );

      m_data[i4*m_nz*m_ny*m_nx + z*m_ny*m_nx + y*m_nx + x] = v;
    }


    void set(long idx, T v)
    {
      assert(idx >= 0 && idx <= m_nx*m_ny*m_nz*m_n4-1);

      m_data[idx] = v;
    }



    long getSizeX() const {return m_nx;}
    long getSizeY() const {return m_ny;}
    long getSizeZ() const {return m_nz;}
    long getSize4() const {return m_n4;}


    double getSpacingX() const {return m_spacingX;}
    double getSpacingY() const {return m_spacingY;}
    double getSpacingZ() const {return m_spacingZ;}
    double getSpacing4() const {return m_spacing4;}

    void setSpacingX(double spX) {m_spacingX = spX;}
    void setSpacingY(double spY) {m_spacingY = spY;}
    void setSpacingZ(double spZ) {m_spacingZ = spZ;}
    void setSpacing4(double sp4) {m_spacing4 = sp4;}


    //    const T* getDataPointer() const {return m_data;}
    T* getDataPointer() const {return m_data;}


    bool setSize(long nx, long ny, long nz, long n4);

    void fill(T val);
    void clear();

    /* data */
  protected:
    T *m_data;
    long m_nx, m_ny, m_nz, m_n4;
    double m_spacingX, m_spacingY, m_spacingZ, m_spacing4;
  };

} // douher

#include "cArray4D.txx"

#endif
/* definition in .txx*/
