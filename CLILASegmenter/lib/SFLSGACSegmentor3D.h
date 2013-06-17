#ifndef SFLSGACSegmentor3D_h_
#define SFLSGACSegmentor3D_h_

#include "SFLSSegmentor3D.h"

//#include <list>


//douher
#include "cArray3D.h"

namespace douher
{
  template< typename TPixel >
  class CSFLSGACSegmentor3D : public CSFLSSegmentor3D< TPixel >
  {
  public:
    typedef CSFLSSegmentor3D< TPixel > SuperClassType;

    typedef CSFLSGACSegmentor3D< TPixel > Self;
    typedef boost::shared_ptr< Self > Pointer;
    typedef Pointer PointerType; // for backward compatability


    /*================================================================================
      ctor */
    CSFLSGACSegmentor3D() : CSFLSSegmentor3D< TPixel >()
    {
      basicInit();
    }

    void basicInit();

    /* New */
    static Pointer New() { return Pointer(new Self); }

    void setCurvatureWeight(double a) 
    {
      m_curvatureWeight = a;
    }

    void setInflation(double a)
    {
      m_inflation = a;
    }

    void setFeatureImage(cArray3D< double >::PointerType fImg);

    // data
    cArray3D< double >::PointerType mp_featureImage;


    double m_curvatureWeight;
    double m_inflation;


    /* ============================================================
     * functions
     * ============================================================*/
    void doGACSegmenation();
    


    /* ============================================================
       computeForce    */
    void computeForce();
    double computeForceAt(long ix, long iy, long iz);

  };

}//douher

#include "SFLSGACSegmentor3D.hpp"

#endif
