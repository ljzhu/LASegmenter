#ifndef SFLSSegmentor3D_h_
#define SFLSSegmentor3D_h_

#include "SFLS.h"

#include <list>

//douher
#include "cArray3D.h"

namespace douher
{
  template< typename TPixel >
  class CSFLSSegmentor3D : public CSFLS
  {
  public:
    virtual ~CSFLSSegmentor3D() {}

    typedef CSFLSSegmentor3D< TPixel > Self;
    typedef boost::shared_ptr< Self > Pointer;

    typedef cArray3D< TPixel > ImageType;
    typedef cArray3D< double > LSImageType;
    typedef cArray3D< unsigned char > MaskImageType;

    CSFLSSegmentor3D();

    /* New */
    static Pointer New() { return Pointer(new Self); }


    /* ============================================================
     * functions         */
    void basicInit();

    void setNumIter(unsigned long n);

    void setImage(typename ImageType::Pointer img);
    void setMask(typename MaskImageType::Pointer mask);

    virtual void computeForce() = 0;

    void normalizeForce();

//     double maxPhi(long ix, long iy, long iz, double level);
//     double minPhi(long ix, long iy, long iz, double level);
    bool getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long ix, long iy, long iz, double& thePhi);

    void oneStepLevelSetEvolution();

    void getSFLSFromPhi();

    void initializeSFLS() { initializeSFLSFromMask(); }
    void initializeSFLSFromMask();

    void initializeLabel();
    void initializePhi();

    void doSegmenation();


    // geometry
    double computeKappa(long ix, long iy, long iz);
    



    /* ============================================================
     * data     */
    // CSFLS::PointerType mp_sfls;

    typename cArray3D< TPixel >::PointerType mp_img;
    typename cArray3D< char >::PointerType mp_label;
    typename cArray3D< unsigned char >::PointerType mp_mask;
    typename cArray3D< double >::PointerType mp_phi;

    long m_nx;
    long m_ny;
    long m_nz;


    std::list< double > m_force;

    double m_timeStep;

    unsigned long m_numIter;

  protected:
    inline bool doubleEqual(double a, double b, double eps = 1e-10)
    {
      return (a-b < eps && b-a < eps);
    }

    /*----------------------------------------------------------------------
     These two record the pts which change status 

    Because they are created and visited sequentially, and when not
    needed, are clear-ed as a whole. No random insertion or removal is
    needed. So use vector is faster than list.  */ 
    CSFLSLayer m_lIn2out;
    CSFLSLayer m_lOut2in;
  };


} // douher


#include "SFLSSegmentor3D.hpp"


#endif
