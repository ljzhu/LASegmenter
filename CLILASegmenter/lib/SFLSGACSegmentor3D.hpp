#ifndef SFLSGACSegmentor3D_hpp_
#define SFLSGACSegmentor3D_hpp_

#include "SFLSGACSegmentor3D.h"

#include <algorithm>


//douher
#include "cArrayOp.h"

namespace douher
{
  /* ============================================================
     basicInit    */
  template< typename TPixel >
  void
  CSFLSGACSegmentor3D< TPixel >
  ::basicInit()
  {
    SuperClassType::basicInit();
    
    m_curvatureWeight = 0.1;
    m_inflation = 1.0; // neg: shrink, pos: dilation
  }  


  /* ============================================================
     setFeatureImage    */
  template< typename TPixel >
  void
  CSFLSGACSegmentor3D< TPixel >
  ::setFeatureImage(cArray3D< double >::PointerType fImg)
  {
    mp_featureImage = fImg;

    if (this->m_nx + this->m_ny + this->m_nz == 0)
      {
        this->m_nx = fImg->getSizeX();
        this->m_ny = fImg->getSizeY();
        this->m_nz = fImg->getSizeZ();
      }
    else if ( this->m_nx != fImg->getSizeX() || this->m_ny != fImg->getSizeY() || this->m_nz != fImg->getSizeZ() )
      {
        std::cerr<<"image sizes do not match. abort\n";
        raise(SIGABRT);
      }
  }


  /* ============================================================
     computeForce    */
  template< typename TPixel >
  void
  CSFLSGACSegmentor3D< TPixel >
  ::computeForce()
  {
    this->m_force.clear();

    double fmax = -1e10;

    long n = this->m_lz->size();
    //    double* kappaOnZeroLS = new double[ n ];
    double* theForce = new double[ n ];

    {
      long i = 0;
      for (CSFLSLayer::iterator itz = this->m_lz->begin(); itz != this->m_lz->end(); ++itz, ++i)
        {
          long ix = (*itz)[0];
          long iy = (*itz)[1];
          long iz = (*itz)[2];

          double f = computeForceAt(ix, iy, iz);
        
          fmax = fabs(f)>fmax?fabs(f):fmax;

          theForce[i] = f;
        }
    }


    for (long i = 0; i < n; ++i)
      {
        this->m_force.push_back(theForce[i]/(fmax + 1e-10));
      }
    
    delete[] theForce;
  }


  /* ============================================================
     doChanVeseSegmenation    */
  template< typename TPixel >
  void
  CSFLSGACSegmentor3D< TPixel >
  ::doGACSegmenation()
  {
    if (!mp_featureImage)
      {
	std::cerr<<"feature image not ready. stop.\n";
	raise(SIGABRT);
      }

    /*============================================================
     * From the initial mask, generate: 1. SFLS, 2. mp_label and
     * 3. mp_phi.      
     */
    this->initializeSFLS();


    //douher::saveAsImage3< double >(mp_phi, "initPhi.nrrd");
    for (unsigned int it = 0; it < this->m_numIter; ++it)
      {
        computeForce();

        this->normalizeForce();

        this->oneStepLevelSetEvolution();
      }
  }

  /*============================================================

  computeUpwindGradAt 

  Compute upwind gradient of the level set fn at a point */
  template< typename TPixel >
  double
  CSFLSGACSegmentor3D< TPixel >
  ::computeForceAt(long ix, long iy, long iz)
  {
    /*----------------------------------------------------------------------
      According to eq(8) in Yezzi GAC TMI paper. */

    double f = 0;
    
    double dPhiDxMinus = 0.0;
    double dPhiDxPlus = 0.0;
    double dPhiDyMinus = 0.0;
    double dPhiDyPlus = 0.0;
    double dPhiDzMinus = 0.0;
    double dPhiDzPlus = 0.0;

    double dFeatureDx = 0.0;
    double dFeatureDy = 0.0;
    double dFeatureDz = 0.0;

    double featureVal = mp_featureImage->get(ix, iy, iz);


    if( ix+1 < this->m_nx && ix-1 >=0 )
      {
	dPhiDxMinus = this->mp_phi->get(ix, iy, iz) - this->mp_phi->get(ix-1, iy, iz);
	dPhiDxPlus = this->mp_phi->get(ix+1, iy, iz) - this->mp_phi->get(ix, iy, iz);

        dFeatureDx = (this->mp_featureImage->get(ix+1, iy, iz) - this->mp_featureImage->get(ix-1, iy, iz))/2.0;
      }

    if( iy+1 < this->m_ny && iy-1 >=0 ) 
      {
	dPhiDyMinus = this->mp_phi->get(ix, iy, iz) - this->mp_phi->get(ix, iy-1, iz);
	dPhiDyPlus = this->mp_phi->get(ix, iy+1, iz) - this->mp_phi->get(ix, iy, iz);

        dFeatureDy = (this->mp_featureImage->get(ix, iy+1, iz) - this->mp_featureImage->get(ix, iy-1, iz))/2.0;
      }

    if( iz+1 < this->m_nz && iz-1 >=0 )
      {
	dPhiDzMinus = this->mp_phi->get(ix, iy, iz) - this->mp_phi->get(ix, iy, iz-1);
	dPhiDzPlus = this->mp_phi->get(ix, iy, iz+1) - this->mp_phi->get(ix, iy, iz);

        dFeatureDz = (this->mp_featureImage->get(ix, iy, iz+1) - this->mp_featureImage->get(ix, iy, iz-1))/2.0;
      }


    // \gradient(feature) dotProd \upwind_gradient(\phi)
    double dGdxdPhidx = dFeatureDx>0?(dFeatureDx*dPhiDxPlus):(dFeatureDx*dPhiDxMinus);
    double dGdydPhidy = dFeatureDy>0?(dFeatureDy*dPhiDyPlus):(dFeatureDy*dPhiDyMinus);
    double dGdzdPhidz = dFeatureDz>0?(dFeatureDz*dPhiDzPlus):(dFeatureDz*dPhiDzMinus);

    double gradFeatureDotGradPhi = dGdxdPhidx*dGdydPhidy*dGdzdPhidz;

    f += gradFeatureDotGradPhi;

    double k = this->computeKappa(ix, iy, iz);

    f += featureVal*(m_curvatureWeight*k - m_inflation);

    return f;
  }


}

#endif
