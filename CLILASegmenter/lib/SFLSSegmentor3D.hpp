#ifndef SFLSSegmentor3D_hpp_
#define SFLSSegmentor3D_hpp_

#include "SFLSSegmentor3D.h"

#include <algorithm>
#include <cmath>


//douher
#include "cArrayOp.h"


namespace douher
{
  template< typename TPixel >
  CSFLSSegmentor3D< TPixel >
  ::CSFLSSegmentor3D() : CSFLS()
  {
    basicInit(); 
  }


  /* ============================================================
     basicInit    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::basicInit()
  {
    m_numIter = 100;
    m_timeStep = 1.0;


    m_nx = 0;
    m_ny = 0;
    m_nz = 0;
  }

  /* ============================================================
     setImage    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::setNumIter(unsigned long n)
  {
    m_numIter = n;
  }

  /* ============================================================
     setImage    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::setImage(typename ImageType::Pointer img)
  {
    mp_img = img;

    if (m_nx + m_ny + m_nz == 0)
      {
        m_nx = img->getSizeX();
        m_ny = img->getSizeY();
        m_nz = img->getSizeZ();
      }
    else if ( m_nx != img->getSizeX() || m_ny != img->getSizeY() || m_nz != img->getSizeZ() )
      {
        std::cerr<<"image sizes do not match. abort\n";
        raise(SIGABRT);
      }
  }

  /* ============================================================
     setMask    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::setMask(typename MaskImageType::Pointer mask)
  {
    mp_mask = mask;

    if (m_nx + m_ny + m_nz == 0)
      {
        m_nx = mask->getSizeX();
        m_ny = mask->getSizeY();
        m_nz = mask->getSizeZ();
      }
    else if ( m_nx != mask->getSizeX() || m_ny != mask->getSizeY() || m_nz != mask->getSizeZ() )
      {
        std::cerr<<"image sizes do not match. abort\n";
        raise(SIGABRT);
      }
  }



  template< typename TPixel >
  bool
  CSFLSSegmentor3D< TPixel >
  ::getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long ix, long iy, long iz, double& thePhi)
  {
    /*--------------------------------------------------
     *
     * Look in all the neighbors, to find the phi value of the nbhd:
     * this nbhd should satisfy: 1. its layer is strictly closer to
     * the zero layer hence its value is thought to be updated. 2. If
     * there are several nbhd's belonging to the same layer, choose
     * the one whose phi value has the smallest abs value.  If (ix,
     * iy) is outside, go through all nbhd who is in the layer of
     * label = mylevel-1 pick the SMALLEST phi. If (ix, iy) is inside,
     * go through all nbhd who is in the layer of label = mylevel+1
     * pick the LARGEST phi.
     */

    char mylevel = mp_label->get(ix, iy, iz);
    bool foundNbhd = false;

    if (mylevel > 0)
      {
        // find the SMALLEST phi
        thePhi = 10000;

        if ( (ix+1 < m_nx) && (mp_label->get(ix+1, iy, iz) == mylevel - 1) )
          {
            double itsPhi = mp_phi->get(ix+1, iy, iz);
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( (ix-1 >= 0 ) && (mp_label->get(ix-1, iy, iz) == mylevel - 1) )
          {
            double itsPhi = mp_phi->get(ix-1, iy, iz);
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( (iy+1 < m_ny) && (mp_label->get(ix, iy+1, iz) == mylevel - 1) )
          {
            double itsPhi = mp_phi->get(ix, iy+1, iz);
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( ((iy-1) >= 0 ) && (mp_label->get(ix, iy-1, iz) == mylevel - 1) ) 
          {
            double itsPhi = mp_phi->get(ix, iy-1, iz);
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( (iz+1 < m_nz) && (mp_label->get(ix, iy, iz+1) == mylevel - 1) )
          {
            double itsPhi = mp_phi->get(ix, iy, iz+1);
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( ((iz-1) >= 0 ) && (mp_label->get(ix, iy, iz-1) == mylevel - 1) ) 
          {
            double itsPhi = mp_phi->get(ix, iy, iz-1);
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }


      }
    else
      {
        // find the LARGEST phi
        thePhi = -10000;

        if ( (ix+1 < m_nx) && (mp_label->get(ix+1, iy, iz) == mylevel + 1) )
          {
            double itsPhi = mp_phi->get(ix+1, iy, iz);
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }


        if ( (ix-1 >= 0  ) && (mp_label->get(ix-1, iy, iz) == mylevel + 1) )
          {
            double itsPhi = mp_phi->get(ix-1, iy, iz);
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( (iy+1 < m_ny) && (mp_label->get(ix, iy+1, iz) == mylevel + 1) )
          {
            double itsPhi = mp_phi->get(ix, iy+1, iz);
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( ((iy-1) >= 0 ) && (mp_label->get(ix, iy-1, iz) == mylevel + 1) ) 
          {
            double itsPhi = mp_phi->get(ix, iy-1, iz);
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( (iz+1 < m_nz) && (mp_label->get(ix, iy, iz+1) == mylevel + 1) )
          {
            double itsPhi = mp_phi->get(ix, iy, iz+1);
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }

        if ( ((iz-1) >= 0 ) && (mp_label->get(ix, iy, iz-1) == mylevel + 1) ) 
          {
            double itsPhi = mp_phi->get(ix, iy, iz-1);
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
          }
      }



    return foundNbhd;
  }


  /* ============================================================
     normalizeForce   
     Normalize m_force s.t. max(abs(m_force)) < 0.5 */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::normalizeForce()
  {
    unsigned long nLz = m_lz->size();

    if (m_force.size() != nLz )
      {
        std::cerr<<"m_force.size() = "<<m_force.size()<<std::endl;
        std::cerr<<"nLz = "<<nLz<<std::endl;

        std::cerr<<"m_force.size() != nLz, abort.\n";
        raise(SIGABRT);
      }

    double fMax = fabs( m_force.front() );

    for (std::list<double>::const_iterator itf = m_force.begin(); itf != m_force.end(); ++itf)
      {
        double v = fabs(*itf);
        fMax = fMax>v?fMax:v;
      }

    fMax /= 0.49;

    for (std::list<double>::iterator itf = m_force.begin(); itf != m_force.end(); ++itf)
      {
        (*itf) /= (fMax + 1e-10);
      }
  }

  /* ============================================================
     oneStepLevelSetEvolution    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::oneStepLevelSetEvolution()
  {
    // create 'changing status' lists
    CSFLSLayer Sz;
    CSFLSLayer Sn1;
    CSFLSLayer Sp1;
    CSFLSLayer Sn2;
    CSFLSLayer Sp2;


    m_lIn2out.clear();
    m_lOut2in.clear();

    /*--------------------------------------------------
      1. add F to phi(Lz), create Sn1 & Sp1 
      scan Lz values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
                                          ========                */
    {
      std::list<double>::const_iterator itf = m_force.begin();
      for (CSFLSLayer::iterator itz = m_lz->begin(); \
           itz != m_lz->end(); \
           ++itf)
        {
          long ix = (*itz)[0];
          long iy = (*itz)[1];
          long iz = (*itz)[2];

          double phi_old = mp_phi->get(ix, iy, iz);
          double phi_new = phi_old + (*itf);

	  /*----------------------------------------------------------------------
	   Update the lists of pt who change the state, for faster
	   energy fnal computation. */
	  if ( phi_old <= 0 && phi_new > 0 )
	    {
	      m_lIn2out.push_back(CSFLSLayer::NodeType(ix, iy, iz));
	    }

	  if( phi_old>0  && phi_new <= 0)
	    {
	      m_lOut2in.push_back(CSFLSLayer::NodeType(ix, iy, iz));
	    }

//           // DEBUG
//           if (phi_new > 3.1 || phi_new < -3.1)
//             {
//               std::cout<<"phi_old = "<<phi_old<<std::endl;
//               std::cout<<"its lbl = "<<(int)mp_label->get(ix, iy)<<std::endl;

//               std::cerr<<"phi_new > 3.1 || phi_new < -3.1\n";
//               raise(SIGABRT);
//             }

          mp_phi->set(ix, iy, iz, phi_new);

          if(phi_new > 0.5)
            {
              Sp1.push_back(*itz);
              itz = m_lz->erase(itz);
            }
          else if (phi_new < -0.5)
            {
              Sn1.push_back(*itz);
              itz = m_lz->erase(itz);
            }
          else
            {
              ++itz;              
            }

          /*--------------------------------------------------
            NOTE, mp_label are (should) NOT update here. They should
            be updated with Sz, Sn/p's
            --------------------------------------------------*/
        }
    }


//     // debug
//     labelsCoherentCheck1();


    /*--------------------------------------------------
      2. update Ln1,Lp1,Lp2,Lp2, ****in that order****

      2.1 scan Ln1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5] 
                                     ==========                     */
    for (CSFLSLayer::iterator itn1 = m_ln1->begin(); itn1 != m_ln1->end(); )
      {
        long ix = (*itn1)[0];
        long iy = (*itn1)[1];
        long iz = (*itn1)[2];

        double thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, thePhi);

        if (found)
          {
            double phi_new = thePhi-1;
            mp_phi->set(ix, iy, iz, phi_new);

            if (phi_new >= -0.5)
              {
                Sz.push_back(*itn1);
                itn1 = m_ln1->erase(itn1);
              }
            else if (phi_new < -1.5)
              {
                Sn2.push_back(*itn1);
                itn1 = m_ln1->erase(itn1);
              }
            else
              {
                ++itn1;
              }
          }
        else 
          {
            /*--------------------------------------------------
              No nbhd in inner (closer to zero contour) layer, so
              should go to Sn2. And the phi shold be further -1
            */
            Sn2.push_back(*itn1);
            itn1 = m_ln1->erase(itn1);

            mp_phi->set(ix, iy, iz, mp_phi->get(ix, iy, iz) - 1);  
          }
      }



//     // debug
//     labelsCoherentCheck1();


    /*--------------------------------------------------
      2.2 scan Lp1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5] 
                                                       ========          */
    for (CSFLSLayer::iterator itp1 = m_lp1->begin(); itp1 != m_lp1->end();)
      {
        long ix = (*itp1)[0];
        long iy = (*itp1)[1];
        long iz = (*itp1)[2];

        double thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, thePhi);

        if (found)
          {
            double phi_new = thePhi+1;
            mp_phi->set(ix, iy, iz, phi_new);

            if (phi_new <= 0.5)
              {
                Sz.push_back(*itp1);
                itp1 = m_lp1->erase(itp1);
              }
            else if (phi_new > 1.5)
              {
                Sp2.push_back(*itp1);
                itp1 = m_lp1->erase(itp1);
              }
            else
              {
                ++itp1;
              }
          }
        else 
          {
            /*--------------------------------------------------
              No nbhd in inner (closer to zero contour) layer, so
              should go to Sp2. And the phi shold be further +1
            */

            Sp2.push_back(*itp1);
            itp1 = m_lp1->erase(itp1);

            mp_phi->set(ix, iy, iz, mp_phi->get(ix, iy, iz) + 1);  
          }
      }


//     // debug
//     labelsCoherentCheck1();



    /*--------------------------------------------------
      2.3 scan Ln2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5] 
      ==========                                      */
    for (CSFLSLayer::iterator itn2 = m_ln2->begin(); itn2 != m_ln2->end(); )
      {
        long ix = (*itn2)[0];
        long iy = (*itn2)[1];
        long iz = (*itn2)[2];

        double thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, thePhi);

        if (found)
          {
            double phi_new = thePhi-1;
            mp_phi->set(ix, iy, iz, phi_new);

            if (phi_new >= -1.5)
              {
                Sn1.push_back(*itn2);
                itn2 = m_ln2->erase(itn2);
              }
            else if (phi_new < -2.5)
              {
                itn2 = m_ln2->erase(itn2);
                mp_phi->set(ix, iy, iz, -3);
                mp_label->set(ix, iy, iz, -3);
              }
            else
              {
                ++itn2;
              }
          }
        else 
          {
            itn2 = m_ln2->erase(itn2);
            mp_phi->set(ix, iy, iz, -3);
            mp_label->set(ix, iy, iz, -3);
          }
      }


//     // debug
//     labelsCoherentCheck1();



    /*--------------------------------------------------
      2.4 scan Lp2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
      ========= */
    for (CSFLSLayer::iterator itp2 = m_lp2->begin(); itp2 != m_lp2->end(); )
      {
        long ix = (*itp2)[0];
        long iy = (*itp2)[1];
        long iz = (*itp2)[2];

        double thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, thePhi);

        if (found)
          {
            double phi_new = thePhi+1;
            mp_phi->set(ix, iy, iz, phi_new);

            if (phi_new <= 1.5)
              {
                Sp1.push_back(*itp2);
                itp2 = m_lp2->erase(itp2);
              }
            else if (phi_new > 2.5)
              {
                itp2 = m_lp2->erase(itp2);
                mp_phi->set(ix, iy, iz, 3);
                mp_label->set(ix, iy, iz, 3);
              }
            else
              {
                ++itp2;
              }
          }
        else 
          {
            itp2 = m_lp2->erase(itp2);
            mp_phi->set(ix, iy, iz, 3);
            mp_label->set(ix, iy, iz, 3);
          }
      }


//     // debug
//     labelsCoherentCheck1();



    /*--------------------------------------------------
      3. Deal with S-lists Sz,Sn1,Sp1,Sn2,Sp2 
      3.1 Scan Sz */
    for (CSFLSLayer::iterator itSz = Sz.begin(); itSz != Sz.end(); ++itSz)
      {
        long ix = (*itSz)[0];
        long iy = (*itSz)[1];
        long iz = (*itSz)[2];
        
        m_lz->push_back(*itSz);
        mp_label->set(ix, iy, iz, 0);
      }


//     // debug
//     labelsCoherentCheck1();



    /*--------------------------------------------------
      3.2 Scan Sn1     */
    for (CSFLSLayer::iterator itSn1 = Sn1.begin(); itSn1 != Sn1.end(); ++itSn1)
      {
        long ix = (*itSn1)[0];
        long iy = (*itSn1)[1];
        long iz = (*itSn1)[2];
        
        m_ln1->push_back(*itSn1);
        // itSn1 = Sn1.erase(itSn1);

        mp_label->set(ix, iy, iz, -1);

        if ( (ix+1 < m_nx) && doubleEqual(mp_phi->get(ix+1, iy, iz), -3.0) )
          {
            Sn2.push_back(CSFLSLayer::NodeType(ix+1, iy, iz));
            mp_phi->set(ix+1, iy, iz, mp_phi->get(ix, iy, iz) - 1 ); 
          }

        if ( (ix-1 >= 0) && doubleEqual(mp_phi->get(ix-1, iy, iz), -3.0) )
          {
            Sn2.push_back(CSFLSLayer::NodeType(ix-1, iy, iz));
            mp_phi->set(ix-1, iy, iz, mp_phi->get(ix, iy, iz) - 1); 
          }

        if ( (iy+1 < m_ny) && doubleEqual(mp_phi->get(ix, iy+1, iz), -3.0) )
          {
            Sn2.push_back(CSFLSLayer::NodeType(ix, iy+1, iz));
            mp_phi->set(ix, iy+1, iz, mp_phi->get(ix, iy, iz) - 1 ); 
          }

	if ( (iy-1>=0) && doubleEqual(mp_phi->get(ix, iy-1, iz), -3.0) )
          {
            Sn2.push_back(CSFLSLayer::NodeType(ix, iy-1, iz) );
            mp_phi->set(ix, iy-1, iz, mp_phi->get(ix, iy, iz) - 1 );
          }

        if ( (iz+1 < m_nz) && doubleEqual(mp_phi->get(ix, iy, iz+1), -3.0) )
          {
            Sn2.push_back(CSFLSLayer::NodeType(ix, iy, iz+1));
            mp_phi->set(ix, iy, iz+1, mp_phi->get(ix, iy, iz) - 1 ); 
          }

        if ( (iz-1>=0) && doubleEqual(mp_phi->get(ix, iy, iz-1), -3.0) )
          {
            Sn2.push_back(CSFLSLayer::NodeType(ix, iy, iz-1) );
            mp_phi->set(ix, iy, iz-1, mp_phi->get(ix, iy, iz) - 1 );
          }
      }


//     // debug
//     labelsCoherentCheck1();


    /*--------------------------------------------------
      3.3 Scan Sp1     */
    for (CSFLSLayer::iterator itSp1 = Sp1.begin(); itSp1 != Sp1.end(); ++itSp1)
      {
        long ix = (*itSp1)[0];
        long iy = (*itSp1)[1];
        long iz = (*itSp1)[2];
        
        m_lp1->push_back(*itSp1);
        mp_label->set(ix, iy, iz, 1);

        if ( (ix+1 < m_nx) && doubleEqual(mp_phi->get(ix+1, iy, iz), 3.0) )
          {
            Sp2.push_back(CSFLSLayer::NodeType(ix+1, iy, iz));
            mp_phi->set(ix+1, iy, iz, mp_phi->get(ix, iy, iz) + 1 ); 
          }

        if ( (ix-1 >= 0) && doubleEqual(mp_phi->get(ix-1, iy, iz), 3.0) )
          {
            Sp2.push_back(CSFLSLayer::NodeType(ix-1, iy, iz));
            mp_phi->set(ix-1, iy, iz, mp_phi->get(ix, iy, iz) + 1); 
          }

        if ( (iy+1 < m_ny) && doubleEqual(mp_phi->get(ix, iy+1, iz), 3.0) )
          {
            Sp2.push_back(CSFLSLayer::NodeType(ix, iy+1, iz));
            mp_phi->set(ix, iy+1, iz, mp_phi->get(ix, iy, iz) + 1 ); 
          }

        if ( (iy-1>=0) && doubleEqual(mp_phi->get(ix, iy-1, iz), 3.0) )
          {
            Sp2.push_back(CSFLSLayer::NodeType(ix, iy-1, iz) );
            mp_phi->set(ix, iy-1, iz, mp_phi->get(ix, iy, iz) + 1 );
          }

        if ( (iz+1 < m_nz) && doubleEqual(mp_phi->get(ix, iy, iz+1), 3.0) )
          {
            Sp2.push_back(CSFLSLayer::NodeType(ix, iy, iz+1));
            mp_phi->set(ix, iy, iz+1, mp_phi->get(ix, iy, iz) + 1 ); 
          }


        if ( (iz-1>=0) && doubleEqual(mp_phi->get(ix, iy, iz-1), 3.0) )
          {
            Sp2.push_back(CSFLSLayer::NodeType(ix, iy, iz-1) );
            mp_phi->set(ix, iy, iz-1, mp_phi->get(ix, iy, iz) + 1 );
          }
      } 


//     // debug
//     labelsCoherentCheck1();


    /*--------------------------------------------------
      3.4 Scan Sn2     */
    {
      //debug 
      int aaa = 0;
      for (CSFLSLayer::iterator itSn2 = Sn2.begin(); itSn2 != Sn2.end(); ++itSn2, ++aaa)
        {
          long ix = (*itSn2)[0];
          long iy = (*itSn2)[1];
          long iz = (*itSn2)[2];
        
          m_ln2->push_back(*itSn2);

          mp_label->set(ix, iy, iz, -2);

//           // debug
//           labelsCoherentCheck1();
        }
    }



    /*--------------------------------------------------
      3.5 Scan Sp2     */
    for (CSFLSLayer::iterator itSp2 = Sp2.begin(); itSp2 != Sp2.end(); ++itSp2)
      {
        long ix = (*itSp2)[0];
        long iy = (*itSp2)[1];
        long iz = (*itSp2)[2];
        
        m_lp2->push_back(*itSp2);

        mp_label->set(ix, iy, iz, 2);
      }


//     // debug
//     labelsCoherentCheck1();

  }

  /*================================================================================
    initializeLabel*/
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::initializeLabel()
  {
    if (m_nx + m_ny + m_nz == 0)
      {
	std::cerr<<"set mp_img first.\n";
	raise(SIGABRT);
      }

    //find interface and mark as 0, create Lz
    unsigned char defaultLabel = 0;
    mp_label.reset(new cArray3D< char >(m_nx, m_ny, m_nz, defaultLabel) );

    return;
  }


  /*================================================================================
    initializePhi*/
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::initializePhi()
  {
    if (m_nx + m_ny + m_nz == 0)
      {
	std::cerr<<"set mp_img first.\n";
	raise(SIGABRT);
      }

    double arbitraryInitPhi = 1000;
    mp_phi.reset(new cArray3D< double >(m_nx, m_ny, m_nz, arbitraryInitPhi) );

    return;
  }


  /* ============================================================
     initializeSFLSFromMask    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::initializeSFLSFromMask()
  {
    if (!mp_mask)
      {
	std::cerr<<"set mp_mask first.\n";
	raise(SIGABRT);
      }


    initializePhi();
    initializeLabel();

    for (long ix = 0; ix < m_nx; ++ix) 
      {
        for (long iy = 0; iy < m_ny; ++iy) 
          {
	    for (long iz = 0; iz < m_nz; ++iz) 
	      {
		//mark the inside and outside of label and phi
		if( mp_mask->get(ix, iy, iz) == 0 )
		  { 
		    mp_label->set(ix, iy, iz, 3); 
		    mp_phi->set(ix, iy, iz, 3); 
		  }
		else
		  { 
		    mp_label->set(ix, iy, iz, -3); 
		    mp_phi->set(ix, iy, iz, -3); 

		    if ( (iy+1 < m_ny && mp_mask->get(ix, iy+1, iz) == 0)	\
			 || (iy-1 >= 0 && mp_mask->get(ix, iy-1, iz) == 0)	\
			 || (ix+1 < m_nx && mp_mask->get(ix+1, iy, iz) == 0) \
			 || (ix-1 >= 0 && mp_mask->get(ix-1, iy, iz) == 0)	\
			 || (iz+1 < m_nz && mp_mask->get(ix, iy, iz+1) == 0) \
			 || (iz-1 >= 0 && mp_mask->get(ix, iy, iz-1) == 0) )
		      {
			m_lz->push_back(CSFLSLayer::NodeType(ix, iy, iz));

			mp_label->set(ix, iy, iz, 0);
			mp_phi->set(ix, iy, iz, 0.0);
		      }
		  }
	      }
	  }
      }


    //scan Lz to create Ln1 and Lp1
    for (CSFLSLayer::const_iterator it = m_lz->begin(); it != m_lz->end(); ++it)
      {
        long ix = (*it)[0];
        long iy = (*it)[1];
        long iz = (*it)[2];
    
        if(ix+1 < m_nx)
          {
            if ( mp_label->get(ix+1, iy, iz) == 3 )
              {
                mp_label->set(ix+1, iy, iz, 1); 
                mp_phi->set(ix+1, iy, iz, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
              }
            else if ( mp_label->get(ix+1, iy, iz) == -3 )
              {
                mp_label->set(ix+1, iy, iz, -1); 
                mp_phi->set(ix+1, iy, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
              }
          }

        if(ix-1 >= 0)
          {
            if ( mp_label->get(ix-1, iy, iz) == 3 )
              {
                mp_label->set(ix-1, iy, iz, 1); 
                mp_phi->set(ix-1, iy, iz, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
              }
            else if ( mp_label->get(ix-1, iy, iz) == -3 )
              {
                mp_label->set(ix-1, iy, iz, -1); 
                mp_phi->set(ix-1, iy, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
              }
          }

        if(iy+1 < m_ny)
          {
            if ( mp_label->get(ix, iy+1, iz) == 3 )
              {
                mp_label->set(ix, iy+1, iz, 1); 
                mp_phi->set(ix, iy+1, iz, 1); 
            
                m_lp1->push_back(CSFLSLayer::NodeType(ix, iy+1, iz) );
              }
            else if ( mp_label->get(ix, iy+1, iz) == -3 )
              {
                mp_label->set(ix, iy+1, iz, -1); 
                mp_phi->set(ix, iy+1, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy+1, iz) );
              }
          }

        if(iy-1 >= 0)
          {
            if ( mp_label->get(ix, iy-1, iz) == 3 )
              {
                mp_label->set(ix, iy-1, iz, 1); 
                mp_phi->set(ix, iy-1, iz, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
              }
            else if ( mp_label->get(ix, iy-1, iz) == -3 )
              {
                mp_label->set(ix, iy-1, iz, -1); 
                mp_phi->set(ix, iy-1, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
              }
          }

        if(iz+1 < m_nz)
          {
            if ( mp_label->get(ix, iy, iz+1) == 3 )
              {
                mp_label->set(ix, iy, iz+1, 1); 
                mp_phi->set(ix, iy, iz+1, 1); 
            
                m_lp1->push_back(CSFLSLayer::NodeType(ix, iy, iz+1) );
              }
            else if ( mp_label->get(ix, iy, iz+1) == -3 )
              {
                mp_label->set(ix, iy, iz+1, -1); 
                mp_phi->set(ix, iy, iz+1, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy, iz+1) );
              }
          }

        if(iz-1 >= 0)
          {
            if ( mp_label->get(ix, iy, iz-1) == 3 )
              {
                mp_label->set(ix, iy, iz-1, 1); 
                mp_phi->set(ix, iy, iz-1, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
              }
            else if ( mp_label->get(ix, iy, iz-1) == -3 )
              {
                mp_label->set(ix, iy, iz-1, -1); 
                mp_phi->set(ix, iy, iz-1, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
              }
          }
      }


    //scan Ln1 to create Ln2
    for (CSFLSLayer::const_iterator it = m_ln1->begin(); it != m_ln1->end(); ++it)
      {
        long ix = (*it)[0];
        long iy = (*it)[1];
        long iz = (*it)[2];

        
        if(ix+1 < m_nx && mp_label->get(ix+1, iy, iz) == -3 )
          {
            mp_label->set(ix+1, iy, iz, -2); 
            mp_phi->set(ix+1, iy, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
          }

        if(ix-1 >= 0 && mp_label->get(ix-1, iy, iz) == -3 )
          {
            mp_label->set(ix-1, iy, iz, -2); 
            mp_phi->set(ix-1, iy, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
          }

        if(iy+1 < m_ny && mp_label->get(ix, iy+1, iz) == -3 )
          {
            mp_label->set(ix, iy+1, iz, -2); 
            mp_phi->set(ix, iy+1, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy+1, iz) );
          }

        if(iy-1 >= 0 && mp_label->get(ix, iy-1, iz) == -3 )
          {
            mp_label->set(ix, iy-1, iz, -2); 
            mp_phi->set(ix, iy-1, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
          }

        if(iz+1 < m_nz && mp_label->get(ix, iy, iz+1) == -3 )
          {
            mp_label->set(ix, iy, iz+1, -2); 
            mp_phi->set(ix, iy, iz+1, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy, iz+1) );
          }

        if(iz-1 >= 0 && mp_label->get(ix, iy, iz-1) == -3 )
          {
            mp_label->set(ix, iy, iz-1, -2); 
            mp_phi->set(ix, iy, iz-1, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
          }
      }

    //scan Lp1 to create Lp2
    for (CSFLSLayer::const_iterator it = m_lp1->begin(); it != m_lp1->end(); ++it)
      {
        long ix = (*it)[0];
        long iy = (*it)[1];
        long iz = (*it)[2];

        
        if(ix+1 < m_nx && mp_label->get(ix+1, iy, iz) == 3 )
          {
            mp_label->set(ix+1, iy, iz, 2); 
            mp_phi->set(ix+1, iy, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
          }

        if(ix-1 >= 0 && mp_label->get(ix-1, iy, iz) == 3 )
          {
            mp_label->set(ix-1, iy, iz, 2); 
            mp_phi->set(ix-1, iy, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
          }

        if(iy+1 < m_ny && mp_label->get(ix, iy+1, iz) == 3 )
          {
            mp_label->set(ix, iy+1, iz, 2); 
            mp_phi->set(ix, iy+1, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy+1, iz) );
          }

        if(iy-1 >= 0 && mp_label->get(ix, iy-1, iz) == 3 )
          {
            mp_label->set(ix, iy-1, iz, 2); 
            mp_phi->set(ix, iy-1, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
          }

        if(iz+1 < m_nz && mp_label->get(ix, iy, iz+1) == 3 )
          {
            mp_label->set(ix, iy, iz+1, 2); 
            mp_phi->set(ix, iy, iz+1, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy, iz+1) );
          }

        if(iz-1 >= 0 && mp_label->get(ix, iy, iz-1) == 3 )
          {
            mp_label->set(ix, iy, iz-1, 2); 
            mp_phi->set(ix, iy, iz-1, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
          }
      }
  }


  /* ============================================================
     getSFLSFromPhi    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::getSFLSFromPhi()
  {
    initializePhi();
    initializeLabel();

    for (long ix = 0; ix < m_nx; ++ix) 
      {
        for (long iy = 0; iy < m_ny; ++iy) 
          {
	    for (long iz = 0; iz < m_nz; ++iz) 
	      {
		//mark the inside and outside of label and phi
		if( mp_mask->get(ix, iy, iz) == 0 )
		  { 
		    mp_label->set(ix, iy, iz, 3); 
		    mp_phi->set(ix, iy, iz, 3); 
		  }
		else
		  { 
		    mp_label->set(ix, iy, iz, -3); 
		    mp_phi->set(ix, iy, iz, -3); 

		    if ( (iy+1 < m_ny && mp_mask->get(ix, iy+1, iz) == 0)	\
			 || (iy-1 >= 0 && mp_mask->get(ix, iy-1, iz) == 0)	\
			 || (ix+1 < m_nx && mp_mask->get(ix+1, iy, iz) == 0) \
			 || (ix-1 >= 0 && mp_mask->get(ix-1, iy, iz) == 0)	\
			 || (iz+1 < m_nz && mp_mask->get(ix, iy, iz+1) == 0) \
			 || (iz-1 >= 0 && mp_mask->get(ix, iy, iz-1) == 0) )
		      {
			m_lz->push_back(CSFLSLayer::NodeType(ix, iy, iz));

			mp_label->set(ix, iy, iz, 0);
			mp_phi->set(ix, iy, iz, 0.0);
		      }
		  }
	      }
	  }
      }


    //scan Lz to create Ln1 and Lp1
    for (CSFLSLayer::const_iterator it = m_lz->begin(); it != m_lz->end(); ++it)
      {
        long ix = (*it)[0];
        long iy = (*it)[1];
        long iz = (*it)[2];
    
        if(ix+1 < m_nx)
          {
            if ( mp_label->get(ix+1, iy, iz) == 3 )
              {
                mp_label->set(ix+1, iy, iz, 1); 
                mp_phi->set(ix+1, iy, iz, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
              }
            else if ( mp_label->get(ix+1, iy, iz) == -3 )
              {
                mp_label->set(ix+1, iy, iz, -1); 
                mp_phi->set(ix+1, iy, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
              }
          }

        if(ix-1 >= 0)
          {
            if ( mp_label->get(ix-1, iy, iz) == 3 )
              {
                mp_label->set(ix-1, iy, iz, 1); 
                mp_phi->set(ix-1, iy, iz, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
              }
            else if ( mp_label->get(ix-1, iy, iz) == -3 )
              {
                mp_label->set(ix-1, iy, iz, -1); 
                mp_phi->set(ix-1, iy, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
              }
          }

        if(iy+1 < m_ny)
          {
            if ( mp_label->get(ix, iy+1, iz) == 3 )
              {
                mp_label->set(ix, iy+1, iz, 1); 
                mp_phi->set(ix, iy+1, iz, 1); 
            
                m_lp1->push_back(CSFLSLayer::NodeType(ix, iy+1, iz) );
              }
            else if ( mp_label->get(ix, iy+1, iz) == -3 )
              {
                mp_label->set(ix, iy+1, iz, -1); 
                mp_phi->set(ix, iy+1, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy+1, iz) );
              }
          }

        if(iy-1 >= 0)
          {
            if ( mp_label->get(ix, iy-1, iz) == 3 )
              {
                mp_label->set(ix, iy-1, iz, 1); 
                mp_phi->set(ix, iy-1, iz, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
              }
            else if ( mp_label->get(ix, iy-1, iz) == -3 )
              {
                mp_label->set(ix, iy-1, iz, -1); 
                mp_phi->set(ix, iy-1, iz, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
              }
          }

        if(iz+1 < m_nz)
          {
            if ( mp_label->get(ix, iy, iz+1) == 3 )
              {
                mp_label->set(ix, iy, iz+1, 1); 
                mp_phi->set(ix, iy, iz+1, 1); 
            
                m_lp1->push_back(CSFLSLayer::NodeType(ix, iy, iz+1) );
              }
            else if ( mp_label->get(ix, iy, iz+1) == -3 )
              {
                mp_label->set(ix, iy, iz+1, -1); 
                mp_phi->set(ix, iy, iz+1, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy, iz+1) );
              }
          }

        if(iz-1 >= 0)
          {
            if ( mp_label->get(ix, iy, iz-1) == 3 )
              {
                mp_label->set(ix, iy, iz-1, 1); 
                mp_phi->set(ix, iy, iz-1, 1); 
            
                m_lp1->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
              }
            else if ( mp_label->get(ix, iy, iz-1) == -3 )
              {
                mp_label->set(ix, iy, iz-1, -1); 
                mp_phi->set(ix, iy, iz-1, -1); 
            
                m_ln1->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
              }
          }
      }


    //scan Ln1 to create Ln2
    for (CSFLSLayer::const_iterator it = m_ln1->begin(); it != m_ln1->end(); ++it)
      {
        long ix = (*it)[0];
        long iy = (*it)[1];
        long iz = (*it)[2];

        
        if(ix+1 < m_nx && mp_label->get(ix+1, iy, iz) == -3 )
          {
            mp_label->set(ix+1, iy, iz, -2); 
            mp_phi->set(ix+1, iy, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
          }

        if(ix-1 >= 0 && mp_label->get(ix-1, iy, iz) == -3 )
          {
            mp_label->set(ix-1, iy, iz, -2); 
            mp_phi->set(ix-1, iy, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
          }

        if(iy+1 < m_ny && mp_label->get(ix, iy+1, iz) == -3 )
          {
            mp_label->set(ix, iy+1, iz, -2); 
            mp_phi->set(ix, iy+1, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy+1, iz) );
          }

        if(iy-1 >= 0 && mp_label->get(ix, iy-1, iz) == -3 )
          {
            mp_label->set(ix, iy-1, iz, -2); 
            mp_phi->set(ix, iy-1, iz, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
          }

        if(iz+1 < m_nz && mp_label->get(ix, iy, iz+1) == -3 )
          {
            mp_label->set(ix, iy, iz+1, -2); 
            mp_phi->set(ix, iy, iz+1, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy, iz+1) );
          }

        if(iz-1 >= 0 && mp_label->get(ix, iy, iz-1) == -3 )
          {
            mp_label->set(ix, iy, iz-1, -2); 
            mp_phi->set(ix, iy, iz-1, -2); 
            
            m_ln2->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
          }
      }

    //scan Lp1 to create Lp2
    for (CSFLSLayer::const_iterator it = m_lp1->begin(); it != m_lp1->end(); ++it)
      {
        long ix = (*it)[0];
        long iy = (*it)[1];
        long iz = (*it)[2];

        
        if(ix+1 < m_nx && mp_label->get(ix+1, iy, iz) == 3 )
          {
            mp_label->set(ix+1, iy, iz, 2); 
            mp_phi->set(ix+1, iy, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix+1, iy, iz) );
          }

        if(ix-1 >= 0 && mp_label->get(ix-1, iy, iz) == 3 )
          {
            mp_label->set(ix-1, iy, iz, 2); 
            mp_phi->set(ix-1, iy, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix-1, iy, iz) );
          }

        if(iy+1 < m_ny && mp_label->get(ix, iy+1, iz) == 3 )
          {
            mp_label->set(ix, iy+1, iz, 2); 
            mp_phi->set(ix, iy+1, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy+1, iz) );
          }

        if(iy-1 >= 0 && mp_label->get(ix, iy-1, iz) == 3 )
          {
            mp_label->set(ix, iy-1, iz, 2); 
            mp_phi->set(ix, iy-1, iz, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy-1, iz) );
          }

        if(iz+1 < m_nz && mp_label->get(ix, iy, iz+1) == 3 )
          {
            mp_label->set(ix, iy, iz+1, 2); 
            mp_phi->set(ix, iy, iz+1, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy, iz+1) );
          }

        if(iz-1 >= 0 && mp_label->get(ix, iy, iz-1) == 3 )
          {
            mp_label->set(ix, iy, iz-1, 2); 
            mp_phi->set(ix, iy, iz-1, 2); 
            
            m_lp2->push_back( CSFLSLayer::NodeType(ix, iy, iz-1) );
          }
      }
  }


  /* ============================================================
     doSegmenation    */
  template< typename TPixel >
  void
  CSFLSSegmentor3D< TPixel >
  ::doSegmenation()
  {
    double arbitraryInitPhi = 1000;
    mp_phi.reset(new cArray3D< double >(m_nx, m_ny, m_nz, arbitraryInitPhi) );

    // douher::saveAsImage3< double >(mp_phi, "init0.nrrd");


    /*============================================================
     * From the initial mask, generate: 1. SFLS, 2. mp_label and
     * 3. mp_phi.      
     */
    initializeSFLS();

    //douher::saveAsImage3< double >(mp_phi, "initPhi.nrrd");

    for (unsigned int it = 0; it < m_numIter; ++it)
      {
	std::cout<<"iteration "<<it<<"\n"<<std::flush;

        /*--------------------------------------------------
          Compute the force on the zero level set, NOT on the whole domain.
          This is NOT implemented in this base class.    

          This function will compute the m_force. m_force has the same
          size as the m_ln, indicating the change at each pixel on the
          zero level set.
        */
        computeForce(); 


        normalizeForce();

//         // debug
//         for (std::list< double >::const_iterator itf = this->m_force.begin(); itf != this->m_force.end(); ++itf)
//           {
//             std::cout<<(*itf)<<", ";
//           }
//         std::cout<<std::endl<<it<<std::endl<<std::endl;



//         //debug//
//         labelsCoherentCheck1();

        oneStepLevelSetEvolution();


//         //debug//
//         std::cout<<"-----------------------"<<it<<"---------------------------"<<std::endl;
//         std::cout<<"lz \t ln1 \t ln2 \t lp1 \t lp2 \n";
//         std::cout<<m_lz->size()<<"\t"<<m_ln1->size()<<"\t"<<m_ln2->size()<<"\t"<<m_lp1->size()<<"\t"<<m_lp2->size()<<std::endl;
//         std::cout<<"--------------------------------------------------"<<std::endl;


//         // debug
//         labelsCoherentCheck1();


        //        douher::saveAsImage3< double >(mp_phi, "temp.nrrd");
      }
  }



  /*============================================================
    computeKappa 

    Compute kappa at a point in the zero level set  */
  template< typename TPixel >
  double
  CSFLSSegmentor3D< TPixel >
  ::computeKappa(long ix, long iy, long iz)
  {
    double kappa = 0;

    double dx = 0;
    double dy = 0;
    double dz = 0;

    double dxx = 0;
    double dyy = 0;
    double dzz = 0;

    double dx2 = 0;
    double dy2 = 0;
    double dz2 = 0;

    double dxy = 0;
    double dxz = 0;
    double dyz = 0;

    char xok = 0;
    char yok = 0;
    char zok = 0;


    if( ix+1 < this->m_nx && ix-1 >=0 )
      { 
        xok = 1;
      }

    if( iy+1 < this->m_ny && iy-1 >=0 ) 
      {
        yok = 1;
      }

    if( iz+1 < this->m_nz && iz-1 >=0 ) 
      {
        zok = 1;
      }

    if (xok)
      {
        dx  = (this->mp_phi->get(ix+1, iy, iz) - this->mp_phi->get(ix-1, iy, iz) )/2.0;
        dxx = this->mp_phi->get(ix+1, iy, iz) - 2.0*(this->mp_phi->get(ix, iy, iz)) + this->mp_phi->get(ix-1, iy, iz);
        dx2 = dx*dx;
      }

    if (yok)
      {
        dy  = (this->mp_phi->get(ix, iy+1, iz) - this->mp_phi->get(ix, iy-1, iz) )/2.0;
        dyy = this->mp_phi->get(ix, iy+1, iz) - 2*(this->mp_phi->get(ix, iy, iz)) + this->mp_phi->get(ix, iy-1, iz);
        dy2 = dy*dy;
      }

    if (zok)
      {
        dz  = (this->mp_phi->get(ix, iy, iz+1) - this->mp_phi->get(ix, iy, iz-1) )/2.0;
        dzz = this->mp_phi->get(ix, iy, iz+1) - 2.0*(this->mp_phi->get(ix, iy, iz)) + this->mp_phi->get(ix, iy, iz-1);
        dz2 = dz*dz;
      }


    if(xok && yok)
      {
        dxy = 0.25*(this->mp_phi->get(ix+1, iy+1, iz) + this->mp_phi->get(ix-1, iy-1, iz) \
                    - this->mp_phi->get(ix+1, iy-1, iz) - this->mp_phi->get(ix-1, iy+1, iz));
      }

    if(xok && zok)
      {
        dxz = 0.25*(this->mp_phi->get(ix+1, iy, iz+1) + this->mp_phi->get(ix-1, iy, iz-1)\
                    - this->mp_phi->get(ix+1, iy, iz-1) - this->mp_phi->get(ix-1, iy, iz+1));
      }

    if(yok && zok)
      {
        dyz = 0.25*(this->mp_phi->get(ix, iy+1, iz+1) + this->mp_phi->get(ix, iy-1, iz-1)\
                    - this->mp_phi->get(ix, iy+1, iz-1) - this->mp_phi->get(ix, iy-1, iz+1));
      }

    kappa = (dxx*(dy2 + dz2) + dyy*(dx2 + dz2) + dzz*(dx2 + dy2) - 2*dx*dy*dxy - 2*dx*dz*dxz - 2*dy*dz*dyz)/(dx2 + dy2 + dz2 + 1e-10);

    return kappa;
  }

} // douher

#endif
