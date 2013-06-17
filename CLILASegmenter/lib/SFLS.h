#ifndef SFLS_h_
#define SFLS_h_




//douher
#include "SFLSLayer.h"


namespace douher
{
  class CSFLS
  {
  public:
    typedef CSFLS Self;
    typedef boost::shared_ptr< Self > Pointer;

    // ctor
    CSFLS()
    {
      m_lz.reset(new CSFLSLayer);
      m_ln1.reset(new CSFLSLayer);
      m_ln2.reset(new CSFLSLayer);
      m_lp1.reset(new CSFLSLayer);
      m_lp2.reset(new CSFLSLayer);
    }
  
    //   static PointerType New() 
    //   {
    //     return PointerType(new CSFLS);
    //   }

    CSFLSLayer::Pointer m_lz;
    CSFLSLayer::Pointer m_ln1;
    CSFLSLayer::Pointer m_ln2;
    CSFLSLayer::Pointer m_lp1;
    CSFLSLayer::Pointer m_lp2;
  };
}

//#include "SFLS.hpp"

#endif
