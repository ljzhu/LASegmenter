#ifndef SFLSLayer_h_
#define SFLSLayer_h_


//douher
#include "basicHeaders.h"
#include "cArrayOp.h"


// itk
#include "vnl/vnl_vector_fixed.h"

// std
#include <list>

namespace douher
{
  class CSFLSLayer : public std::list< vnl_vector_fixed< long, 3 > >
  {
  public:
    typedef CSFLSLayer Self;
    typedef boost::shared_ptr< Self > Pointer;

    typedef std::list< vnl_vector_fixed< long, 3 > > SuperClassType;


    typedef vnl_vector_fixed< long, 3 > NodeType;

    // ctor 
    CSFLSLayer() : SuperClassType() {}


    // New
    static Pointer New() { return Pointer(new Self); }
  };
}


#endif
