#ifndef SEISSOL_QUADFAULTFACEREFINER_HPP
#define SEISSOL_QUADFAULTFACEREFINER_HPP

#include "FaultRefinerInterface.hpp"

namespace seissol {
  namespace dr {
    namespace output {
      class QuadFaultFaceRefiner;
    }
  }
}

class seissol::dr::output::QuadFaultFaceRefiner : public seissol::dr::output::FaultRefinerInterface {
  void refine(const int LocalSideId) override {};
  int getNumSubTriangles() {return 4;}
};


#endif //SEISSOL_QUADFAULTFACEREFINER_HPP
