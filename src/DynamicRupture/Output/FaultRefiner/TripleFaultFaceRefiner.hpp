#ifndef SEISSOL_TRIPLEFAULTFACEREFINER_HPP
#define SEISSOL_TRIPLEFAULTFACEREFINER_HPP

#include "FaultRefinerInterface.hpp"

namespace seissol {
  namespace dr {
    namespace output {
      class TripleFaultFaceRefiner;
    }
  }
}

class seissol::dr::output::TripleFaultFaceRefiner : public seissol::dr::output::FaultRefinerInterface {
  void refine(const int LocalSideId) override {};
  int getNumSubTriangles() {return 3;}
};


#endif //SEISSOL_TRIPLEFAULTFACEREFINER_HPP
