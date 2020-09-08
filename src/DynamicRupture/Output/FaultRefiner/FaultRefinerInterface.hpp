#ifndef SEISSOL_INTERFACE_HPP
#define SEISSOL_INTERFACE_HPP

#include "DynamicRupture/Output/DrDataTypes.hpp"
#include <memory>

namespace seissol {
  namespace dr {
    namespace output {
      class FaultRefinerInterface;
    }
  }
}

class seissol::dr::output::FaultRefinerInterface {
public:
  virtual void refine(const int LocalSideId) = 0;
  virtual int getNumSubTriangles() = 0;
};

#endif //SEISSOL_INTERFACE_HPP
