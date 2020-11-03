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
  virtual int getNumSubTriangles() = 0;
  virtual void refineAndAccumulate(int refinementLevel,
                                   int faultFaceIndex,
                                   int localFaceSideId,
                                   ExtTriangle referenceFace,
                                   ExtTriangle globalFace) = 0;

  ReceiverPointsT&& moveAllReceiverPoints() {return std::move(points);}
  ReceiverPointsT getAllReceiverPoints() {return points;}

protected:
  ReceiverPointsT points{};
};

#endif //SEISSOL_INTERFACE_HPP
