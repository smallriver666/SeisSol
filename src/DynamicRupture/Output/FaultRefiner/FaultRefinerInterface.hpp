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
  virtual void refineAndAccumulate(int RefinementLevel,
                                   int FaultFaceIndex,
                                   int LocalFaceSideId,
                                   ExtTriangle ReferenceFace,
                                   ExtTriangle GlobalFace) = 0;

  ReceiverPointsT&& moveAllReceiverPoints() {return std::move(m_Points);}
  ReceiverPointsT getAllReceiverPoints() {return m_Points;}

protected:
  ReceiverPointsT m_Points{};
};

#endif //SEISSOL_INTERFACE_HPP
