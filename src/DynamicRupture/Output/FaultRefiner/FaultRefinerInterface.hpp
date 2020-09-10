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
                                   TriangleT ReferenceFace,
                                   TriangleT GlobalFace) = 0;

  ReceiverPointsT&& moveAllReceiverPoints() {return std::move(m_Points);}
  ReceiverPointsT getAllReceiverPoints() {return m_Points;}

protected:

  static GeoCoordsT getMidTrianglePoint(const TriangleT& Triangle) {
    GeoCoordsT AvgPoint{};
    for (int Axis = 0; Axis < 3; ++Axis) {
      AvgPoint.Values[Axis] = (Triangle.p1.Values[Axis] + Triangle.p2.Values[Axis] + Triangle.p3.Values[Axis]) / 3.0;
    }
    return AvgPoint;
  };

  static GeoCoordsT getMidPoint(const GeoCoordsT& p1, const GeoCoordsT& p2) {
    GeoCoordsT MidPoint{};
    for (int Axis = 0; Axis < 3; ++Axis) {
      MidPoint.Values[Axis] = 0.5 * (p1.Values[Axis] + p2.Values[Axis]);
    }
    return MidPoint;
  };

  ReceiverPointsT m_Points{};
};

#endif //SEISSOL_INTERFACE_HPP
