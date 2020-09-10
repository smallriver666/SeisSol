#ifndef SEISSOL_QUADFAULTFACEREFINER_HPP
#define SEISSOL_QUADFAULTFACEREFINER_HPP

#include "FaultRefinerInterface.hpp"
#include "Numerical_aux/Transformation.h"

namespace seissol {
  namespace dr {
    namespace output {
      class QuadFaultFaceRefiner;
    }
  }
}

class seissol::dr::output::QuadFaultFaceRefiner : public seissol::dr::output::FaultRefinerInterface {
  int getNumSubTriangles() override {return 4;}
  void refineAndAccumulate(int RefinementLevel,
                           int FaultFaceIndex,
                           TriangleT ReferenceFace,
                           TriangleT GlobalFace) override {

    if (RefinementLevel == 0) {
      ReceiverPointT Receiver{};
      Receiver.IsInside = true;
      Receiver.FaultFaceIndex = FaultFaceIndex;
      Receiver.GlobalReceiverIndex = m_Points.size();
      Receiver.Global = getMidTrianglePoint(GlobalFace);
      Receiver.Referece = getMidTrianglePoint(ReferenceFace);
      Receiver.GlobalSubTet = GlobalFace;

      m_Points.push_back(Receiver);
      return;
    }

    GeoCoordsT RefMidPoint1 = getMidPoint(ReferenceFace.p1, ReferenceFace.p2);
    GeoCoordsT RefMidPoint2 = getMidPoint(ReferenceFace.p2, ReferenceFace.p3);
    GeoCoordsT RefMidPoint3 = getMidPoint(ReferenceFace.p3, ReferenceFace.p1);

    GeoCoordsT GlbMidPoint1 = getMidPoint(GlobalFace.p1, GlobalFace.p2);
    GeoCoordsT GlbMidPoint2 = getMidPoint(GlobalFace.p2, GlobalFace.p3);
    GeoCoordsT GlbMidPoint3 = getMidPoint(GlobalFace.p3, GlobalFace.p1);

    {
      // First sub-face (-triangle)
      TriangleT SubReferenceFace(ReferenceFace.p1, RefMidPoint1, RefMidPoint3);
      TriangleT SubGlobalFace(GlobalFace.p1, GlbMidPoint1, GlbMidPoint3);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, SubReferenceFace, SubGlobalFace);
    }

    {
      // Second sub-face (-triangle)
      TriangleT SubReferenceFace(RefMidPoint1, ReferenceFace.p2, RefMidPoint2);
      TriangleT SubGlobalFace(GlbMidPoint1, GlobalFace.p2, GlbMidPoint2);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, SubReferenceFace, SubGlobalFace);
    }

    {
      // Third sub-face (-triangle)
      TriangleT SubReferenceFace(RefMidPoint1, RefMidPoint2, RefMidPoint3);
      TriangleT SubGlobalFace(GlbMidPoint1, GlbMidPoint2, GlbMidPoint3);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, SubReferenceFace, SubGlobalFace);
    }

    {
      // Fourth sub-face (-triangle)
      TriangleT SubReferenceFace(RefMidPoint3, RefMidPoint2, ReferenceFace.p3);
      TriangleT SubGlobalFace(GlbMidPoint3, GlbMidPoint2, GlobalFace.p3);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, SubReferenceFace, SubGlobalFace);
    }
  }
};


#endif //SEISSOL_QUADFAULTFACEREFINER_HPP
