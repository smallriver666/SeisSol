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
  int getNumSubTriangles() override {return 3;}
  void refineAndAccumulate(int RefinementLevel,
                           int FaultFaceIndex,
                           int LocalFaceSideId,
                           ExtTriangle ReferenceFace,
                           ExtTriangle GlobalFace) override {

    if (RefinementLevel == 0) {
      ReceiverPointT Receiver{};
      Receiver.IsInside = true;
      Receiver.FaultFaceIndex = FaultFaceIndex;
      Receiver.LocalFaceSideId = LocalFaceSideId;
      Receiver.GlobalReceiverIndex = m_Points.size();
      Receiver.Global = getMidTrianglePoint(GlobalFace);
      Receiver.Referece = getMidTrianglePoint(ReferenceFace);
      Receiver.GlobalSubTet = GlobalFace;

      m_Points.push_back(Receiver);
      return;
    }

    ExtVrtxCoords RefMidTriangle = getMidTrianglePoint(ReferenceFace);
    ExtVrtxCoords GlbMidTriangle = getMidTrianglePoint(GlobalFace);

    {
      // First sub-face (-triangle)
      ExtTriangle SubReferenceFace(ReferenceFace.p1, ReferenceFace.p2, RefMidTriangle);
      ExtTriangle SubGlobalFace(GlobalFace.p1, GlobalFace.p2, GlbMidTriangle);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, LocalFaceSideId, SubReferenceFace, SubGlobalFace);
    }

    {
      // Second sub-face (-triangle)
      ExtTriangle SubReferenceFace(RefMidTriangle, ReferenceFace.p2, ReferenceFace.p3);
      ExtTriangle SubGlobalFace(GlbMidTriangle, ReferenceFace.p2, ReferenceFace.p3);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, LocalFaceSideId, SubReferenceFace, SubGlobalFace);
    }

    {
      // Third sub-face (-triangle)
      ExtTriangle SubReferenceFace(ReferenceFace.p1, RefMidTriangle, ReferenceFace.p3);
      ExtTriangle SubGlobalFace(GlobalFace.p1, GlbMidTriangle, GlobalFace.p3);
      this->refineAndAccumulate(RefinementLevel - 1, FaultFaceIndex, LocalFaceSideId, SubReferenceFace, SubGlobalFace);
    }
  }
};


#endif //SEISSOL_TRIPLEFAULTFACEREFINER_HPP
