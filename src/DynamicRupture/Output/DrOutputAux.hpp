#ifndef SEISSOL_DROUTPUTAUX_HPP
#define SEISSOL_DROUTPUTAUX_HPP

#include "DrDataTypes.hpp"
#include "Geometry/MeshReader.h"
#include <memory>

namespace seissol {

  template<int N, typename T>
  auto reshape(T* Ptr) -> T (*)[N] {
    return reinterpret_cast<T(*)[N]>(Ptr);
  }

  namespace dr {
    int getElementVertexId(int LocalSideId, int LocalFaceVertexId);

    ExtTriangle getReferenceFace(int LocalSideId);

    void computeStrikeAndDipVectors(const VrtxCoords Normal, VrtxCoords Strike, VrtxCoords Dip);

    ExtVrtxCoords getMidTrianglePoint(const ExtTriangle& Triangle);

    ExtVrtxCoords getMidPoint(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2);

    std::tuple<unsigned, std::shared_ptr<double []>, std::shared_ptr<double []>>
    generateTriangleQuadrature(unsigned PolyDegree);

    void assignNearestGaussianPoints(ReceiverPointsT& GeoPoints);

    std::pair<int, double> getNearestFacePoint(const double TargetPoint[2],
                                               const double (*FacePoints)[2],
                                               unsigned NumFacePoints);

    void projectPointToFace(ExtVrtxCoords& Point, const ExtTriangle& Face, const VrtxCoords FaceNormal);

    PlusMinusBasisFunctionsT getPlusMinusBasisFunctions(const VrtxCoords Point,
                                                        const VrtxCoords PlusElementCoords[4],
                                                        const VrtxCoords MinusElementCoords[4]);
  }
}

#endif //SEISSOL_DROUTPUTAUX_HPP
