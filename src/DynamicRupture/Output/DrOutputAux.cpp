#include "DrOutputAux.hpp"
#include "DynamicRupture/DR_math.h"
#include "Geometry/MeshTools.h"
#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/Transformation.h"
#include "Numerical_aux/BasisFunction.h"
#include <Eigen/Dense>
#include <limits>

namespace seissol {
  namespace dr {

    int getElementVertexId(int LocalSideId, int LocalFaceVertexId) {
      return MeshTools::FACE2NODES[LocalSideId][LocalFaceVertexId];
    }


    ExtTriangle getReferenceFace(int LocalSideId) {
      ExtTriangle ReferenceFace;
      switch (LocalSideId) {
        case 0:
          ReferenceFace.p1.xi = 0.0; ReferenceFace.p1.eta = 0.0; ReferenceFace.p1.zeta = 0.0;
          ReferenceFace.p2.xi = 0.0; ReferenceFace.p2.eta = 1.0; ReferenceFace.p2.zeta = 0.0;
          ReferenceFace.p3.xi = 1.0; ReferenceFace.p3.eta = 0.0; ReferenceFace.p3.zeta = 0.0;
          break;
        case 1:
          ReferenceFace.p1.xi = 0.0; ReferenceFace.p1.eta = 0.0; ReferenceFace.p1.zeta = 0.0;
          ReferenceFace.p2.xi = 1.0; ReferenceFace.p2.eta = 0.0; ReferenceFace.p2.zeta = 0.0;
          ReferenceFace.p3.xi = 0.0; ReferenceFace.p3.eta = 0.0; ReferenceFace.p3.zeta = 1.0;
          break;
        case 2:
          ReferenceFace.p1.xi = 0.0; ReferenceFace.p1.eta = 0.0; ReferenceFace.p1.zeta = 0.0;
          ReferenceFace.p2.xi = 0.0; ReferenceFace.p2.eta = 0.0; ReferenceFace.p2.zeta = 1.0;
          ReferenceFace.p3.xi = 0.0; ReferenceFace.p3.eta = 1.0; ReferenceFace.p3.zeta = 0.0;
          break;
        case 3:
          ReferenceFace.p1.xi = 1.0; ReferenceFace.p1.eta = 0.0; ReferenceFace.p1.zeta = 0.0;
          ReferenceFace.p2.xi = 0.0; ReferenceFace.p2.eta = 1.0; ReferenceFace.p2.zeta = 0.0;
          ReferenceFace.p3.xi = 0.0; ReferenceFace.p3.eta = 0.0; ReferenceFace.p3.zeta = 1.0;
          break;
        default:
          throw std::runtime_error("Unknown Local Side Id. Must be 0, 1, 2 or 3");
      }

      return ReferenceFace;
    }


    void computeStrikeAndDipVectors(const VrtxCoords Normal, VrtxCoords Strike, VrtxCoords Dip) {
      // Note: equations are explained in documentation -> left-lateral-right-lateral-normal-reverse

      // compute normalized strike vector
      auto StrikeInvLength = 1.0 / std::sqrt(Normal[0] * Normal[0] + Normal[1] * Normal[1]);
      Strike[0] = Normal[1] * StrikeInvLength;
      Strike[1] = -Normal[0] * StrikeInvLength;
      Strike[2] = 0.0;

      // compute normalized dip vector
      Dip[0] = -1.0 * Strike[1] * Normal[2];
      Dip[1] = Strike[0] * Normal[2];
      Dip[2] = Strike[1] * Normal[0] - Strike[0] * Normal[1];
      auto DipInvLength = 1.0 / std::sqrt(Dip[0] * Dip[0] + Dip[1] * Dip[1] + Dip[2] * Dip[2]);
      Dip[0] *= DipInvLength;
      Dip[1] *= DipInvLength;
      Dip[2] *= DipInvLength;
    }


    ExtVrtxCoords getMidTrianglePoint(const ExtTriangle& Triangle) {
      ExtVrtxCoords AvgPoint{};
      for (int Axis = 0; Axis < 3; ++Axis) {
        AvgPoint.Coords[Axis] = (Triangle.p1.Coords[Axis] + Triangle.p2.Coords[Axis] + Triangle.p3.Coords[Axis]) / 3.0;
      }
      return AvgPoint;
    }


    ExtVrtxCoords getMidPoint(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2) {
      ExtVrtxCoords MidPoint{};
      for (int Axis = 0; Axis < 3; ++Axis) {
        MidPoint.Coords[Axis] = 0.5 * (p1.Coords[Axis] + p2.Coords[Axis]);
      }
      return MidPoint;
    }


    std::tuple<unsigned, std::shared_ptr<double []>, std::shared_ptr<double []>>
    generateTriangleQuadrature(unsigned PolyDegree) {

      // allocate data
      unsigned numQuadraturePoints = PolyDegree * PolyDegree;
      std::shared_ptr<double []> Weights( new double[numQuadraturePoints], std::default_delete<double[]>());
      std::shared_ptr<double []> Points( new double[2 * numQuadraturePoints], std::default_delete<double[]>());

      // Generate triangle quadrature points and weights (Factory Method)
      seissol::quadrature::TriangleQuadrature(reshape<2>(&Points[0]), &Weights[0], PolyDegree);
      return std::make_tuple(numQuadraturePoints, Weights, Points);
    }


    double distance(const double V1[2], const double V2[2]) {
      Eigen::Vector2d Vector1(V1[0],V1[1]), Vector2(V2[0],V2[1]);
      return (Vector1 - Vector2).norm();
    }


    std::pair<int, double> getNearestFacePoint(const double TargetPoint[2],
                                               const double (*FacePoints)[2],
                                               const unsigned NumFacePoints) {

      int NearestPoint{-1};
      double ShortestDistance = std::numeric_limits<double>::max();

      for (unsigned Index = 0; Index < NumFacePoints; ++Index) {
        double NextPoint[2] = {FacePoints[Index][0], FacePoints[Index][1]};

        auto CurrentDistance = distance(TargetPoint, NextPoint);
        if (ShortestDistance > CurrentDistance) {
          ShortestDistance = CurrentDistance;
          NearestPoint = Index;
        }
      }
      return std::make_pair(NearestPoint, ShortestDistance);
    }


    void assignNearestGaussianPoints(ReceiverPointsT& GeoPoints) {
      std::shared_ptr<double []> Weights = nullptr;
      std::shared_ptr<double []> PointsData = nullptr;
      unsigned NumPoints{};

      std::tie(NumPoints, Weights, PointsData) = generateTriangleQuadrature(CONVERGENCE_ORDER + 1);
      double (*TrianglePoints2D)[2] = reshape<2>(&PointsData[0]);

      for (auto& GeoPoint: GeoPoints) {

        double TargetPoint2D[2];
        transformations::XiEtaZeta2chiTau(GeoPoint.LocalFaceSideId, GeoPoint.Referece.Coords, TargetPoint2D);

        int NearestPoint{-1};
        double ShortestDistance = std::numeric_limits<double>::max();
        std::tie(NearestPoint, ShortestDistance) = getNearestFacePoint(TargetPoint2D,
                                                                       TrianglePoints2D,
                                                                       NumPoints);
        GeoPoint.NearestGpIndex = NearestPoint;
        GeoPoint.DistanceToNearestGp = ShortestDistance;
      }
    }


    void projectPointToFace(ExtVrtxCoords& Point, const ExtTriangle& Face, const VrtxCoords FaceNormal) {
      using namespace Eigen;

      Vector3d OriginalPoint(Point.x, Point.y, Point.z);

      Vector3d R = OriginalPoint - Vector3d(Face.p1.x, Face.p1.y, Face.p1.z);
      Vector3d Direction(FaceNormal[0], FaceNormal[1], FaceNormal[2]);
      Direction.normalize();

      auto Displacement = -Direction.dot(R);
      Vector3d TargetPoint = OriginalPoint + Displacement * Direction;

      for (int i = 0; i < 3; ++i) {
        Point.Coords[i] = TargetPoint(i);
      }
    }

    PlusMinusBasisFunctionsT getPlusMinusBasisFunctions(const VrtxCoords PointCoords,
                                                        const VrtxCoords PlusElementCoords[4],
                                                        const VrtxCoords MinusElementCoords[4]) {

      PlusMinusBasisFunctionsT BasisFunctions{};
      Eigen::Vector3d Point(PointCoords[0], PointCoords[1], PointCoords[2]);

      {
        auto PlusXiEtaZeta = transformations::tetrahedronGlobalToReference(PlusElementCoords[0],
                                                                           PlusElementCoords[1],
                                                                           PlusElementCoords[2],
                                                                           PlusElementCoords[3],
                                                                           Point);

        basisFunction::SampledBasisFunctions<real> Sampler(CONVERGENCE_ORDER,
                                                           PlusXiEtaZeta[0],
                                                           PlusXiEtaZeta[1],
                                                           PlusXiEtaZeta[2]);
        BasisFunctions.PlusSide = std::move(Sampler.m_data);
      }

      std::vector<real> MinusBasisFunctions;
      {
        auto MinusXiEtaZeta = transformations::tetrahedronGlobalToReference(MinusElementCoords[0],
                                                                            MinusElementCoords[1],
                                                                            MinusElementCoords[2],
                                                                            MinusElementCoords[3],
                                                                            Point);

        basisFunction::SampledBasisFunctions<real> Sampler(CONVERGENCE_ORDER,
                                                           MinusXiEtaZeta[0],
                                                           MinusXiEtaZeta[1],
                                                           MinusXiEtaZeta[2]);
        BasisFunctions.MinusSide = std::move(Sampler.m_data);
      }

      return BasisFunctions;
    }
  }
}