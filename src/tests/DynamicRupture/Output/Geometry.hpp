#ifndef SEISSOL_GEOMETRY_HPP
#define SEISSOL_GEOMETRY_HPP

#include <cxxtest/TestSuite.h>
#include "DynamicRupture/Output/DrOutputAux.hpp"
#include "Numerical_aux/BasisFunction.h"
#include "Numerical_aux/Transformation.h"
#include "Initializer/PointMapper.h"
#include "Geometry/MeshReader.h"
#include "tests/Geometry/MockReader.h"
#include <iostream>
#include <Eigen/Dense>

namespace seissol {
  namespace unit_test {
    namespace dr {
      class Geometry;
    }
  }
}

using namespace seissol;
using namespace seissol::dr;

class seissol::unit_test::dr::Geometry : public CxxTest::TestSuite {
public:
  void testProjection() {

    // Given a reference triangle in the first octant
    // TargetPoint - intersection of a line (which starts from the origin and goes along [1,1,1] vector)
    // and the inclined face (4th face)
    ExtVrtxCoords TargetPoint{1.0/3.0, 1.0/3.0, 1.0/3.0};

    // 4th face
    ExtTriangle Face(ExtVrtxCoords{1.0, 0.0, 0.0},
                     ExtVrtxCoords{0.0, 1.0, 0.0},
                     ExtVrtxCoords{0.0, 0.0, 1.0});

    const double EPS = 1e-6;
    {
      ExtVrtxCoords TestPoint{0.0, 0.0, 0.0};
      VrtxCoords NormalDirection{1.0, 1.0, 1.0};

      projectPointToFace(TestPoint, Face, NormalDirection);


      TS_ASSERT_DELTA(TestPoint.x, TargetPoint.x, EPS);
      TS_ASSERT_DELTA(TestPoint.y, TargetPoint.y, EPS);
      TS_ASSERT_DELTA(TestPoint.z, TargetPoint.z, EPS);

    }
    {
      ExtVrtxCoords TestPoint{1.0, 1.0, 1.0};
      VrtxCoords NormalDirection{-1.0, -1.0, -1.0};

      projectPointToFace(TestPoint, Face, NormalDirection);

      TS_ASSERT_DELTA(TestPoint.x, TargetPoint.x, EPS);
      TS_ASSERT_DELTA(TestPoint.y, TargetPoint.y, EPS);
      TS_ASSERT_DELTA(TestPoint.z, TargetPoint.z, EPS);
    }
  }


  void testClosestPoint() {
    double TargetPoint[2] = {-0.25, -0.25};
    double FacePoints[4][2] = {{1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}};

    int TestPointId{-1};
    double TestDistance{-1.0};
    std::tie(TestPointId, TestDistance) = getNearestFacePoint(TargetPoint, FacePoints, 4);

    const  double EPS = 1e-6;
    TS_ASSERT_EQUALS(TestPointId, 2);
    TS_ASSERT_DELTA(TestDistance, std::sqrt(2 * 0.75 * 0.75), EPS);
  }


  void testMiddlePoint() {
    ExtVrtxCoords Point1{1.0, 2.0, 3.0};
    ExtVrtxCoords Point2{-3.0, -2.0, -1.0};

    auto TestMiddle = getMidPoint(Point1, Point2);

    const  double EPS = 1e-6;
    TS_ASSERT_DELTA(TestMiddle.x, -1.0, EPS);
    TS_ASSERT_DELTA(TestMiddle.y, 0.0, EPS);
    TS_ASSERT_DELTA(TestMiddle.z, 1.0, EPS);
  }


  void testMidTrianglePoint() {
    ExtVrtxCoords Point1{0.5, 0.0, 2.0};
    ExtVrtxCoords Point2{-0.5, 0.0, 3.0};
    ExtVrtxCoords Point3{3.0, 1.0, -2.0};
    ExtTriangle Triangle(Point1, Point2, Point3);

    auto TestMiddle = getMidTrianglePoint(Triangle);

    const  double EPS = 1e-6;
    TS_ASSERT_DELTA(TestMiddle.x, 1.0, EPS);
    TS_ASSERT_DELTA(TestMiddle.y, 1 / 3.0, EPS);
    TS_ASSERT_DELTA(TestMiddle.z, 1.0, EPS);
  }


  void testTriangleQuadraturePoints() {
    std::shared_ptr<double []> Weights = nullptr;
    std::shared_ptr<double []> PointsData = nullptr;
    unsigned NumPoints{};

    // Coordinates are taken from the Fortran implementation
    double ChiFortran[] = {0.94373743946307787,0.94373743946307787,0.94373743946307787,0.94373743946307787,
                           0.94373743946307787,0.94373743946307787,0.94373743946307787,0.81975930826310761,
                           0.81975930826310761,0.81975930826310761,0.81975930826310761,0.81975930826310761,
                           0.81975930826310761,0.81975930826310761,0.64737528288683033,0.64737528288683033,
                           0.64737528288683033,0.64737528288683033,0.64737528288683033,0.64737528288683033,
                           0.64737528288683033,0.45284637366944464,0.45284637366944464,0.45284637366944464,
                           0.45284637366944464,0.45284637366944464,0.45284637366944464,0.45284637366944464,
                           0.26578982278458946,0.26578982278458946,0.26578982278458946,0.26578982278458946,
                           0.26578982278458946,0.26578982278458946,0.26578982278458946,0.11467905316090421,
                           0.11467905316090421,0.11467905316090421,0.11467905316090421,0.11467905316090421,
                           0.11467905316090421,0.11467905316090421,2.2479386438712501E-002,2.2479386438712501E-002,
                           2.2479386438712501E-002,2.2479386438712501E-002,2.2479386438712501E-002,
                           2.2479386438712501E-002, 2.2479386438712501E-002};

    double TauFortran[] = {5.4830900955589179E-002,4.8991501878361855E-002,3.9548223967454631E-002,
                           2.8131280268461067E-002,1.6714336569467501E-002,7.2710586585602805E-003,
                           1.4316595813329493E-003,0.17565427919525450,0.15694739278690259,0.12669525127960912,
                           9.0120345868446194E-002,5.3545440457283260E-002,2.3293298949989799E-002,
                           4.5864125416378871E-003,0.34365181310645293,0.30705347083287471,0.24786787440468791,
                           0.17631235855658484,0.10475684270848173,4.5571246280294943E-002,8.9729040067167108E-003,
                           0.53323073117395925,0.47644255178423012,0.38460663631768571,0.27357681316527771,
                           0.16254699001286968,7.0711074546325303E-002,1.3922895156596097E-002,0.71552743286656784,
                           0.63932496020254781,0.51609290886511228,0.36710508860770530,0.21811726835029835,
                           9.4885217012862830E-002,1.8682744348842751E-002,0.86279303122343209,0.77090701909233450,
                           0.62231208026329454,0.44266047341954790,0.26300886657580119,0.11441392774676130,
                           2.2527915615663658E-002,0.95264658118522672,0.85119131654161828,0.68712130747329714,
                           0.48876030678064375,0.29039930608799031,0.12632929701966925,2.4874032376060777E-002};

    std::tie(NumPoints, Weights, PointsData) = generateTriangleQuadrature(7);
    double (*TestTrianglePoints)[2] = reshape<2>(&PointsData[0]);

    const  double EPS = 1e-6;
    for (unsigned i = 0; i < NumPoints; ++i) {
      TS_ASSERT_DELTA(TestTrianglePoints[i][0], ChiFortran[i], EPS);
      TS_ASSERT_DELTA(TestTrianglePoints[i][1], TauFortran[i], EPS);
    }
  }

  void testStrikeAndDipVectors() {
    VrtxCoords TestNormal{-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    VrtxCoords TestStrike{0.0, 0.0, 0.0};
    VrtxCoords TestDip{0.0, 0.0, 0.0};
    computeStrikeAndDipVectors(TestNormal, TestStrike, TestDip);

    // compute expected Strike results
    Eigen::Vector3d e3(0.0, 0.0, -1.0);
    Eigen::Vector3d Normal(TestNormal[0], TestNormal[1], TestNormal[2]);
    Eigen::Vector3d ResultStrike = e3.cross(Normal).normalized();


    const  double EPS = 1e-6;
    for (unsigned i = 0; i < 3; ++i) {
      TS_ASSERT_DELTA(TestStrike[i], ResultStrike(i), EPS);
    }
    // compute expected Dip results
    Eigen::Vector3d ResultDip = Normal.cross(ResultStrike);
    for (unsigned i = 0; i < 3; ++i) {
      TS_ASSERT_DELTA(TestDip[i], ResultDip(i), EPS);
    }
  }


  void testXiEtaZeta2chiTau() {
    const  double EPS = 1e-6;
    double TestChiTau[2] = {0.0, 0.0};
    {
      unsigned Face = 0;
      VrtxCoords xiEtaZeta{0.25, 0.1, 0.0};
      transformations::XiEtaZeta2chiTau(Face, xiEtaZeta, TestChiTau);
      TS_ASSERT_DELTA(TestChiTau[0], 0.1, EPS);
      TS_ASSERT_DELTA(TestChiTau[1], 0.25, EPS);
    }
    {
      unsigned Face = 1;
      VrtxCoords xiEtaZeta{0.1, 0.0, 0.25};
      transformations::XiEtaZeta2chiTau(Face, xiEtaZeta, TestChiTau);
      TS_ASSERT_DELTA(TestChiTau[0], 0.1, EPS);
      TS_ASSERT_DELTA(TestChiTau[1], 0.25, EPS);
    }
    {
      unsigned Face = 2;
      VrtxCoords xiEtaZeta{0.0, 0.1, 0.25};
      transformations::XiEtaZeta2chiTau(Face, xiEtaZeta, TestChiTau);
      TS_ASSERT_DELTA(TestChiTau[0], 0.25, EPS);
      TS_ASSERT_DELTA(TestChiTau[1], 0.1, EPS);
    }
    {
      unsigned Face = 3;
      VrtxCoords xiEtaZeta{1/3.0, 1/3.0, 1/3.0}; // center of the 4th face (triangle in 3D space)
      transformations::XiEtaZeta2chiTau(Face, xiEtaZeta, TestChiTau);
      TS_ASSERT_DELTA(TestChiTau[0], 1 / 3.0, EPS);
      TS_ASSERT_DELTA(TestChiTau[1], 1 / 3.0, EPS);
    }
    {
      unsigned Face = 3;
      ExtVrtxCoords xiEtaZeta{0.0, -0.15, 0.15};
      ExtTriangle FourthFace(ExtVrtxCoords{1.0, 0.0, 0.0},
                             ExtVrtxCoords{0.0, 1.0, 0.0},
                             ExtVrtxCoords{0.0, 0.0, 1.0});
      VrtxCoords NormalDirection{1.0, 1.0, 1.0};
      projectPointToFace(xiEtaZeta, FourthFace, NormalDirection);

      transformations::XiEtaZeta2chiTau(Face, xiEtaZeta.Coords, TestChiTau);
      TS_ASSERT_DELTA(TestChiTau[0], xiEtaZeta.eta, EPS);
      TS_ASSERT_DELTA(TestChiTau[1], xiEtaZeta.zeta, EPS);

    }
  }


  void testBasisFunctions() {

    VrtxCoords Point{0.25, 0.25, 0.0};

    // placing two elements in such a way that basis functions on both sides end up being the same
    VrtxCoords PlusElementCoords[4]{{2.0, 0.0, 0.0},
                                    {0.0, 2.0, 0.0},
                                    {0.0, 0.0, 0.0},
                                    {0.0, 0.0, 2.0}};

    VrtxCoords MinusElementCoords[4]{{2.0, 0.0, 0.0},
                                     {0.0, 2.0, 0.0},
                                     {0.0, 0.0, 0.0},
                                     {0.0, 0.0, -2.0}};

    auto BasisFunctions = getPlusMinusBasisFunctions(Point, PlusElementCoords, MinusElementCoords);

    const  double EPS = 1e-6;
    for (unsigned i = 0; i < BasisFunctions.PlusSide.size(); ++i) {
      TS_ASSERT_DELTA(BasisFunctions.PlusSide[i], BasisFunctions.MinusSide[i], EPS);
    }
  }


  void testIsElementInside() {

    Eigen::Vector3d Points[3] = {{0.25, 0.25, 0.25},
                                 {0.5, 0.5, 0.5},
                                 {0.75, 0.75, 0.1}};
    unsigned numPoints = 3;
    short contained[3] = {0, 0, 0};
    unsigned meshId[3] = {std::numeric_limits<unsigned>::max(),
                          std::numeric_limits<unsigned>::max(),
                          std::numeric_limits<unsigned>::max()};

    std::vector<Vertex> Vertices;
    Vertices.push_back({{0.0, 0.0, 0.0}, {0}});
    Vertices.push_back({{1.0, 0.0, 0.0}, {0, 1}});
    Vertices.push_back({{0.0, 1.0, 0.0}, {0, 1}});
    Vertices.push_back({{0.0, 0.0, 1.0}, {0}});
    Vertices.push_back({{1.0, 1.0, 0.0}, {1}});
    Vertices.push_back({{1.0, 1.0, 1.0}, {1}});

    std::vector<Element> Elements;
    Element E1;
    E1.localId = 0;
    E1.vertices[0] = 0; E1.vertices[1] = 1; E1.vertices[2] = 2; E1.vertices[3] = 3;
    Elements.push_back(E1);

    Element E2;
    E2.localId = 1;
    E2.vertices[0] = 1; E2.vertices[1] = 4; E2.vertices[2] = 2; E2.vertices[3] = 5;
    Elements.push_back(E2);


    initializers::findMeshIds(Points,
                              Vertices,
                              Elements,
                              numPoints,
                              contained,
                              meshId);

    TS_ASSERT_EQUALS(contained[0], 1);
    TS_ASSERT_EQUALS(contained[1], 0);
    TS_ASSERT_EQUALS(contained[2], 1);

    TS_ASSERT_EQUALS(meshId[0], 0);
    TS_ASSERT_EQUALS(meshId[1], std::numeric_limits<unsigned>::max());
    TS_ASSERT_EQUALS(meshId[2], 1);
  }
};

#endif //SEISSOL_GEOMETRY_HPP