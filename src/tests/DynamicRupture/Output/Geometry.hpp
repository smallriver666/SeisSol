#ifndef SEISSOL_GEOMETRY_HPP
#define SEISSOL_GEOMETRY_HPP

#include <cxxtest/TestSuite.h>
#include "DynamicRupture/Output/DrOutputAux.hpp"
#include <iostream>

namespace seissol {
  namespace unit_test {
    namespace dr {
      class Geometry;
    }
  }
}

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
};

#endif //SEISSOL_GEOMETRY_HPP