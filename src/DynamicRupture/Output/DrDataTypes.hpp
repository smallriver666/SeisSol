#ifndef SEISSOL_DROUTOUT_DRDATATYPES_HPP
#define SEISSOL_DROUTOUT_DRDATATYPES_HPP

#include "Kernels/precision.hpp"
#include "Geometry/MeshDefinition.h"
#include <vector>
#include <array>
#include <cassert>
#include <limits>

namespace seissol {
  namespace dr {
    namespace output {
      enum class OutputType : int {None = 0,
                                   AtPickpoint = 3,
                                   Elementwise = 4,
                                   AtPickpointAndElementwise = 5};


      struct GeneralParamsT {
        OutputType OutputPointType{OutputType::AtPickpoint};
        int SlipRateOutputType{1};
        int FrictionLawType{0};
        int BackgroundType{0};
        bool IsRfOutputOn{false};
        bool IsDsOutputOn{false};
        bool IsMagnitudeOutputOn{false};
        bool IsEnergyRateOutputOn{false};
        bool IsGpWiseOutput{false};
        bool IsTermalPressureOn{false};
        int EnergyRatePrintTimeInterval{1};
        bool IsRfTimeOn{false};
        bool FaultOutputFlag {false};
        std::string OutputFilePrefix{"data"};
        std::string XdmfWriterBackend{"hdf5"};
        std::string CheckPointBackend{"none"};
      };


      struct PickpointParamsT {
        std::array<bool, 12> OutputMask{true, true, true}; // the rest is false by default
        int PrintTimeInterval{1};
        int NumOutputPoints{0};
        std::string PPFileName{};
      };


      struct ElementwiseFaultParamsT {
        int PrintTimeInterval{2};
        double PrintTimeIntervalSec{1.0};
        int PrintIntervalCriterion{1};
        int MaxPickStore{50};
        std::array<bool, 12> OutputMask{true, true, true, true};
        int RefinementStrategy{2};
        int Refinement{2};
      };
    }


    struct ExtVrtxCoords {
      ExtVrtxCoords() = default;
      ~ExtVrtxCoords() = default;
      ExtVrtxCoords(const ExtVrtxCoords& Other) {
        for (int i = 0; i < 3; ++i)
          Coords[i] = Other.Coords[i];
      }
      ExtVrtxCoords& operator=(const ExtVrtxCoords& Other) {
        for (int i = 0; i < 3; ++i)
          Coords[i] = Other.Coords[i];
        return *this;
      }
      ExtVrtxCoords(std::initializer_list<double> InputCoords) {
        assert(InputCoords.size() == 3 && "ExtVrtxCoords must get initialized with 3 values");
        auto Begin = InputCoords.begin();
        for (int i = 0; i < 3; ++i, ++Begin)
          Coords[i] = *Begin;
      }

      VrtxCoords Coords = {0.0, 0.0, 0.0};
      double& x = Coords[0];
      double& y = Coords[1];
      double& z = Coords[2];

      double& xi = Coords[0];
      double& eta = Coords[1];
      double& zeta = Coords[2];
  };


    struct ExtTriangle {
      ExtTriangle() = default;
      ~ExtTriangle() = default;
      explicit ExtTriangle(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2, const ExtVrtxCoords& p3) {
        Points[0] = p1;
        Points[1] = p2;
        Points[2] = p3;
      }

      ExtTriangle(const ExtTriangle& Other) {
        for (int i = 0; i < 3; ++i)
          Points[i] = Other.Points[i];
      }
      ExtTriangle& operator=(const ExtTriangle& Other) {
        for (int i = 0; i < 3; ++i)
          Points[i] = Other.Points[i];
        return *this;
      }

      std::array<ExtVrtxCoords, 3> Points{};
      ExtVrtxCoords& p1 = Points[0];
      ExtVrtxCoords& p2 = Points[1];
      ExtVrtxCoords& p3 = Points[2];
    };


    struct ReceiverPointT {
      ExtVrtxCoords Global{};    // physical coords of a receiver
      ExtVrtxCoords Referece{};  // reference coords of a receiver
      ExtTriangle GlobalSubTet{};// (subtet) vertices coordinates (of a surrounding triangle)
      int FaultFaceIndex{-1};    // Face Fault index which the receiver belongs to
      int LocalFaceSideId{-1};   // Side ID of a reference element
      int ElementIndex{-1};      // Element which the receiver belongs to
      int GlobalReceiverIndex{-1};  // receiver index of global list
      bool IsInside{};  // If a point is inside the mesh or not
      int NearestGpIndex{-1};
      double DistanceToNearestGp{std::numeric_limits<double>::max()};
    };
    using ReceiverPointsT = std::vector<ReceiverPointT>;

    struct PlusMinusBasisFunctionsT {
      std::vector<real> PlusSide;
      std::vector<real> MinusSide;
    };

    struct ConstantT {
      real P0{0.0};
      real Ts0{0.0};
      real Td0{0.0};
    };
    using ConstantsT = std::vector<ConstantT>;
  }
}


#endif //SEISSOL_DROUTOUT_DRDATATYPES_HPP
