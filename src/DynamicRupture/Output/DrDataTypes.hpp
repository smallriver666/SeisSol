#ifndef SEISSOL_DROUTOUT_DRDATATYPES_HPP
#define SEISSOL_DROUTOUT_DRDATATYPES_HPP

#include "Kernels/precision.hpp"
#include <vector>
#include <array>

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
        real PrintTimeIntervalSec{1.0};
        int PrintIntervalCriterion{1};
        int MaxPickStore{50};
        std::array<bool, 12> OutputMask{true, true, true, true};
        int RefinementStrategy{2};
        int Refinement{2};
      };
    }


    struct GeoCoordsT {
      GeoCoordsT() = default;
      ~GeoCoordsT() = default;
      GeoCoordsT(const GeoCoordsT& Other) {
        for (int i = 0; i < 3; ++i)
          Values[i] = Other.Values[i];
      }
      GeoCoordsT& operator=(const GeoCoordsT& Other) {
        for (int i = 0; i < 3; ++i)
          Values[i] = Other.Values[i];
        return *this;
      }

      std::array<real, 3> Values = {0.0, 0.0, 0.0};
      double& x = Values[0];
      double& y = Values[1];
      double& z = Values[2];

      double& xi = Values[0];
      double& eta = Values[1];
      double& zeta = Values[2];
    };


    struct TriangleT {
      TriangleT() = default;
      ~TriangleT() = default;
      explicit TriangleT(GeoCoordsT p1, GeoCoordsT p2, GeoCoordsT p3) {
        Points[0] = p1;
        Points[1] = p2;
        Points[3] = p3;
      }

      TriangleT(const TriangleT& Other) {
        for (int i = 0; i < 3; ++i)
          Points[i] = Other.Points[i];
      }
      TriangleT& operator=(const TriangleT& Other) {
        for (int i = 0; i < 3; ++i)
          Points[i] = Other.Points[i];
        return *this;
      }

      std::array<GeoCoordsT, 3> Points{};
      GeoCoordsT& p1 = Points[0];
      GeoCoordsT& p2 = Points[1];
      GeoCoordsT& p3 = Points[2];
    };


    struct ReceiverPointT {
      GeoCoordsT Global{};   // physical coords of a receiver
      GeoCoordsT Referece{}; // reference coords of a receiver
      TriangleT GlobalSubTet{};// (subtet) vertices coordinates (of a surrounding triangle)
      int FaultFaceIndex{-1};  // Face Fault index which the receiver belongs to
      int ElementIndex{-1};    // Element which the receiver belongs to
      int GlobalReceiverIndex{-1};  // receiver index of global list
      bool IsInside{};  // If a point is inside the mesh or not
    };
    using ReceiverPointsT = std::vector<ReceiverPointT>;


    struct ConstantT {
      real P0{0.0};
      real Ts0{0.0};
      real Td0{0.0};
    };
    using ConstantsT = std::vector<ConstantT>;
  }
}


#endif //SEISSOL_DROUTOUT_DRDATATYPES_HPP
