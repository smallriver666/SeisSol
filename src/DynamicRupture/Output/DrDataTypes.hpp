#ifndef SEISSOL_DROUTOUT_DRDATATYPES_HPP
#define SEISSOL_DROUTOUT_DRDATATYPES_HPP

#include "Kernels/precision.hpp"
#include "Geometry/MeshDefinition.h"
#include "Initializer/tree/Layer.hpp"
#include <vector>
#include <array>
#include <tuple>
#include <cassert>
#include <limits>

namespace seissol {
  namespace dr {
    namespace output {

      template <int Dim>
      constexpr bool isMultiDimVar() {
        return (Dim > 1);
      }

      template<int DIM>
      struct VarT {
        constexpr int dim() {return DIM;}


        real* operator[](int index) {
          assert(index < DIM && "access is out of the bounds");
          return data[index];
        }


        void releaseData() {
          if (isActive) {
            for (auto item: data) {
              delete[] item;
            }
          }
        }

        std::array<real*, DIM> data = {nullptr};
        bool isActive{false};
        size_t size{};
      };

      using Var1D = VarT<1>;
      using Var2D = VarT<2>;
      using Var3D = VarT<3>;
      using DrVarsT = std::tuple<Var2D, Var3D, Var1D, Var2D, Var3D, Var2D, Var1D, Var1D, Var1D, Var1D, Var1D, Var2D>;

      enum DirectionID {STRIKE = 0, DIP = 1, NORMAL = 2};
      enum ThermoID {PRESSURE = 0, TEMPERATURE = 1};
      enum ParamID {FUNCTION = 0, STATE = 1};

      enum VariableID {SlipRate = 0,
                       TransientShierStress,
                       NormalVelocity,
                       FunctionAndState,
                       TotalStresses,
                       Slip,
                       RuptureVelocity,
                       AbsoluteSlip,
                       PeakSlipsRate,
                       RuptureTime,
                       Ds,
                       Thermo,
                       Size};

      enum class OutputType : int {None = 0,
                                   AtPickpoint = 3,
                                   Elementwise = 4,
                                   AtPickpointAndElementwise = 5};


      struct GeneralParamsT {
        OutputType outputPointType{OutputType::AtPickpoint};
        int slipRateOutputType{1};
        int frictionLawType{0};
        int backgroundType{0};
        bool isRfOutputOn{false};
        bool isDsOutputOn{false};
        bool isMagnitudeOutputOn{false};
        bool isEnergyRateOutputOn{false};
        bool isGpWiseOutput{false};
        bool isTermalPressureOn{false};
        int energyRatePrintTimeInterval{1};
        bool isRfTimeOn{false};
        bool faultOutputFlag {false};
        std::string outputFilePrefix{"data"};
        std::string xdmfWriterBackend{"hdf5"};
        std::string checkPointBackend{"none"};
        real endTime{0.0};
        size_t maxIteration{10000000};
      };


      struct PickpointParamsT {
        std::array<bool, std::tuple_size<DrVarsT>::value> outputMask{true, true, true}; // the rest is false by default
        int printTimeInterval{1};
        int numOutputPoints{0};
        int maxPickStore{50};
        std::string ppFileName{};
      };


      struct ElementwiseFaultParamsT {
        int printTimeInterval{2};
        double printTimeIntervalSec{1.0};
        int printIntervalCriterion{1};
        int maxPickStore{50};
        std::array<bool, std::tuple_size<DrVarsT>::value> outputMask{true, true, true, true};
        int refinementStrategy{2};
        int refinement{2};
      };
    }


    struct ExtVrtxCoords {
      ExtVrtxCoords() = default;
      ~ExtVrtxCoords() = default;
      ExtVrtxCoords(const ExtVrtxCoords& other) {
        for (int i = 0; i < 3; ++i)
          coords[i] = other.coords[i];
      }
      ExtVrtxCoords& operator=(const ExtVrtxCoords& other) {
        for (int i = 0; i < 3; ++i)
          coords[i] = other.coords[i];
        return *this;
      }
      ExtVrtxCoords(std::initializer_list<double> inputCoords) {
        assert(inputCoords.size() == 3 && "ExtVrtxCoords must get initialized with 3 values");
        auto begin = inputCoords.begin();
        for (int i = 0; i < 3; ++i, ++begin)
          coords[i] = *begin;
      }

      VrtxCoords coords = {0.0, 0.0, 0.0};
      double& x = coords[0];
      double& y = coords[1];
      double& z = coords[2];

      double& xi = coords[0];
      double& eta = coords[1];
      double& zeta = coords[2];
  };


    struct ExtTriangle {
      ExtTriangle() = default;
      ~ExtTriangle() = default;
      explicit ExtTriangle(const ExtVrtxCoords& p1, const ExtVrtxCoords& p2, const ExtVrtxCoords& p3) {
        points[0] = p1;
        points[1] = p2;
        points[2] = p3;
      }

      ExtTriangle(const ExtTriangle& other) {
        for (int i = 0; i < 3; ++i)
          points[i] = other.points[i];
      }
      ExtTriangle& operator=(const ExtTriangle& other) {
        for (int i = 0; i < 3; ++i)
          points[i] = other.points[i];
        return *this;
      }

      std::array<ExtVrtxCoords, 3> points{};
      ExtVrtxCoords& p1 = points[0];
      ExtVrtxCoords& p2 = points[1];
      ExtVrtxCoords& p3 = points[2];
    };


    struct ReceiverPointT {
      ExtVrtxCoords global{};    // physical coords of a receiver
      ExtVrtxCoords referece{};  // reference coords of a receiver
      ExtTriangle globalSubTet{};// (subtet) vertices coordinates (of a surrounding triangle)
      int faultFaceIndex{-1};    // Face Fault index which the receiver belongs to
      int localFaceSideId{-1};   // Side ID of a reference element
      int elementIndex{-1};      // Element which the receiver belongs to
      int globalReceiverIndex{-1};  // receiver index of global list
      bool isInside{};  // If a point is inside the mesh or not
      int nearestGpIndex{-1};
      double distanceToNearestGp{std::numeric_limits<double>::max()};
    };
    using ReceiverPointsT = std::vector<ReceiverPointT>;

    struct FaultDirectionsT {
      const double* faceNormal{};
      const double* tangent1{};
      const double* tangent2{};
      VrtxCoords strike{0.0, 0.0, 0.0};
      VrtxCoords dip{0.0, 0.0, 0.0};
    };

    struct PlusMinusBasisFunctionsT {
      std::vector<real> plusSide;
      std::vector<real> minusSide;
    };

    struct ConstantT {
      real p0{0.0};
      real ts0{0.0};
      real td0{0.0};
    };
    using ConstantsT = std::vector<ConstantT>;

    struct OutputData {
      output::DrVarsT vars;
      std::vector<PlusMinusBasisFunctionsT> basisFunctions;
      std::vector<ReceiverPointT> receiverPoints;
      std::vector<std::vector<real>> rotationMatrices;
      std::vector<FaultDirectionsT> faultDirections{};
      std::vector<ConstantT> constrains;
      std::vector<size_t> currentPick{};
      std::vector<real> cachedTime{};
      size_t maxPickStore;
    };
  }
}


#endif //SEISSOL_DROUTOUT_DRDATATYPES_HPP
