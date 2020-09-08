#ifndef SEISSOL_DROUTOUT_DRDATATYPES_HPP
#define SEISSOL_DROUTOUT_DRDATATYPES_HPP

#include <vector>

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
        std::array<bool, 12> OutputMask{true, true, true, true};
        int RefinementStrategy{2};
        int Refinement{2};
      };

    }

    struct PointT {
      real x{}, y{}, z{}; // physical coordinates
      real xi{}, eta{}, zeta{};  // reference coordinates
      std::array<real, 3> coordX{0.0, 0.0, 0.0}; // (subtet) vertices coordinates in x
      std::array<real, 3> coordY{0.0, 0.0, 0.0}; // (subtet) vertices coordinates in y
      std::array<real, 3> coordZ{0.0, 0.0, 0.0}; // (subtet) vertices coordinates in z
      int Index{-1};  // Element index
      int GlobalReceiverIndex{-1};  // receiver index of global list
      bool IsInside{};  // If a point is inside the mesh or not
    };

    using PointsT = std::vector<PointT>;

  }
}


#endif //SEISSOL_DROUTOUT_DRDATATYPES_HPP
