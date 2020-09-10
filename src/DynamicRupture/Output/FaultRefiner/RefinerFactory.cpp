#include "FaultRefiner.hpp"

namespace seissol {
  namespace dr {
    namespace output {

      std::unique_ptr<FaultRefinerInterface> getRefiner(const int Strategy) {
        switch (Strategy) {
          case 1:
            return std::unique_ptr<FaultRefinerInterface>(new TripleFaultFaceRefiner);
          case 2:
            return std::unique_ptr<FaultRefinerInterface>(new QuadFaultFaceRefiner);
          default:
            throw ("Unknown refinement strategy for Fault Face Refiner");
        }
      }

    }
  }
}