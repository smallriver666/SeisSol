#include "FaultRefiner.hpp"

namespace seissol {
  namespace dr {
    namespace output {

      std::unique_ptr<RefinerInterface> getRefiner(const int Strategy) {
        switch (Strategy) {
          case 1:
            return std::unique_ptr<RefinerInterface>(new TripleFaultFaceRefiner);
            break;
          case 2:
            return std::unique_ptr<RefinerInterface>(new QuadFaultFaceRefiner);
            break;
          default:
            throw ("Unknown refinement strategy for Fault Face Refiner");
        }
      }

    }
  }
}