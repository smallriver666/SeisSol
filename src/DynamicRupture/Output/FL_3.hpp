#ifndef SEISSOL_DROUTOUT_FL_3_HPP
#define SEISSOL_DROUTOUT_FL_3_HPP

#include "DynamicRupture/Output/DrBase.hpp"

namespace seissol {
  namespace dr {
    namespace output {
      class FL_3;
    }
  }
}

class seissol::dr::output::FL_3 : public seissol::dr::output::Base {
public:
  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) override {
    seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
    std::cout << "tie ptr for FL_33\n";
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
    seissol::initializers::DR_FL_33 &ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 &>(DynRup);
    std::cout << "output vars for FL_33\n";
  }
};

#endif //SEISSOL_DROUTOUT_FL_3_HPP
