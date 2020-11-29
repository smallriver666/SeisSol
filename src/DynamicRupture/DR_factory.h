//
// Created by adrian on 08.07.20.
//

#ifndef SEISSOL_DR_FACTORY_H
#define SEISSOL_DR_FACTORY_H

#include <iostream>
#include <tuple>
#include <stdexcept>
#include <Initializer/DynamicRupture.h>
#include "DR_initializer_base.h"
#include "DR_solver_base.h"
#include "DR_output.h"
#include "DR_solver_rate_and_state.h"
#include "DR_solver_linear.h"
#include "DR_solver_legacy_RS.h"


namespace seissol {
  namespace dr {
    namespace factory {
      class AbstractFactory;
      struct Factory_FL_0;
      struct Factory_FL_2;
      struct Factory_FL_3; //aging law
      struct Factory_FL_4; //slip law
      struct Factory_FL_6;
      struct Factory_FL_7;
      struct Factory_FL_16;
      struct Factory_FL_33;
      struct Factory_FL_103;
      struct Factory_FL_103_Thermal;
      seissol::dr::factory::AbstractFactory* getFactory(dr::DrParameterT *DynRupParameter);
    }
  }
}


using products = std::tuple<seissol::initializers::DynamicRupture*,seissol::initializers::BaseDrInitializer*, seissol::dr::fr_law::BaseFrictionSolver*, seissol::dr::output::Output_Base*>;
class seissol::dr::factory::AbstractFactory {
public:
  virtual ~AbstractFactory() {}
  virtual products produce() = 0;
};
class seissol::dr::factory::Factory_FL_0 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DynamicRupture,
                           new seissol::initializers::Init_FL_0,
                           new seissol::dr::fr_law::SolverNoFaultFL0,
                           new seissol::dr::output::Output_FL_0);
  }
};

class seissol::dr::factory::Factory_FL_2 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_linear,
                           new seissol::initializers::Init_linear,
                           new seissol::dr::fr_law::LinearSlipWeakeningSolverFL2,
                           new seissol::dr::output::Output_FL_2);
  }
};
class seissol::dr::factory::Factory_FL_3 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_3,
                           new seissol::initializers::Init_FL_3,
                           new seissol::dr::fr_law::Solver_FL_3,
                           new seissol::dr::output::Output_FL_3);
  }
};

class seissol::dr::factory::Factory_FL_4 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_3,
                           new seissol::initializers::Init_FL_3,
                           new seissol::dr::fr_law::Solver_FL_4,
                           new seissol::dr::output::Output_FL_3);
  }
};

class seissol::dr::factory::Factory_FL_6 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_6,
                           new seissol::initializers::Init_FL_6,
                           new seissol::dr::fr_law::LinearSlipWeakSolverBimaterialFL6,
                           new seissol::dr::output::Output_FL_6);
  }
};


class seissol::dr::factory::Factory_FL_7 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_3,
                           new seissol::initializers::Init_FL_3,
                           new seissol::dr::fr_law::SolverRateAndStateVwFL7,
                           new seissol::dr::output::Output_FL_3);
  }
};

class seissol::dr::factory::Factory_FL_16 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_16,
                           new seissol::initializers::Init_FL_16,
                           new seissol::dr::fr_law::LinearSlipWeakeningSolverFL16,
                           new seissol::dr::output::Output_FL_2);
  }
};

class seissol::dr::factory::Factory_FL_33 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_33,
                           new seissol::initializers::Init_FL_33,
                           new seissol::dr::fr_law::SolverImposedSlipRatesFL33,
                           new seissol::dr::output::Output_FL_33);
  }
};

class seissol::dr::factory::Factory_FL_103 : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_103,
                           new seissol::initializers::Init_FL_103,
                           new seissol::dr::fr_law::RateAndStateNucFL103,
                           new seissol::dr::output::Output_FL_103);
  }
};

class seissol::dr::factory::Factory_FL_103_Thermal : public seissol::dr::factory::AbstractFactory {
  virtual products produce() override {
    return std::make_tuple(new seissol::initializers::DR_FL_103_Thermal,
                           new seissol::initializers::Init_FL_103_Thermal,
                           new seissol::dr::fr_law::RateAndStateThermalFL103,
                           new seissol::dr::output::Output_FL_103);
  }
};

seissol::dr::factory::AbstractFactory* seissol::dr::factory::getFactory(dr::DrParameterT *DynRupParameter);

int temporary_main();

#endif //SEISSOL_DR_FACTORY_H