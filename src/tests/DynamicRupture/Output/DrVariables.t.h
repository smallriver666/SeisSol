#ifndef SEISSOL_DRVARIABLES_T_H
#define SEISSOL_DRVARIABLES_T_H

#include <cxxtest/TestSuite.h>
#include <DynamicRupture/Output/DrDataTypes.hpp>
#include <DynamicRupture/DR_math.h>
#include <Initializer/MemoryAllocator.h>

namespace seissol {
  namespace unit_test {
    namespace dr {
      class DrVariables;
    }
  }
}



using namespace seissol;
using namespace seissol::dr;

class seissol::unit_test::dr::DrVariables : public CxxTest::TestSuite {
public:
  void testGeneralVariablesCount() {
    output::DrVarsT drVars;

    unsigned variableCounter = 0;
    auto countVariables = [&variableCounter](auto& var, int) {
      ++variableCounter;
    };

    aux::forEach(drVars, countVariables);
    TS_ASSERT_EQUALS(variableCounter, 12);
  }

  void testTotalVariablesCount() {
    output::DrVarsT drVars;

    std::array<bool, std::tuple_size<output::DrVarsT>::value> mask;
    for (size_t i = 0; i < std::tuple_size<output::DrVarsT>::value; ++i)
      mask[i] = true;

    auto assignMask = [&mask](auto& var, int index) {
      var.isActive = mask[index];
    };

    aux::forEach(drVars, assignMask);

    unsigned variableCounter = 0;
    auto countVariables = [&variableCounter](auto& var, int) {
      if (var.isActive) {
        variableCounter += var.dim();
      }
    };

    aux::forEach(drVars, countVariables);
    TS_ASSERT_EQUALS(variableCounter, 20);
  }

  void testAllocationDeallocationOfVariables() {
    output::DrVarsT drVars;

    std::array<bool, std::tuple_size<output::DrVarsT>::value> mask;
    for (size_t i = 0; i < std::tuple_size<output::DrVarsT>::value; ++i)
      mask[i] = true;

    auto assignMask = [&mask](auto& var, int index) {
      var.isActive = mask[index];
    };

    aux::forEach(drVars, assignMask);

    seissol::memory::ManagedAllocator allocator;
    const unsigned numElements = 1024;
    auto allocateVariables = [numElements, &allocator](auto& var, int) {
      if (var.isActive) {
        var.size = numElements;
        for (size_t dim = 0; dim < var.data.size(); ++dim) {
          var.data[dim] = reinterpret_cast<real*>(allocator.allocateMemory(numElements * sizeof(real),
                                                                           1,
                                                                           memory::Memkind::Standard));
        }
      }
      else {
        var.size = 0;
        for (size_t dim = 0; dim < var.data.size(); ++dim) {
          var.data[dim] = nullptr;
        }
      }
    };
    aux::forEach(drVars, allocateVariables);


    real initValue = 0.0;
    auto initVariables = [initValue](auto& var, int) {
      if (var.isActive) {
        for (size_t dim = 0; dim < var.data.size(); ++dim) {
          for (size_t i = 0; i < var.size; ++i) {
            var[dim][i] = initValue;
          }
        }
      }
    };
    aux::forEach(drVars, initVariables);
  }
};

#endif //SEISSOL_DRVARIABLES_T_H
