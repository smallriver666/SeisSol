#ifndef SEISSOL_DROUTOUT_DRINITIALIZER_HPP
#define SEISSOL_DROUTOUT_DRINITIALIZER_HPP

#include "DrDataTypes.hpp"
#include <yaml-cpp/yaml.h>

namespace seissol {
  namespace dr {
    namespace output {
      class ParametersInitializer;
    }
  }
}

class seissol::dr::output::ParametersInitializer {
public:
  explicit ParametersInitializer(const YAML::Node& Data) : m_Data(Data) {}


  GeneralParamsT getDrGeneralParams() {
    using namespace initializers;
    GeneralParamsT Params{};

    if (!m_Data["dynamicrupture"]) {
      throw std::runtime_error("dynamic rupture params. is not provided");
    }

    const YAML::Node& DrData = m_Data["dynamicrupture"];

    Params.OutputPointType = static_cast<OutputType>(getParamIfExists(DrData, "outputpointtype", 3));
    Params.SlipRateOutputType = getParamIfExists(DrData, "sliprateoutputtype", 1);
    Params.FrictionLawType = getParamIfExists(DrData, "fl", 0);
    Params.BackgroundType = getParamIfExists(DrData, "backgroundtype", 0);
    Params.IsRfOutputOn = getParamIfExists(DrData, "rf_output_on", false);
    Params.IsDsOutputOn = getParamIfExists(DrData, "ds_output_on", false);
    Params.IsMagnitudeOutputOn = getParamIfExists(DrData, "magnitude_output_on", false);
    Params.IsEnergyRateOutputOn = getParamIfExists(DrData, "energy_rate_output_on", false);
    Params.IsGpWiseOutput = getParamIfExists(DrData, "gpwise", false);
    Params.IsTermalPressureOn = getParamIfExists(DrData, "thermalpress", false);
    Params.BackgroundType = getParamIfExists(DrData, "energy_rate_printtimeinterval", 1);

    return Params;
  }


  PickpointParamsT getPickPointParams() {
    using namespace initializers;
    PickpointParamsT PpParams{};


    if (!m_Data["pickpoint"]) {
      throw std::runtime_error("pickpoint output parameters for dynamic rupture is not provided");
    }

    const YAML::Node& PpData = m_Data["pickpoint"];
    PpParams.PrintTimeInterval = getParamIfExists(PpData, "printtimeinterval", 1);
    PpParams.NumOutputPoints = getParamIfExists(PpData, "noutpoints", 0);
    PpParams.PPFileName = getParamIfExists(PpData, "ppfilename", std::string());

    if (m_Data["outputmask"]) {
      convertStringToMask(m_Data["outputmask"].as<std::string>(), PpParams.OutputMask);
    }

    return PpParams;
  }


  ElementwiseFaultParamsT getElementwiseFaultParams() {
    using namespace initializers;

    ElementwiseFaultParamsT EwParams{};

    if (!m_Data["elementwise"]) {
      throw std::runtime_error("elementwise fault output parameters for dynamic rupture is not provided");
    }

    const YAML::Node& EwData = m_Data["elementwise"];

    EwParams.PrintTimeInterval = getParamIfExists(EwData, "printtimeinterval", 2);
    EwParams.PrintTimeIntervalSec = getParamIfExists(EwData, "printtimeinterval_sec", 1.0);
    EwParams.PrintIntervalCriterion = getParamIfExists(EwData, "printintervalcriterion", 1);
    EwParams.RefinementStrategy = getParamIfExists(EwData, "refinement_strategy", 2);
    EwParams.Refinement = getParamIfExists(EwData, "refinement", 2);

    if (EwData["outputmask"]) {
      convertStringToMask(EwData["outputmask"].as<std::string>(), EwParams.OutputMask);
    }

    return EwParams;
  }


private:
  const YAML::Node& m_Data;
};

#endif //SEISSOL_DRINITIALIZER_HPP
