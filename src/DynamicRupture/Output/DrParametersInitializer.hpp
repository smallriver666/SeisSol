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
  explicit ParametersInitializer(const YAML::Node& userData) : data(userData) {}


  GeneralParamsT getDrGeneralParams() {
    using namespace initializers;
    GeneralParamsT params{};

    if (!data["dynamicrupture"]) {
      throw std::runtime_error("dynamic rupture params. is not provided");
    }

    const YAML::Node& DrData = data["dynamicrupture"];
    params.outputPointType = static_cast<OutputType>(getParamIfExists(DrData, "outputpointtype", 3));
    params.slipRateOutputType = getParamIfExists(DrData, "sliprateoutputtype", 1);
    params.frictionLawType = getParamIfExists(DrData, "fl", 0);
    params.backgroundType = getParamIfExists(DrData, "backgroundtype", 0);
    params.isRfOutputOn = getParamIfExists(DrData, "rf_output_on", false);
    params.isDsOutputOn = getParamIfExists(DrData, "ds_output_on", false);
    params.isMagnitudeOutputOn = getParamIfExists(DrData, "magnitude_output_on", false);
    params.isEnergyRateOutputOn = getParamIfExists(DrData, "energy_rate_output_on", false);
    params.isGpWiseOutput = getParamIfExists(DrData, "gpwise", false);
    params.isTermalPressureOn = getParamIfExists(DrData, "thermalpress", false);
    params.backgroundType = getParamIfExists(DrData, "energy_rate_printtimeinterval", 1);


    const YAML::Node& OutputParams = data["output"];
    params.faultOutputFlag = getParamIfExists(OutputParams, "faultoutputflag", false);
    params.outputFilePrefix = getParamIfExists(OutputParams, "outputfile", std::string("data"));
    params.checkPointBackend = getParamIfExists(OutputParams, "checkpointbackend", std::string("none"));

#ifdef USE_HDF
    params.xdmfWriterBackend = getParamIfExists(OutputParams, "xdmfwriterbackend", std::string("hdf5"));
#else
    params.xdmfWriterBackend = getParamIfExists(OutputParams, "xdmfwriterbackend", std::string("posix"));
#endif

    return params;
  }


  PickpointParamsT getPickPointParams() {
    using namespace initializers;
    PickpointParamsT ppParams{};


    if (!data["pickpoint"]) {
      throw std::runtime_error("pickpoint output parameters for dynamic rupture is not provided");
    }

    const YAML::Node& ppData = data["pickpoint"];
    ppParams.printTimeInterval = getParamIfExists(ppData, "printtimeinterval", 1);
    ppParams.numOutputPoints = getParamIfExists(ppData, "noutpoints", 0);
    ppParams.ppFileName = getParamIfExists(ppData, "ppfilename", std::string());

    if (data["outputmask"]) {
      convertStringToMask(data["outputmask"].as<std::string>(), ppParams.OutputMask);
    }

    return ppParams;
  }


  ElementwiseFaultParamsT getElementwiseFaultParams() {
    using namespace initializers;

    ElementwiseFaultParamsT ewParams{};

    if (!data["elementwise"]) {
      throw std::runtime_error("elementwise fault output parameters for dynamic rupture is not provided");
    }

    const YAML::Node& ewData = data["elementwise"];

    ewParams.printTimeInterval = getParamIfExists(ewData, "printtimeinterval", 2);
    ewParams.printTimeIntervalSec = getParamIfExists(ewData, "printtimeinterval_sec", 1.0);
    ewParams.printIntervalCriterion = getParamIfExists(ewData, "printintervalcriterion", 1);
    ewParams.maxPickStore = getParamIfExists(ewData, "maxpickstore", 50);
    ewParams.refinementStrategy = getParamIfExists(ewData, "refinement_strategy", 2);
    ewParams.refinement = getParamIfExists(ewData, "refinement", 2);

    if (ewData["outputmask"]) {
      convertStringToMask(ewData["outputmask"].as<std::string>(), ewParams.outputMask);
    }

    return ewParams;
  }


private:
  const YAML::Node& data;
};

#endif //SEISSOL_DRINITIALIZER_HPP
