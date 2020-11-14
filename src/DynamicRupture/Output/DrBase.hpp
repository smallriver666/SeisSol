#ifndef SEISSOL_DROUTOUT_DRBASE_HPP
#define SEISSOL_DROUTOUT_DRBASE_HPP

#include "Initializer/InputAux.hpp"
#include "DynamicRupture/Output/DrParametersInitializer.hpp"
#include "DynamicRupture/Output/Builders/DrElementWiseOutput.hpp"
#include "DynamicRupture/Output/Builders/DrPickpointOutput.hpp"
#include <iostream>
#include <memory>

namespace seissol {
  namespace dr {
    namespace output {
      class Base;
    }
  }
}

class seissol::dr::output::Base {
public:
  virtual ~Base() {}

  void setInputParam(const YAML::Node& InputData, const MeshReader& Mesher) {
    using namespace initializers;

    ParametersInitializer Reader(InputData);
    generalParams = Reader.getDrGeneralParams();

    // adjust general output parameters
    generalParams.isRfTimeOn = generalParams.isRfOutputOn;
    if (generalParams.isDsOutputOn && !generalParams.isRfOutputOn) {
      generalParams.isRfOutputOn = true;
      generalParams.isRfTimeOn = true;
    }

    PickpointParamsT PpParams;
    ElementwiseFaultParamsT EwParams;
    switch (generalParams.outputPointType) {
      case OutputType::None:
        break;

      case OutputType::AtPickpoint:
        ppOutputBuilder = std::make_unique<PickpointOutput>();
        ppOutputBuilder->setParams(Reader.getPickPointParams(), &Mesher);
        break;

      case OutputType::Elementwise:
        ewOutputBuilder = std::make_unique<ElementWiseOutput>();
        ewOutputBuilder->setParams(Reader.getElementwiseFaultParams(), &Mesher);
        break;

      case OutputType::AtPickpointAndElementwise:
        ppOutputBuilder = std::make_unique<PickpointOutput>();
        ppOutputBuilder->setParams(Reader.getPickPointParams(), &Mesher);

        ewOutputBuilder = std::make_unique<ElementWiseOutput>();
        ewOutputBuilder->setParams(Reader.getElementwiseFaultParams(), &Mesher);
        break;

      default:
        throw std::runtime_error("Unknown fault output type (not 3,4,5)");
    }
  }
  void setDrData(initializers::LTSTree *userDrTree,
                 initializers::DynamicRupture *drDescription) {
    drTree = userDrTree;
    dynRup = drDescription;
  }

  void init(const std::unordered_map<std::string, double*>& FaultParams);
  void initFaceToLtsMap();

  void writePickpointOutput(double time, double dt);
  bool isAtPickpoint(double time, double dt);
  void updateElementwiseOutput();


  virtual void tiePointers(seissol::initializers::Layer&  layerData,
                           seissol::initializers::DynamicRupture *dynRup,
                           seissol::Interoperability &e_interoperability) {

    real  (*mu)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->mu);
    real  (*slip)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip);
    real  (*slip1)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip1);
    real  (*slip2)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slip2);
    real  (*slipRate1)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRate1);
    real  (*slipRate2)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->slipRate2);
    real  (*rupture_time)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->rupture_time);
    real  (*peakSR)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->peakSR);
    real  (*tracXY)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->tracXY);
    real  (*tracXZ)[ init::QInterpolated::Stop[0] ] = layerData.var(dynRup->tracXZ);



    DRFaceInformation* faceInformation = layerData.var(dynRup->faceInformation);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      e_interoperability.copyFrictionOutputToFortran(ltsFace,  meshFace, mu, slip, slip1, slip2, slipRate1, slipRate2, rupture_time, peakSR, tracXY, tracXZ);
    }
  }

  virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) = 0;

protected:
  void initEwOutput(const std::unordered_map<std::string, double*>& FaultParams);
  void initPickpointOutput();

  std::string constructPpReveiverFileName(int receiverGlobalIndex) const;
  void calcFaultOutput(const OutputType type, OutputData& state, double time = 0.0);

  GeneralParamsT generalParams;

  std::unique_ptr<ElementWiseOutput> ewOutputBuilder{nullptr};
  OutputData ewOutputData{};
  //std::vector<std::pair<initializers::Layer*, size_t>> faceToLtsMap{};

  std::unique_ptr<PickpointOutput>  ppOutputBuilder{nullptr};
  OutputData ppOutputState{};

  initializers::LTSTree *drTree{nullptr};
  seissol::initializers::DynamicRupture *dynRup{nullptr};

  std::vector<std::pair<initializers::Layer*, size_t>> faceToLtsMap{};
  size_t iterationStep{0};
};

#endif //SEISSOL_DROUTOUT_DRBASE_HPP
