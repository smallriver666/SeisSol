#ifndef SEISSOL_DROUTOUT_DRBASE_HPP
#define SEISSOL_DROUTOUT_DRBASE_HPP

#include "Initializer/InputAux.hpp"
#include "DynamicRupture/Output/DrParametersInitializer.hpp"
#include "DynamicRupture/Output/DrElementWiseOutput.hpp"
#include "DynamicRupture/Output/DrPickpointOutput.hpp"
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
    m_GeneralParams = Reader.getDrGeneralParams();

    // adjust general output parameters
    m_GeneralParams.IsRfTimeOn = m_GeneralParams.IsRfOutputOn;
    if (m_GeneralParams.IsDsOutputOn && !m_GeneralParams.IsRfOutputOn) {
      m_GeneralParams.IsRfOutputOn = true;
      m_GeneralParams.IsRfTimeOn = true;
    }

    PickpointParamsT PpParams;
    ElementwiseFaultParamsT EwParams;
    switch (m_GeneralParams.OutputPointType) {
      case OutputType::None:
        break;

      case OutputType::AtPickpoint:
        // readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
        m_PpOutput.reset(new PickpointOutput);
        m_PpOutput->setParams(Reader.getPickPointParams(), &Mesher);
        break;

      case OutputType::Elementwise:
        // readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
        m_EwOutput.reset(new ElementWiseOutput);
        m_EwOutput->setParams(Reader.getElementwiseFaultParams(), &Mesher);
        break;

      case OutputType::AtPickpointAndElementwise:
        // readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
        // readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
        m_PpOutput.reset(new PickpointOutput);
        m_PpOutput->setParams(Reader.getPickPointParams(), &Mesher);

        m_EwOutput.reset(new ElementWiseOutput);
        m_EwOutput->setParams(Reader.getElementwiseFaultParams(), &Mesher);
        break;

      default:
        throw std::runtime_error("Unkown fault output type (e.g.3,4,5)");
    }
  }

  void init(const std::unordered_map<std::string, double*>& FaultParams);

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
  GeneralParamsT m_GeneralParams;

  std::unique_ptr<ElementWiseOutput> m_EwOutput{nullptr};
  std::unique_ptr<PickpointOutput>  m_PpOutput{nullptr};
};

#endif //SEISSOL_DROUTOUT_DRBASE_HPP
