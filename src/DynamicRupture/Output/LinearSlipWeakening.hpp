#ifndef SEISSOL_DR_OUTPUT_LSW_HPP
#define SEISSOL_DR_OUTPUT_LSW_HPP

#include "DynamicRupture/Output/ReceiverBasedOutput.hpp"

namespace seissol::dr::output {
class LinearSlipWeakening : public ReceiverBasedOutput {
  public:
  void tiePointers(seissol::initializers::Layer& layerData,
                   seissol::initializers::DynamicRupture* drDescr,
                   seissol::Interoperability& eInteroperability) override {
    ReceiverBasedOutput::tiePointers(layerData, drDescr, eInteroperability);

    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(drDescr);

    DRFaceInformation* faceInformation = layerData.var(concreteLts->faceInformation);
    real* averagedSlip = layerData.var(concreteLts->averagedSlip);
    constexpr auto size = init::QInterpolated::Stop[0];
    real(*slipRate1)[size] = layerData.var(concreteLts->slipRate1);
    real(*slipRate2)[size] = layerData.var(concreteLts->slipRate2);
    real(*mu)[size] = layerData.var(concreteLts->mu);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
      eInteroperability.copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRate1, slipRate2, mu);
    }
  }

  protected:
  real computeLocalStrength() override {
    using DrLtsDescrT = seissol::initializers::LTS_LinearSlipWeakening;
    auto* cohesions = local.layer->var(static_cast<DrLtsDescrT*>(drDescr)->cohesion);
    auto cohesion = cohesions[local.ltsId][local.nearestGpIndex];

    auto resultingPressure = local.pressure + local.iniPressure - local.internalPressure;
    return -1.0 * local.frictionCoefficient * std::min(resultingPressure, static_cast<real>(0.0)) -
           cohesion;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_OUTPUT_LSW_HPP
