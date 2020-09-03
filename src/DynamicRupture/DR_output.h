#ifndef SEISSOL_DR_OUTPUT_H
#define SEISSOL_DR_OUTPUT_H

#include "Initializer/InputAux.hpp"

namespace seissol {
    namespace dr {
        namespace output {
            class Base;
            class FL_2;
            class FL_3;
            class FL_33;

            enum class OutputType : int {None = 0,
                                         AtPickpoint = 3,
                                         Elementwise = 4,
                                         AtPickpointAndElementwise = 5};
        }
    }
}

class seissol::dr::output::Base{
public:
    virtual ~Base() {}

    void setInputParam(const YAML::Node& InputParams) {
      using namespace initializers;

      if (!InputParams["dynamicrupture"]) {
        throw std::runtime_error("dynamic rupture params. is not provided");
      }

      const YAML::Node& DrParams = InputParams["dynamicrupture"];

      m_Params.OutputPointType = static_cast<OutputType>(getParamIfExists(DrParams, "outputpointtype", 3));
      m_Params.SlipRateOutputType = getParamIfExists(DrParams, "sliprateoutputtype", 1);
      m_Params.FrictionLawType = getParamIfExists(DrParams, "fl", 0);
      m_Params.BackgroundType = getParamIfExists(DrParams, "backgroundtype", 0);
      m_Params.IsRfOutputOn = getParamIfExists(DrParams, "rf_output_on", false);
      m_Params.IsDsOutputOn = getParamIfExists(DrParams, "ds_output_on", false);
      m_Params.IsMagnitudeOutputOn = getParamIfExists(DrParams, "magnitude_output_on", false);
      m_Params.IsEnergyRateOutputOn = getParamIfExists(DrParams, "energy_rate_output_on", false);
      m_Params.IsGpWiseOutput = getParamIfExists(DrParams, "gpwise", false);
      m_Params.IsTermalPressureOn = getParamIfExists(DrParams, "thermalpress", false);
      m_Params.BackgroundType = getParamIfExists(DrParams, "energy_rate_printtimeinterval", 1);

      m_Params.IsRfTimeOn = m_Params.IsRfOutputOn;

      if (m_Params.IsDsOutputOn && !m_Params.IsRfOutputOn) {
        m_Params.IsRfOutputOn = true;
        m_Params.IsRfTimeOn = true;
      }

      switch (m_Params.OutputPointType) {
        case OutputType::None:
          break;

        case OutputType::AtPickpoint:
          // readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
          initPickpointParams(InputParams);
          break;

        case OutputType::Elementwise:
          // readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
          initElementWiseParams(InputParams);
          break;

        case OutputType::AtPickpointAndElementwise:
          // readpar_faultElementwise(EQN,BND,IC,DISC,IO,CalledFromStructCode)
          // readpar_faultAtPickpoint(EQN,BND,IC,DISC,IO,CalledFromStructCode)
          initPickpointParams(InputParams);
          initElementWiseParams(InputParams);
          break;

        default:
          throw std::runtime_error("Unkown fault output type (e.g.3,4,5)");
      }
    }

    void initPickpointParams(const YAML::Node& InputParams) {
      using namespace initializers;

      if (!InputParams["pickpoint"]) {
        throw std::runtime_error("pickpoint output parameters for dynamic rupture is not provided");
      }

      const YAML::Node& Params = InputParams["pickpoint"];
      m_PickpointParams.PrintTimeInterval = getParamIfExists(Params, "printtimeinterval", 1);
      m_PickpointParams.NumOutputPoints = getParamIfExists(Params, "noutpoints", 0);
      m_PickpointParams.PPFileName = getParamIfExists(Params, "ppfilename", std::string());

      if (Params["outputmask"]) {
        convertStringToMask(Params["outputmask"].as<std::string>(), m_PickpointParams.OutputMask);
      }
    }

    void initElementWiseParams(const YAML::Node& InputParams) {
      using namespace initializers;
      if (!InputParams["elementwise"]) {
        throw std::runtime_error("elementwise fault output parameters for dynamic rupture is not provided");
      }

      const YAML::Node& Params = InputParams["elementwise"];

      m_ElementwiseParams.PrintTimeInterval = getParamIfExists(Params, "printtimeinterval", 2);
      m_ElementwiseParams.PrintTimeIntervalSec = getParamIfExists(Params, "printtimeinterval_sec", 1.0);
      m_ElementwiseParams.PrintIntervalCriterion = getParamIfExists(Params, "printintervalcriterion", 1);
      m_ElementwiseParams.RefinementStrategy = getParamIfExists(Params, "refinement_strategy", 2);
      m_ElementwiseParams.Refinement = getParamIfExists(Params, "refinement", 2);

      if (Params["outputmask"]) {
        convertStringToMask(Params["outputmask"].as<std::string>(), m_ElementwiseParams.OutputMask);
      }
    }

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
  struct ParamsT {
    OutputType OutputPointType{OutputType::AtPickpoint};
    int SlipRateOutputType{1};
    int FrictionLawType{0};
    int BackgroundType{0};
    bool IsRfOutputOn{false};
    bool IsDsOutputOn{false};
    bool IsMagnitudeOutputOn{false};
    bool IsEnergyRateOutputOn{false};
    bool IsGpWiseOutput{false};
    bool IsTermalPressureOn{false};
    int EnergyRatePrintTimeInterval{1};
    bool IsRfTimeOn{false};
  } m_Params;

  struct PickpointParamsT {
    std::array<bool, 12> OutputMask{true, true, true}; // the rest is false by default
    int PrintTimeInterval{1};
    int NumOutputPoints{0};
    std::string PPFileName{};

  } m_PickpointParams;

  struct ElementwiseFaultParamsT {
    int PrintTimeInterval{2};
    real PrintTimeIntervalSec{1.0};
    int PrintIntervalCriterion{1};
    std::array<bool, 12> OutputMask{true, true, true, true};
    int RefinementStrategy{2};
    int Refinement{2};
  } m_ElementwiseParams;
};

class seissol::dr::output::FL_2 : public seissol::dr::output::Base {
public:
    virtual void tiePointers(seissol::initializers::Layer&  layerData,
            seissol::initializers::DynamicRupture *dynRup,
            seissol::Interoperability &e_interoperability) override {

        Base::tiePointers(layerData, dynRup, e_interoperability);

        seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);

        DRFaceInformation*                    faceInformation = layerData.var(ConcreteLts->faceInformation);
        real  *averaged_Slip = layerData.var(ConcreteLts->averaged_Slip);
        real  (*dynStress_time)[ init::QInterpolated::Stop[0] ] = layerData.var(ConcreteLts->dynStress_time);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
            unsigned meshFace = static_cast<int>(faceInformation[ltsFace].meshFace);
            e_interoperability.copyFrictionOutputToFortranFL2(ltsFace,  meshFace,
                    averaged_Slip,  dynStress_time);
        }
    }

    virtual void postCompute(seissol::initializers::DynamicRupture &DynRup) override {
        std::cout << "output vars for FL_2\n";
    }
};

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



class seissol::dr::output::FL_33 : public seissol::dr::output::Base {
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

#endif //SEISSOL_DR_OUTPUT_H
