#include <unordered_map>
#include "Solver/Interoperability.h"
#include "Initializer/tree/Layer.hpp"
#include "Initializer/DynamicRupture.h"
#include "SeisSol.h"
#include "DrBase.hpp"
#include "ResultWriter/common.hpp"


std::vector<double> getAllVertices(const seissol::dr::ReceiverPointsT& receiverPoints) {
  std::vector<double> vertices(3 * (3  * receiverPoints.size()), 0.0);

  for (size_t pointIndex{0}; pointIndex < receiverPoints.size(); ++pointIndex) {
    for (int vertexIndex{0}; vertexIndex < 3; ++vertexIndex) {
      const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;

      vertices[3 * globalVertexIndex] = receiverPoints[pointIndex].globalSubTet.points[vertexIndex].x;
      vertices[3 * globalVertexIndex + 1] = receiverPoints[pointIndex].globalSubTet.points[vertexIndex].y;
      vertices[3 * globalVertexIndex + 2] = receiverPoints[pointIndex].globalSubTet.points[vertexIndex].z;
    }
  }
  return vertices;
}


std::vector<unsigned int> getCellConnectivity(const seissol::dr::ReceiverPointsT& receiverPoints) {
  std::vector<unsigned int> cells(3 * receiverPoints.size());

  for (size_t pointIndex{0}; pointIndex < receiverPoints.size(); ++pointIndex) {
    for (int vertexIndex{0}; vertexIndex < 3; ++vertexIndex) {
      const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;
      cells[globalVertexIndex] = globalVertexIndex;
    }
  }
  return cells;
}


std::unique_ptr<int []> convertMaskFromBoolToInt(const std::array<bool, 12>& boolMask) {
  auto intMask = std::unique_ptr<int []>(new int[boolMask.size()]);

  for (size_t i = 0; i < boolMask.size(); ++i) {
    intMask[i] = static_cast<int>(boolMask[i]);
  }

  return intMask;
}


void seissol::dr::output::Base::init(const std::unordered_map<std::string, double*>& FaultParams) {
  if (ewOutputBuilder) {
    ewOutputBuilder->init(FaultParams);

    const auto& receiverPoints = ewOutputBuilder->receiverPoints;
    auto cellConnectivity = getCellConnectivity(receiverPoints);
    auto vertices = getAllVertices(receiverPoints);
    auto intMask = convertMaskFromBoolToInt(ewOutputBuilder->elementwiseParams.outputMask);

    std::string newFaultFilePrefix = generalParams.outputFilePrefix + std::string("-new");
    double printTime = ewOutputBuilder->elementwiseParams.printTimeIntervalSec;
    auto backendType = seissol::writer::backendType(generalParams.xdmfWriterBackend.c_str());

   std::vector<real*> dataPointers;
   auto recordPointers = [&dataPointers](auto& var, int) {
     if (var.isActive) {
       for (int dim = 0; dim < var.dim(); ++dim)
        dataPointers.push_back(var.data[dim]);
     }
   };
   aux::forEach(ewOutputBuilder->drVars, recordPointers);


   seissol::SeisSol::main.secondFaultWriter().init(cellConnectivity.data(),
                                                    vertices.data(),
                                                    static_cast<unsigned int>(receiverPoints.size()),
                                                    static_cast<unsigned int>(3 * receiverPoints.size()),
                                                    &intMask[0],
                                                    const_cast<const real**>(dataPointers.data()),
                                                    newFaultFilePrefix.data(),
                                                    printTime,
                                                    backendType);

    seissol::SeisSol::main.secondFaultWriter().setupCallbackObject(this);
  }

  if (ppOutputBuilder) {
    ppOutputBuilder->init();
  }
}

void seissol::dr::output::Base::initFaceToLtsMap() {

  faceToLtsMap.resize(drTree->getNumberOfCells(Ghost));
  for (auto it = drTree->beginLeaf(initializers::LayerMask(Ghost));
       it != drTree->endLeaf(); ++it) {

    DRFaceInformation* faceInformation = it->var(dynRup->faceInformation);
    for (size_t ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      faceToLtsMap[faceInformation[ltsFace].meshFace] = std::make_pair(&(*it), ltsFace);
    }
  }
}

void seissol::dr::output::Base::writePickpointOutput() {

}


void seissol::dr::output::Base::updateElementwiseOutput() {
  calcFaultOutput(ewOutputState);
}


using DrPaddedArrayT = real (*)[seissol::init::QInterpolated::Stop[0]];
void seissol::dr::output::Base::calcFaultOutput(const OutputState& state) {

  for (size_t i = 0; i < ewOutputBuilder->receiverPoints.size(); ++i) {
    auto ltsMap = faceToLtsMap[ewOutputBuilder->receiverPoints[i].faultFaceIndex];
    const auto layer = ltsMap.first;
    const auto ltsId = ltsMap.second;

    auto mu = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->mu));
    auto rt = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->rupture_time));
    auto slip = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->slip));
    auto peakSR = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->peakSR));

    const auto nearestGaussPoint = ewOutputBuilder->receiverPoints[i].nearestGpIndex;
    const auto faceIndex = ewOutputBuilder->receiverPoints[i].faultFaceIndex;

    const auto normal = ewOutputBuilder->faultDirections[faceIndex].faceNormal;
    const auto tangent1 = ewOutputBuilder->faultDirections[faceIndex].tangent1;
    const auto tangent2 = ewOutputBuilder->faultDirections[faceIndex].tangent2;
    const auto strike = ewOutputBuilder->faultDirections[faceIndex].strike;
    const auto dip = ewOutputBuilder->faultDirections[faceIndex].dip;





    auto &functionAndState = std::get<VariableID::FunctionAndState>(ewOutputBuilder->drVars);
    if (functionAndState.isActive) {
      functionAndState[ParamID::FUNCTION][i] = mu[ltsId][nearestGaussPoint];
    }


    auto &ruptureTime = std::get<VariableID::RuptureTime>(ewOutputBuilder->drVars);
    if (ruptureTime.isActive) {
      ruptureTime[0][i] = rt[ltsId][nearestGaussPoint];
    }


    auto &absoluteSlip = std::get<VariableID::AbsoluteSlip>(ewOutputBuilder->drVars);
    if (absoluteSlip.isActive) {
      absoluteSlip[0][i] = slip[ltsId][nearestGaussPoint];
    }


    auto &peakSlipsRate = std::get<VariableID::PeakSlipsRate>(ewOutputBuilder->drVars);
    if (peakSlipsRate.isActive) {
      peakSlipsRate[0][i] = peakSR[ltsId][nearestGaussPoint];
    }


    auto &slipVectors = std::get<VariableID::Slip>(ewOutputBuilder->drVars);
    if (slipVectors.isActive) {
      VrtxCoords crossProduct = {0.0, 0.0, 0.0};
      MeshTools::cross(strike, tangent1, crossProduct);

      double cos1 = MeshTools::dot(strike, tangent1);
      double scalarProd = MeshTools::dot(crossProduct, normal);

      // Note: cos1**2 can be greater than 1.0 because of rounding errors -> min
      double sin1 = std::sqrt(1.0 - std::min(1.0, cos1 * cos1));
      sin1 = (scalarProd > 0) ? sin1 : - sin1;

      auto slip1 = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->slip1));
      auto slip2 = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->slip2));

      slipVectors[DirectionID::STRIKE][i] = cos1 * slip1[ltsId][nearestGaussPoint] -
                                            sin1 * slip2[ltsId][nearestGaussPoint];

      slipVectors[DirectionID::DIP][i] = sin1 * slip1[ltsId][nearestGaussPoint] +
                                         cos1 * slip2[ltsId][nearestGaussPoint];
    }
  }
}