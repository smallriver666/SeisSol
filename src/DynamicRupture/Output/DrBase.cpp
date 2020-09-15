#include <unordered_map>
#include "Solver/Interoperability.h"
#include "Initializer/tree/Layer.hpp"
#include "Initializer/DynamicRupture.h"
#include "SeisSol.h"
#include "DrBase.hpp"
#include "ResultWriter/common.hpp"


std::vector<double> getAllVertices(const seissol::dr::ReceiverPointsT& ReceiverPoints) {
  std::vector<double> vertices(3 * (3  * ReceiverPoints.size()), 0.0);

  for (size_t pointIndex{0}; pointIndex < ReceiverPoints.size(); ++pointIndex) {
    for (int vertexIndex{0}; vertexIndex < 3; ++vertexIndex) {
      const size_t globalVertexIndex = 3 * pointIndex + vertexIndex;

      vertices[3 * globalVertexIndex] = ReceiverPoints[pointIndex].GlobalSubTet.Points[vertexIndex].x;
      vertices[3 * globalVertexIndex + 1] = ReceiverPoints[pointIndex].GlobalSubTet.Points[vertexIndex].y;
      vertices[3 * globalVertexIndex + 2] = ReceiverPoints[pointIndex].GlobalSubTet.Points[vertexIndex].z;
    }
  }
  return vertices;
}


std::vector<unsigned int> getCellConnectivity(const seissol::dr::ReceiverPointsT& ReceiverPoints) {
  std::vector<unsigned int> cells(3 * ReceiverPoints.size());

  for (size_t pointIndex{0}; pointIndex < ReceiverPoints.size(); ++pointIndex) {
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
  if (m_EwOutput) {
    m_EwOutput->init(FaultParams);

    const auto& ReceiverPoints = m_EwOutput->m_ReceiverPoints;
    auto cellConnectivity = getCellConnectivity(ReceiverPoints);
    auto vertices = getAllVertices(ReceiverPoints);
    auto intMask = convertMaskFromBoolToInt(m_EwOutput->m_ElementwiseParams.OutputMask);

    std::string NewFaultFilePrefix = m_GeneralParams.OutputFilePrefix.data() + std::string("-new");
    double printTime = m_EwOutput->m_ElementwiseParams.PrintTimeIntervalSec;
    auto backendType = seissol::writer::backendType(m_GeneralParams.XdmfWriterBackend.data());

    seissol::SeisSol::main.secondFaultWriter().init(cellConnectivity.data(),
                                                    vertices.data(),
                                                    static_cast<unsigned int>(ReceiverPoints.size()),
                                                    static_cast<unsigned int>(3 * ReceiverPoints.size()),
                                                    &intMask[0],
                                                    const_cast<const real**>(m_EwOutput->m_TmpState),
                                                    NewFaultFilePrefix.data(),
                                                    printTime,
                                                    backendType);

    for (size_t i = 0; i < ReceiverPoints.size(); ++i) {
      m_EwOutput->m_TmpState[5][i] = 1.0;
    }

    seissol::SeisSol::main.secondFaultWriter().write(1.0);
  }

  if (m_PpOutput) {
    m_PpOutput->init();
  }
}