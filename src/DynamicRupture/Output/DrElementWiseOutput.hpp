#ifndef SEISSOL_DRELEMENTWISEOUTPUT_HPP
#define SEISSOL_DRELEMENTWISEOUTPUT_HPP

#include "DrDataTypes.hpp"
#include "FaultRefiner/FaultRefiner.hpp"
#include "Geometry/MeshReader.h"

namespace seissol {
  namespace dr {
    namespace output {
      class ElementWiseOutput;
      struct FaultGeomParamsT;
    }
  }
}

struct seissol::dr::output::FaultGeomParamsT{
  int NumSubTriangles{};
  int NumSubElements{};
  size_t NumSides{};
};

class seissol::dr::output::ElementWiseOutput {
public:
  void setParams(const ElementwiseFaultParamsT& Params, const MeshReader* Reader) {
    m_ElementwiseParams = Params;
    m_FaultRefiner = getRefiner(Params.RefinementStrategy);
    m_MeshReader = Reader;
  }


  void init() {
    m_GeomParam.NumSides = m_MeshReader->getFault().size();
    m_GeomParam.NumSubTriangles = m_FaultRefiner->getNumSubTriangles();
    m_GeomParam.NumSubElements = std::pow(m_GeomParam.NumSubTriangles, m_ElementwiseParams.Refinement);
    //m_Points.resize(m_GeomParam.TotalNumReceivers);

    const auto& Fault = m_MeshReader->getFault();

    // get arrays of elements and vertices from the mesher
    const auto& ElementsInfo = m_MeshReader->getElements();
    const auto& VerticesInfo = m_MeshReader->getVertices();


    size_t TotalNumReceivers = 0;
    // iterate through each fault side
    for (size_t i = 0; i < m_GeomParam.NumSides; ++i) {

      // get a Global Element ID for the current fault face
      auto ElementIndex = Fault[i].element;
      if (ElementIndex > 0) {

        // store coords of vertices of the current ELEMENT
        std::array<const double *, 4> ElementVerticesCoords{};
        for (int ElementVertexId = 0; ElementVertexId < 4; ++ElementVertexId) {
          auto GlobalVertexId = ElementsInfo[ElementIndex].vertices[ElementVertexId];
          ElementVerticesCoords[ElementVertexId] = VerticesInfo[GlobalVertexId].coords;
        }

        // compute rupture coordinates inside of the reference element
        std::array<real, 3> Xi{};
        std::array<real, 3> Eta{};
        std::array<real, 3> Zeta{};
        auto LocalFaceSideId = Fault[i].side;
        std::tie(Xi, Eta, Zeta) = getXiEtaZeta(LocalFaceSideId);

        // store coords of vertices of the current FACE
        std::array<const double *, 3> FaultFaceVerticesCoords{};

        for (int FaceVertexId = 0; FaceVertexId < 3; ++FaceVertexId) {
          auto ElementVertexId = getElementVertexId(LocalFaceSideId, FaceVertexId);
          auto GlobalVertexId = ElementsInfo[ElementIndex].vertices[ElementVertexId];
          FaultFaceVerticesCoords[FaceVertexId] = VerticesInfo[GlobalVertexId].coords;
        }

        TotalNumReceivers += m_GeomParam.NumSubElements;
        std::cout << ElementIndex << " | " << LocalFaceSideId << std::endl;
      }
    }

    // TODO: check ElementVerticesCoords and FaultFaceVerticesCoords
    // TODO: refineFaultOutput

    std::cout << "Hello world\n";


    // TODO: eval_faultreceiver
    // TODO: create_fault_rotationmatrix
  }

private:
  static int getElementVertexId(int LocalSideId, int LocalFaceVertexId) {
    // 4 - number of faces of an element
    // 3 - number of vertices of a face
    static int LocalVertexMap[4][3] = {{0, 2, 1},  // Local tet. vertices of tet. side I
                                       {0, 1, 3},  // Local tet. vertices of tet. side II
                                       {0, 3, 2},  // Local tet. vertices of tet. side III
                                       {1, 2, 3}}; // Local tet. vertices of tet. side IV
    return LocalVertexMap[LocalSideId][LocalFaceVertexId];
  }

  std::tuple<std::array<real, 3>, std::array<real, 3>, std::array<real, 3>> getXiEtaZeta(int LocalSideId) {
    std::array<real, 3> Xi{};
    std::array<real, 3> Eta{};
    std::array<real, 3> Zeta{};
    switch (LocalSideId) {
      case 0:
        Xi = {0.0, 0.0, 1.0};
        Eta = {0.0, 1.0, 0.0};
        Zeta = {0.0, 0.0, 0.0};
        break;
      case 1:
        Xi = {0.0, 1.0, 0.0};
        Eta = {0.0, 0.0, 0.0};
        Zeta = {0.0, 0.0, 1.0};
        break;
      case 2:
        Xi = {0.0, 0.0, 0.0};
        Eta = {0.0, 0.0, 1.0};
        Zeta = {0.0, 1.0, 0.0};
        break;
      case 3:
        Xi = {1.0, 0.0, 0.0};
        Eta = {0.0, 1.0, 0.0};
        Zeta = {0.0, 0.0, 1.0};
        break;
      default:
        throw std::runtime_error("Unknown Local Side Id. Must be 0, 1, 2 or 3");
    }

    return std::make_tuple(Xi, Eta, Zeta);
  }


  ElementwiseFaultParamsT m_ElementwiseParams;
  FaultGeomParamsT m_GeomParam;
  std::unique_ptr<FaultRefinerInterface> m_FaultRefiner{nullptr};

  PointsT m_Points{};
  const MeshReader* m_MeshReader;
};

#endif //SEISSOL_DRELEMENTWISEOUTPUT_HPP
