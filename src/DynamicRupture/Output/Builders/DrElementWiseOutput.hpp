#ifndef SEISSOL_DRELEMENTWISEOUTPUT_HPP
#define SEISSOL_DRELEMENTWISEOUTPUT_HPP

#include "DynamicRupture/Output/DrDataTypes.hpp"
#include "DynamicRupture/Output/FaultRefiner/FaultRefiner.hpp"
#include "Geometry/MeshReader.h"
#include "DynamicRupture/Output/DrOutputAux.hpp"
#include "Parallel/MPI.h"

namespace seissol {
  namespace dr {
    namespace output {
      class Base;
      class ElementWiseOutput;
      struct FaultGeomParamsT;
    }
  }
}

struct seissol::dr::output::FaultGeomParamsT {
  int numSubTriangles{};
  int numSubElements{};
  size_t numSides{};
};

class seissol::dr::output::ElementWiseOutput {
public:
  friend Base;

  ~ElementWiseOutput() {
    auto deallocateVars = [](auto& var, int) {
      var.releaseData();
    };
    aux::forEach(drVars, deallocateVars);
  }

  void setParams(const ElementwiseFaultParamsT& params, const MeshReader* reader) {
    elementwiseParams = params;
    meshReader = reader;
    const int localRank = MPI::mpi.rank();
  }


  void init(const std::unordered_map<std::string, double*>& faultParams) {
    initReceiverLocations();
    assignNearestGaussianPoints(receiverPoints);
    initOutputVariables();

    /*
    allocateOutputVariables();
    initOutputVariables();
    initRotationMatrices(FaultParams);
    */
  }

  void initOutputVariables() {
    auto assignMask = [this](auto& var, int index) {
      var.isActive = this->elementwiseParams.outputMask[index];
    };
    aux::forEach(drVars, assignMask);

    auto allocateVariables = [this](auto& var, int) {
      var.size = var.isActive ? this->receiverPoints.size() : 0;
      if (var.isActive) {
        for (int dim = 0; dim < var.dim(); ++dim)
          var.data[dim] = new real[var.size];
      }
      else {
        for (int dim = 0; dim < var.dim(); ++dim)
          var.data[dim] = nullptr;
      }
    };
    aux::forEach(drVars, allocateVariables);

    real initialValue = 0.0;
    auto initVars = [initialValue](auto& var, int) {
      if (var.isActive) {
        for (int dim = 0; dim < var.dim(); ++dim) {
          for (size_t i = 0; i < var.size; ++i)
            var[dim][i] = initialValue;
        }
      }
    };
    aux::forEach(drVars, initVars);
  }

  void initReceiverLocations() {
    std::unique_ptr<FaultRefinerInterface> faultRefiner{nullptr};
    faultRefiner = getRefiner(elementwiseParams.refinementStrategy);

    geomParam.numSides = meshReader->getFault().size();
    geomParam.numSubTriangles = faultRefiner->getNumSubTriangles();
    geomParam.numSubElements = std::pow(geomParam.numSubTriangles, elementwiseParams.refinement);
    //m_Points.resize(geomParam.TotalNumReceivers);


    logInfo(localRank) << "CPP: Initialising Fault output. Refinement strategy: "
                    << elementwiseParams.refinementStrategy
                  << " Number of sub-triangles: "
                  << geomParam.numSubTriangles;


    // get arrays of elements and vertices from the mesher
    const auto &faultInfo = meshReader->getFault();
    const auto &elementsInfo = meshReader->getElements();
    const auto &verticesInfo = meshReader->getVertices();


    // iterate through each fault side
    for (size_t faceIndex = 0; faceIndex < geomParam.numSides; ++faceIndex) {

      // get a Global Element ID for the current fault face
      auto elementIndex = faultInfo[faceIndex].element;
      if (elementIndex > 0) {

        // store coords of vertices of the current ELEMENT
        std::array<const double *, 4> elementVerticesCoords{};
        for (int ElementVertexId = 0; ElementVertexId < 4; ++ElementVertexId) {
          auto globalVertexId = elementsInfo[elementIndex].vertices[ElementVertexId];
          elementVerticesCoords[ElementVertexId] = verticesInfo[globalVertexId].coords;
        }

        auto localFaceSideId = faultInfo[faceIndex].side;

        // init reference coordinates of the fault face
        ExtTriangle referenceFace = getReferenceFace(localFaceSideId);

        // init global coordinates of the fault face
        ExtTriangle globalFace = getReferenceFace(localFaceSideId);
        for (int faceVertexId = 0; faceVertexId < 3; ++faceVertexId) {
          auto elementVertexId = getElementVertexId(localFaceSideId, faceVertexId);
          auto globalVertexId = elementsInfo[elementIndex].vertices[elementVertexId];

          globalFace.points[faceVertexId].x = verticesInfo[globalVertexId].coords[0];
          globalFace.points[faceVertexId].y = verticesInfo[globalVertexId].coords[1];
          globalFace.points[faceVertexId].z = verticesInfo[globalVertexId].coords[2];
        }

        faultRefiner->refineAndAccumulate(elementwiseParams.refinement,
                                          faceIndex,
                                          localFaceSideId,
                                          referenceFace,
                                          globalFace);
      }
    }

    // retrieve all receivers from a fault face refiner
    receiverPoints = faultRefiner->moveAllReceiverPoints();
    faultRefiner.reset(nullptr);

    isDrPickOutput = !receiverPoints.empty();
    nDrPick = receiverPoints.size();
    nOutPoints = receiverPoints.size();
  }


  void allocateOutputVariables() {
    std::vector<real> currentPick(nDrPick);
    std::vector<real> tmpTime(elementwiseParams.maxPickStore);
    // std::vector<real> TmpState;
    // OutVal
    std::vector<real> rotationMatrix(nDrPick / geomParam.numSubTriangles, 0);
    std::vector<ConstantT> constant(nDrPick);

    // TODO: alloc DynRup_Constants
    // TODO: alloc DynRup_Constants_GlobInd

    // TODO: alloc CurrentPick
    // TODO: alloc TmpTime
    // TODO: alloc TmpState
    // TODO: alloc rotmat

  }

  void initRotationMatrices(const std::unordered_map<std::string, double*>& faultParams) {
    using namespace seissol::transformations;
    using RotationMatrixViewT = yateto::DenseTensorView<2, double, unsigned>;

    // allocate Rotation Matrices
    // Note: several receiver can share the same rotation matrix
    m_RotationMatrices.resize(geomParam.numSides);

    // init Rotation Matrices
    const auto &faultInfo = meshReader->getFault();
    for (size_t index = 0; index < geomParam.numSides; ++index) {
      const auto faceNormal = faultInfo[index].normal;
      VrtxCoords strike = {0.0, 0.0, 0.0};
      VrtxCoords dip = {0.0, 0.0, 0.0};

      computeStrikeAndDipVectors(faceNormal, strike, dip);

      std::vector<real> rotationMatrix(36, 0.0);
      RotationMatrixViewT rotationMatrixView(const_cast<real*>(rotationMatrix.data()), {6, 6});

      symmetricTensor2RotationMatrix(faceNormal, strike, dip, rotationMatrixView, 0, 0);
      m_RotationMatrices[index] = std::move(rotationMatrix);
    }
  }

  void initConstrains() {
    /*
    for (const auto& Point: m_ReceiverPoints) {
      const auto& RotationMatrix = m_RotationMatrices[Point.FaultFaceIndex];

    }
    */
  }


  void computeBasisFunctionsAtReceiver() {
    // TODO: compute self Basis functions
    // TODO: compute neighbor Basis functions
  }

  void evaluateInitialStressInFaultCS() {
    // Compute initialStressInFaultCS
  }

  /*
  void initOutputVariables() {
    // TODO: eval_faultreceiver
    // TODO: create_fault_rotationmatrix
  }
  */

private:

  ElementwiseFaultParamsT elementwiseParams;
  FaultGeomParamsT geomParam;

  ReceiverPointsT receiverPoints{};
  ConstantsT constants{};
  std::vector<int> outputLabels{};
  std::vector<std::vector<real>> m_RotationMatrices{};

  const MeshReader* meshReader;
  bool isDrPickOutput{};
  size_t nDrPick;
  size_t nOutPoints;
  int localRank{-1};
  //OutputState state;
  DrVarsT drVars;
};

#endif //SEISSOL_DRELEMENTWISEOUTPUT_HPP
