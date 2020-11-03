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
  int NumSubTriangles{};
  int NumSubElements{};
  size_t NumSides{};
};

class seissol::dr::output::ElementWiseOutput {
public:
  friend Base;

  ~ElementWiseOutput() {
    auto deallocateVars = [](auto& var, int) {
      if (var.isActive) {
        for (int dim = 0; dim < var.dim(); ++dim) {
          delete [] var.data[dim];
        }
      }
    };
    aux::forEach(drVars, deallocateVars);
  }

  void setParams(const ElementwiseFaultParamsT& Params, const MeshReader* Reader) {
    m_ElementwiseParams = Params;
    m_MeshReader = Reader;
    const int m_rank = MPI::mpi.rank();
  }


  void init(const std::unordered_map<std::string, double*>& FaultParams) {
    initReceiverLocations();
    //assignNearestGaussianPoints(m_ReceiverPoints);
    initOutputVariables();

    /*
    allocateOutputVariables();
    initOutputVariables();
    initRotationMatrices(FaultParams);
    */
  }

  void initOutputVariables() {
    auto assignMask = [this](auto& var, int index) {
      var.isActive = this->m_ElementwiseParams.OutputMask[index];
    };
    aux::forEach(drVars, assignMask);

    auto allocateVariables = [this](auto& var, int) {
      var.size = var.isActive ? this->m_ReceiverPoints.size() : 0;
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
    std::unique_ptr<FaultRefinerInterface> FaultRefiner{nullptr};
    FaultRefiner = getRefiner(m_ElementwiseParams.RefinementStrategy);

    m_GeomParam.NumSides = m_MeshReader->getFault().size();
    m_GeomParam.NumSubTriangles = FaultRefiner->getNumSubTriangles();
    m_GeomParam.NumSubElements = std::pow(m_GeomParam.NumSubTriangles, m_ElementwiseParams.Refinement);
    //m_Points.resize(m_GeomParam.TotalNumReceivers);


    logInfo(m_rank) << "CPP: Initialising Fault output. Refinement strategy: "
                  << m_ElementwiseParams.RefinementStrategy
                  << " Number of sub-triangles: "
                  << m_GeomParam.NumSubTriangles;


    // get arrays of elements and vertices from the mesher
    const auto &FaultInfo = m_MeshReader->getFault();
    const auto &ElementsInfo = m_MeshReader->getElements();
    const auto &VerticesInfo = m_MeshReader->getVertices();


    // iterate through each fault side
    for (size_t FaceIndex = 0; FaceIndex < m_GeomParam.NumSides; ++FaceIndex) {

      // get a Global Element ID for the current fault face
      auto ElementIndex = FaultInfo[FaceIndex].element;
      if (ElementIndex > 0) {

        // store coords of vertices of the current ELEMENT
        std::array<const double *, 4> ElementVerticesCoords{};
        for (int ElementVertexId = 0; ElementVertexId < 4; ++ElementVertexId) {
          auto GlobalVertexId = ElementsInfo[ElementIndex].vertices[ElementVertexId];
          ElementVerticesCoords[ElementVertexId] = VerticesInfo[GlobalVertexId].coords;
        }

        auto LocalFaceSideId = FaultInfo[FaceIndex].side;

        // init reference coordinates of the fault face
        ExtTriangle ReferenceFace = getReferenceFace(LocalFaceSideId);

        // init global coordinates of the fault face
        ExtTriangle GlobalFace = getReferenceFace(LocalFaceSideId);
        for (int FaceVertexId = 0; FaceVertexId < 3; ++FaceVertexId) {
          auto ElementVertexId = getElementVertexId(LocalFaceSideId, FaceVertexId);
          auto GlobalVertexId = ElementsInfo[ElementIndex].vertices[ElementVertexId];

          GlobalFace.Points[FaceVertexId].x = VerticesInfo[GlobalVertexId].coords[0];
          GlobalFace.Points[FaceVertexId].y = VerticesInfo[GlobalVertexId].coords[1];
          GlobalFace.Points[FaceVertexId].z = VerticesInfo[GlobalVertexId].coords[2];
        }

        FaultRefiner->refineAndAccumulate(m_ElementwiseParams.Refinement,
                                          FaceIndex,
                                          LocalFaceSideId,
                                          ReferenceFace,
                                          GlobalFace);
      }
    }

    // retrieve all receivers from a fault face refiner
    m_ReceiverPoints = FaultRefiner->moveAllReceiverPoints();
    FaultRefiner.reset(nullptr);

    DR_pick_output = !m_ReceiverPoints.empty();
    nDR_pick = m_ReceiverPoints.size();
    nOutPoints = m_ReceiverPoints.size();
  }


  void allocateOutputVariables() {
    std::vector<real> CurrentPick(nDR_pick);
    std::vector<real> TmpTime(m_ElementwiseParams.MaxPickStore);
    // std::vector<real> TmpState;
    // OutVal
    std::vector<real> RotationMatrix(nDR_pick / m_GeomParam.NumSubTriangles, 0);
    std::vector<ConstantT> Constanta(nDR_pick);

    // TODO: alloc DynRup_Constants
    // TODO: alloc DynRup_Constants_GlobInd

    // TODO: alloc CurrentPick
    // TODO: alloc TmpTime
    // TODO: alloc TmpState
    // TODO: alloc rotmat

  }

  void initRotationMatrices(const std::unordered_map<std::string, double*>& FaultParams) {
    using namespace seissol::transformations;
    using RotationMatrixViewT = yateto::DenseTensorView<2, double, unsigned>;

    // allocate Rotation Matrices
    // Note: several receiver can share the same rotation matrix
    m_RotationMatrices.resize(m_GeomParam.NumSides);

    // init Rotation Matrices
    const auto &FaultInfo = m_MeshReader->getFault();
    for (size_t Index = 0; Index < m_GeomParam.NumSides; ++Index) {
      const auto FaceNormal = FaultInfo[Index].normal;
      VrtxCoords Strike = {0.0, 0.0, 0.0};
      VrtxCoords Dip = {0.0, 0.0, 0.0};

      computeStrikeAndDipVectors(FaceNormal, Strike, Dip);

      std::vector<real> RotationMatrix(36, 0.0);
      RotationMatrixViewT RotationMatrixView(const_cast<real*>(RotationMatrix.data()), {6, 6});

      symmetricTensor2RotationMatrix(FaceNormal, Strike, Dip, RotationMatrixView, 0, 0);
      m_RotationMatrices[Index] = std::move(RotationMatrix);
    }

    /*
    // check whether the key exists. Otherwise, init with 0s
    auto IniBulk_xx = FaultParams.at("s_xx");
    auto IniBulk_yy = FaultParams.at("s_yy");
    auto IniBulk_zz = FaultParams.at("s_zz");
    auto IniShear_xy = FaultParams.at("s_xy");
    auto IniShear_yz = FaultParams.at("s_yz");
    auto IniShear_xz = FaultParams.at("s_xz");
    */
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
    // Compute InitialStressInFaultCS
  }

  /*
  void initOutputVariables() {
    // TODO: eval_faultreceiver
    // TODO: create_fault_rotationmatrix
  }
  */

private:

  ElementwiseFaultParamsT m_ElementwiseParams;
  FaultGeomParamsT m_GeomParam;

  ReceiverPointsT m_ReceiverPoints{};
  ConstantsT m_Constants{};
  std::vector<int> m_OutputLabels{};
  std::vector<std::vector<real>> m_RotationMatrices{};

  const MeshReader* m_MeshReader;

  bool DR_pick_output;
  size_t nDR_pick;
  size_t nOutPoints;
  int m_rank{-1};
  DrVarsT drVars;
};

#endif //SEISSOL_DRELEMENTWISEOUTPUT_HPP
