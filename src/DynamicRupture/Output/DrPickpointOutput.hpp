#ifndef SEISSOL_PICKPOINTOUTPUT_HPP
#define SEISSOL_PICKPOINTOUTPUT_HPP

#include "DrDataTypes.hpp"
#include "Geometry/MeshReader.h"

namespace seissol {
  namespace dr {
    namespace output {
      class PickpointOutput;
    }
  }
}

class seissol::dr::output::PickpointOutput {
public:
  void setParams(const PickpointParamsT& Params, const MeshReader* Reader) {
    m_PickpointParams = Params;
    m_meshReader = Reader;
  }

  void init() {
    readCoordsFromFile();
    // findElementsContainingPoints();
    // InitPointsIndices();
    // ProjectPointsToFaces();
    // FindClosestGpPoint();
  }

  void readCoordsFromFile() {
    using namespace initializers;
    StringsT Content = FileProcessor::getFileAsStrings(m_PickpointParams.PPFileName);
    FileProcessor::removeEmptyLines(Content);

    if (m_PickpointParams.NumOutputPoints > static_cast<int>(Content.size()))
      throw std::runtime_error("requested num. of fault pick-points is more than the file contains");

    // iterate line by line and initialize DrRecordPoints
    for (const auto &Line: Content) {
      std::array<real, 3> Coords{};
      convertStringToMask(Line, Coords);

      ReceiverPointT Point{};
      for (int i = 0; i < 3; ++i) {
        Point.Global.Coords[0] = Coords[0];
      }

      m_Points.push_back(Point);
    }
  }

private:
  PickpointParamsT m_PickpointParams;
  ReceiverPointsT m_Points{};
  const MeshReader* m_meshReader;
};

#endif //SEISSOL_PICKPOINTOUTPUT_HPP
