#ifndef SEISSOL_PICKPOINTOUTPUT_HPP
#define SEISSOL_PICKPOINTOUTPUT_HPP

#include "DynamicRupture/Output/DrDataTypes.hpp"
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
  void setParams(const PickpointParamsT& params, const MeshReader* reader) {
    pickpointParams = params;
    meshReader = reader;
  }

  void init() {
    readCoordsFromFile();
    // findElementsContainingPoints();
    // initPointsIndices();
    // projectPointsToFaces();
    // findClosestGpPoint();
  }

  void readCoordsFromFile() {
    using namespace initializers;
    StringsT content = FileProcessor::getFileAsStrings(pickpointParams.ppFileName);
    FileProcessor::removeEmptyLines(content);

    if (pickpointParams.numOutputPoints > static_cast<int>(content.size()))
      throw std::runtime_error("requested num. of fault pick-points is more than the file contains");

    // iterate line by line and initialize DrRecordPoints
    for (const auto &line: content) {
      std::array<real, 3> coords{};
      convertStringToMask(line, coords);

      ReceiverPointT point{};
      for (int i = 0; i < 3; ++i) {
        point.global.coords[0] = coords[0];
      }

      points.push_back(point);
    }
  }

private:
  PickpointParamsT pickpointParams;
  ReceiverPointsT points{};
  const MeshReader* meshReader;
};

#endif //SEISSOL_PICKPOINTOUTPUT_HPP
