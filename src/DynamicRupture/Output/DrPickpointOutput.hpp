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

  }

private:
  PickpointParamsT m_PickpointParams;
  PointsT m_Points{};
  const MeshReader* m_meshReader;
};

#endif //SEISSOL_PICKPOINTOUTPUT_HPP

/*
 *
    // read provided text file if it is needed
    if (PpParams.NumOutputPoints > 0) {

      StringsT Content = FileProcessor::getFileAsStrings(PpParams.PPFileName);
      FileProcessor::removeEmptyLines(Content);

      if (PpParams.NumOutputPoints > static_cast<int>(Content.size()))
        throw std::runtime_error("requested num. pickpoints is more that the file contains");

      // iterate line by line and initialize DrRecordPoints
      for (const auto &Line: Content) {
        std::array<real, 3> Coords{};
        convertStringToMask(Line, Coords);
        PpParams.DrRecordPoints.emplace_back(Coords);
      }
    }
 * */
