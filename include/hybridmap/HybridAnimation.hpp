//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "NeutralAtomArchitecture.hpp"
#include "NeutralAtomDefinitions.hpp"
#include "operations/Operation.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace qc {
class AnimationAtoms {
  using axesId   = std::uint32_t;
  using marginId = std::uint32_t;

protected:
  uint32_t                  colorSlm    = 0;
  uint32_t                  colorAod    = 1;
  uint32_t                  colorLocal  = 2;
  [[maybe_unused]] uint32_t colorGlobal = 3;
  uint32_t                  colorCz     = 4;

  std::map<CoordIndex, HwQubit>        coordIdxToId;
  std::map<HwQubit, std::pair<fp, fp>> idToCoord;
  std::map<HwQubit, uint32_t>          axesIds;
  std::map<HwQubit, uint32_t>          marginIds;
  uint32_t                             axesIdCounter   = 0;
  uint32_t                             marginIdCounter = 0;

  axesId   addAxis(HwQubit id);
  void     removeAxis(HwQubit id) { axesIds.erase(id); }
  marginId addMargin(HwQubit id);
  void     removeMargin(HwQubit id) { marginIds.erase(id); }

public:
  AnimationAtoms(const std::map<HwQubit, HwQubit>& initHwPos,
                 const NeutralAtomArchitecture&    arch);

  std::string        getInitString();
  std::string        getEndString(fp endTime);
  static std::string createCsvLine(fp startTime, HwQubit id, fp x, fp y,
                                   uint32_t size = 1, uint32_t color = 0,
                                   bool axes = false, axesId axId = 0,
                                   bool margin = false, marginId marginId = 0,
                                   fp marginSize = 0);
  std::string createCsvOp(const std::unique_ptr<Operation>& op, fp startTime,
                          fp endTime, const qc::NeutralAtomArchitecture& arch);
};

} // namespace qc
