/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#ifndef CS_SYNTHESISDATA_HPP
#define CS_SYNTHESISDATA_HPP

#include "LogicBlock/LogicBlock.hpp"
#include "utils.hpp"

namespace cs {
struct SynthesisData {
  std::uint16_t                           nqubits;
  std::size_t                             timesteps;
  const CouplingMap&                      reducedCM;
  const QubitSubset&                      qubitChoice;
  std::unique_ptr<logicbase::LogicBlock>& lb;
  const logicbase::LogicMatrix&           x;
  const logicbase::LogicMatrix&           z;
  const logicbase::LogicVector&           r;
  const logicbase::LogicMatrix3D&         gS;
  const logicbase::LogicMatrix3D&         gTwoQubit;
};
} // namespace cs

#endif // CS_SYNTHESISDATA_HPP
