/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#ifndef CS_SYNTHESISDATA_HPP
#define CS_SYNTHESISDATA_HPP

#include "LogicBlock/LogicBlock.hpp"
namespace cs {
    struct SynthesisData {
        std::uint32_t                                            nqubits;
        std::uint32_t                                            timesteps;
        const std::set<std::pair<std::uint16_t, std::uint16_t>>& reducedCM;
        const std::vector<std::uint16_t>&                        qubitChoice;
        std::unique_ptr<logicbase::LogicBlock>&                  lb;
        const logicbase::LogicMatrix&                            x;
        const logicbase::LogicMatrix&                            z;
        const logicbase::LogicVector&                            r;
        const logicbase::LogicMatrix3D&                          gS;
        const logicbase::LogicMatrix3D&                          gTwoQubit;
    };
} // namespace cs

#endif //CS_SYNTHESISDATA_HPP
