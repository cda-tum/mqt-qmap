/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/GateEncoding.hpp"
void cs::GateEncoding::makeSingleGateEncoding(const SynthesisData& data) {
    auto changes = logicbase::LogicTerm(true);
    // CONSISTENCY
    makeSingleGateConsistency(data);

    // GATE CONSTRAINTS
    for (std::size_t gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        auto aIt = data.qubitChoice.cbegin();
        for (std::size_t a = 0; a < data.nqubits; ++a) {
            // NO GATE
            changes = logicbase::LogicTerm(true);
            makeIConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes && (data.r[gateStep] == data.r[gateStep - 1]);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][0][a], changes);
            data.lb->assertFormula(changes);

            // H
            changes = logicbase::LogicTerm(true);
            makeHConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & data.z[gateStep - 1][a])));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][1][a], changes);

            data.lb->assertFormula(changes);

            // S
            changes = logicbase::LogicTerm(true);
            makeSConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & data.z[gateStep - 1][a])));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][2][a], changes);
            data.lb->assertFormula(changes);

            // Sdag
            changes = logicbase::LogicTerm(true);
            makeSdagConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & (data.x[gateStep - 1][a] ^ data.z[gateStep - 1][a]))));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            data.lb->assertFormula(changes);

            // X
            changes = logicbase::LogicTerm(true);
            makeXConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ data.z[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            data.lb->assertFormula(changes);

            // Y
            changes = logicbase::LogicTerm(true);
            makeYConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ (data.z[gateStep - 1][a]) ^ data.x[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            data.lb->assertFormula(changes);

            // Z
            changes = logicbase::LogicTerm(true);
            makeZConstraints(data, a, gateStep, changes);
            makeNotChangingSet(data, a, gateStep, changes);

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ data.x[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            data.lb->assertFormula(changes);

            auto bIt = data.qubitChoice.cbegin();
            for (std::size_t b = 0; b < data.nqubits; ++b) {
                const auto q0 = *aIt;
                const auto q1 = *bIt;
                if (data.reducedCM.find({q0, q1}) == data.reducedCM.end()) {
                    data.lb->assertFormula(!data.gTwoQubit[gateStep][a][b]);
                } else {
                    // CNOT
                    changes =
                            (data.r[gateStep] == (data.r[gateStep - 1] ^
                                                  ((data.x[gateStep - 1][a] & data.z[gateStep - 1][b]) &
                                                   ((data.x[gateStep - 1][b] ^ data.z[gateStep - 1][a]) ^
                                                    logicbase::LogicTerm((1 << data.nqubits) - 1, data.nqubits)))));
                    makeCNOTConstraints(data, a, b, gateStep, changes);
                    for (std::size_t c = 0; c < data.nqubits; ++c) { // All other entries do not change
                        if (a == c || b == c) {
                            continue;
                        }
                        changes = changes && (data.x[gateStep][c] == data.x[gateStep - 1][c]);
                        changes = changes && (data.z[gateStep][c] == data.z[gateStep - 1][c]);
                    }

                    changes = logicbase::LogicTerm::implies(data.gTwoQubit[gateStep][a][b], changes);
                    data.lb->assertFormula(changes);
                }
                ++bIt;
            }
            ++aIt;
        }
    }
}
void cs::GateEncoding::makeMultiGateEncoding(const SynthesisData& data) {
    auto changes = logicbase::LogicTerm(true);
    // CONSISTENCY
    makeMultiGateConsistency(data);

    // GATE CONSTRAINTS
    for (std::size_t gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        logicbase::LogicTerm rChanges = data.r[gateStep - 1];
        auto                 aIt      = data.qubitChoice.cbegin();
        for (std::size_t a = 0; a < data.nqubits; ++a) {
            // NO GATE
            changes = logicbase::LogicTerm(true);
            makeIConstraints(data, a, gateStep, changes);

            changes = logicbase::LogicTerm::implies(data.gS[gateStep][0][a], changes);
            data.lb->assertFormula(changes);

            // H
            changes = logicbase::LogicTerm(true);
            makeHConstraints(data, a, gateStep, changes);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][1][a],
                    rChanges ^ (data.x[gateStep - 1][a] & data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][1][a], changes);

            data.lb->assertFormula(changes);

            // S
            changes = logicbase::LogicTerm(true);
            makeSConstraints(data, a, gateStep, changes);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][2][a],
                    rChanges ^ (data.x[gateStep - 1][a] & data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][2][a], changes);
            data.lb->assertFormula(changes);

            // Sdag
            changes = logicbase::LogicTerm(true);
            makeSdagConstraints(data, a, gateStep, changes);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a],
                    rChanges ^ (data.x[gateStep - 1][a] & (data.x[gateStep - 1][a] ^ data.z[gateStep - 1][a])), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            data.lb->assertFormula(changes);

            // X
            changes = logicbase::LogicTerm(true);
            makeXConstraints(data, a, gateStep, changes);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a],
                    rChanges ^ (data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            data.lb->assertFormula(changes);

            // Y
            changes = logicbase::LogicTerm(true);
            makeYConstraints(data, a, gateStep, changes);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a],
                    rChanges ^ (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            data.lb->assertFormula(changes);

            // Z
            changes = logicbase::LogicTerm(true);
            makeZConstraints(data, a, gateStep, changes);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a],
                    rChanges ^ (data.x[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            data.lb->assertFormula(changes);

            // CNOT
            auto bIt = data.qubitChoice.cbegin();
            for (std::size_t b = 0; b < data.nqubits; ++b) {
                const auto q0 = *aIt;
                const auto q1 = *bIt;
                if (data.reducedCM.find({q0, q1}) == data.reducedCM.end()) {
                    data.lb->assertFormula(!data.gTwoQubit[gateStep][a][b]);
                } else {
                    changes  = logicbase::LogicTerm(true);
                    rChanges = logicbase::LogicTerm::ite(
                            data.gTwoQubit[gateStep][a][b],
                            (rChanges ^ ((data.x[gateStep - 1][a] & data.z[gateStep - 1][b]) &
                                         ((data.x[gateStep - 1][b] ^ data.z[gateStep - 1][a]) ^
                                          logicbase::LogicTerm((1 << data.nqubits) - 1, data.nqubits)))),
                            rChanges);
                    makeCNOTConstraints(data, a, b, gateStep, changes);
                    changes = logicbase::LogicTerm::implies(data.gTwoQubit[gateStep][a][b], changes);
                    data.lb->assertFormula(changes);
                }
                ++bIt;
            }
            ++aIt;
        }
        data.lb->assertFormula(data.r[gateStep] == rChanges);
    }
}
void cs::GateEncoding::makeGateEncoding(const SynthesisData& data, const cs::Configuration& configuration) {
    switch (configuration.target) {
        case cs::TargetMetric::DEPTH:
        case cs::TargetMetric::FIDELITY:
            makeMultiGateEncoding(data);
            break;
        default:
            makeSingleGateEncoding(data);
            break;
    }
}
void cs::GateEncoding::makeNotChangingSet(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    for (std::size_t b = 0; b < data.nqubits; ++b) {
        if (a == b) {
            continue;
        }
        changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
        changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
    }
}
void cs::GateEncoding::makeIConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
    changes = changes && (data.z[gateStep][a] == data.z[gateStep - 1][a]);
}
void cs::GateEncoding::makeHConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes = changes && (data.z[gateStep][a] == data.x[gateStep - 1][a]);
    changes = changes && (data.x[gateStep][a] == data.z[gateStep - 1][a]);
}
void cs::GateEncoding::makeSConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes = changes && (data.z[gateStep][a] ==
                          (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
}
void cs::GateEncoding::makeSdagConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes =
            (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
}
void cs::GateEncoding::makeXConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes =
            (data.z[gateStep][a] == data.z[gateStep - 1][a]);
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
}
void cs::GateEncoding::makeYConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes =
            (data.z[gateStep][a] == data.z[gateStep - 1][a]);
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
}
void cs::GateEncoding::makeZConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes =
            (data.z[gateStep][a] == data.z[gateStep - 1][a]);
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
}
void cs::GateEncoding::makeSingleGateConsistency(const cs::SynthesisData& data) {
    logicbase::LogicTerm changes = logicbase::LogicTerm(true);
    // One gate per qubit, per step
    for (std::size_t gateStep = 1U; gateStep < data.timesteps + 1U; ++gateStep) {
        std::vector<logicbase::LogicTerm> vars{};
        auto                              aIt = data.qubitChoice.cbegin();
        for (std::size_t a = 0U; a < data.nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(data.gS[gateStep][Gates::toIndex(gate)][a]);
            }
            auto bIt = data.qubitChoice.cbegin();
            for (std::size_t b = 0; b < data.nqubits; ++b) {
                const auto q0 = *aIt;
                const auto q1 = *bIt;
                ++bIt;
                if (a == b || data.reducedCM.find({q0, q1}) == data.reducedCM.end()) {
                    continue;
                }
                vars.emplace_back(data.gTwoQubit[gateStep][a][b]);
            }
            ++aIt;
        }
        data.lb->assertFormula(encodings::exactlyOneCmdr(
                encodings::groupVars(vars, static_cast<std::size_t>(vars.size() / 2U)),
                logicbase::LogicTerm::noneTerm(), data.lb.get()));
    }
}
void cs::GateEncoding::makeMultiGateConsistency(const cs::SynthesisData& data) {
    auto changes = logicbase::LogicTerm(true);
    // One gate per qubit, per step
    for (std::size_t gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        for (std::size_t a = 0; a < data.nqubits; ++a) {
            std::vector<logicbase::LogicTerm> vars{};
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(data.gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (std::size_t b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                vars.emplace_back(data.gTwoQubit[gateStep][a][b]);
                vars.emplace_back(data.gTwoQubit[gateStep][b][a]);
            }
            data.lb->assertFormula(encodings::exactlyOneCmdr(
                    encodings::groupVars(vars, static_cast<std::size_t>(vars.size() / 2)),
                    logicbase::LogicTerm::noneTerm(), data.lb.get()));
        }
    }
    // Maximum any combination of 1 and 2 qubit gates adding up to n
    for (std::size_t gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        changes = logicbase::LogicTerm(0);
        for (std::size_t a = 0; a < data.nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                changes = changes + data.gS[gateStep][Gates::toIndex(gate)][a];
            }
            for (std::size_t b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes + data.gTwoQubit[gateStep][a][b] + data.gTwoQubit[gateStep][a][b];
            }
        }
        changes = changes < logicbase::LogicTerm(static_cast<int>(data.nqubits + 1));
        data.lb->assertFormula(changes);
    }
}
void cs::GateEncoding::makeCNOTConstraints(const cs::SynthesisData& data, std::size_t a, std::size_t b, std::size_t gateStep, encodings::LogicTerm& changes) {
    changes = changes && (data.x[gateStep][b] ==
                          (data.x[gateStep - 1][b] ^ data.x[gateStep - 1][a]));
    changes = changes && (data.z[gateStep][a] ==
                          (data.z[gateStep - 1][a] ^ data.z[gateStep - 1][b]));
    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
    changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
}
