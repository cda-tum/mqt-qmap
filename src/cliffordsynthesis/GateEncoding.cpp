/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
*/

#include "cliffordsynthesis/GateEncoding.hpp"
void cs::GateEncoding::makeSingleGateEncoding(const SynthesisData& data) {
    logicbase::LogicTerm changes = logicbase::LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        std::vector<logicbase::LogicTerm> vars{};
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(data.gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (int b = 0; b < data.nqubits; ++b) {
                if (a == b || data.reducedCM.find({data.qubitChoice.at(a), data.qubitChoice.at(b)}) ==
                                      data.reducedCM.end()) {
                    continue;
                }
                vars.emplace_back(data.gC[gateStep][a][b]);
            }
        }
        data.lb->assertFormula(encodings::exactlyOneCmdr(
                encodings::groupVars(vars, static_cast<std::size_t>(vars.size() / 2U)),
                logicbase::LogicTerm::noneTerm(), data.lb.get()));
    }

    // GATE CONSTRAINTS
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            // NO GATE
            changes = (data.x[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.z[gateStep][a] == data.z[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes && (data.r[gateStep] == data.r[gateStep - 1]);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][0][a], changes);
            data.lb->assertFormula(changes);

            // H
            changes = (data.z[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.z[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & data.z[gateStep - 1][a])));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][1][a], changes);

            data.lb->assertFormula(changes);

            // S
            changes =
                    (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & data.z[gateStep - 1][a])));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][2][a], changes);
            data.lb->assertFormula(changes);

            // Sdag
            changes =
                    (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^
                                            (data.x[gateStep - 1][a] & (data.x[gateStep - 1][a] ^ data.z[gateStep - 1][a]))));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            data.lb->assertFormula(changes);

            // Z
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ data.x[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            data.lb->assertFormula(changes);

            // X
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ data.z[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            data.lb->assertFormula(changes);

            // Y
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes && (data.x[gateStep][b] == data.x[gateStep - 1][b]);
                changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
            }

            changes = changes &&
                      (data.r[gateStep] == (data.r[gateStep - 1] ^ (data.z[gateStep - 1][a]) ^ data.x[gateStep - 1][a]));
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            data.lb->assertFormula(changes);

            // CNOT
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (data.reducedCM.find({data.qubitChoice.at(a), data.qubitChoice.at(b)}) ==
                    data.reducedCM.end()) {
                    data.lb->assertFormula(!data.gC[gateStep][a][b]);
                } else {
                    changes =
                            (data.r[gateStep] == (data.r[gateStep - 1] ^
                                                  ((data.x[gateStep - 1][a] & data.z[gateStep - 1][b]) &
                                                   ((data.x[gateStep - 1][b] ^ data.z[gateStep - 1][a]) ^
                                                    logicbase::LogicTerm((1 << data.nqubits) - 1, data.nqubits)))));
                    changes = changes && (data.x[gateStep][b] ==
                                          (data.x[gateStep - 1][b] ^ data.x[gateStep - 1][a]));
                    changes = changes && (data.z[gateStep][a] ==
                                          (data.z[gateStep - 1][a] ^ data.z[gateStep - 1][b]));
                    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
                    changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);

                    for (unsigned int c = 0; c < data.nqubits; ++c) { // All other entries do not change
                        if (a == c || b == c) {
                            continue;
                        }
                        changes = changes && (data.x[gateStep][c] == data.x[gateStep - 1][c]);
                        changes = changes && (data.z[gateStep][c] == data.z[gateStep - 1][c]);
                    }

                    changes = logicbase::LogicTerm::implies(data.gC[gateStep][a][b], changes);
                    data.lb->assertFormula(changes);
                }
            }
        }
    }
}
void cs::GateEncoding::makeMultiGateEncoding(const SynthesisData& data) {
    logicbase::LogicTerm changes = logicbase::LogicTerm(true);
    // CONSISTENCY
    // One gate per qubit, per step
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            std::vector<logicbase::LogicTerm> vars{};
            for (auto gate: Gates::SINGLE_QUBIT) {
                vars.emplace_back(data.gS[gateStep][Gates::toIndex(gate)][a]);
            }
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                vars.emplace_back(data.gC[gateStep][a][b]);
                vars.emplace_back(data.gC[gateStep][b][a]);
            }
            data.lb->assertFormula(encodings::exactlyOneCmdr(
                    encodings::groupVars(vars, static_cast<std::size_t>(vars.size() / 2)),
                    logicbase::LogicTerm::noneTerm(), data.lb.get()));
        }
    }
    // Maximum any combination of 1 and 2 qubit gates adding up to n
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        changes = logicbase::LogicTerm(0);
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            for (auto gate: Gates::SINGLE_QUBIT) {
                changes = changes + data.gS[gateStep][Gates::toIndex(gate)][a];
            }
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (a == b) {
                    continue;
                }
                changes = changes + data.gC[gateStep][a][b] + data.gC[gateStep][a][b];
            }
        }
        changes = changes < logicbase::LogicTerm(static_cast<int>(data.nqubits + 1));
        data.lb->assertFormula(changes);
    }

    // GATE CONSTRAINTS
    for (unsigned int gateStep = 1; gateStep < data.timesteps + 1; ++gateStep) {
        logicbase::LogicTerm rChanges = data.r[gateStep - 1];
        for (unsigned int a = 0; a < data.nqubits; ++a) {
            // NO GATE
            changes = logicbase::LogicTerm(true);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][0][a], changes);
            data.lb->assertFormula(changes);

            // H
            changes = logicbase::LogicTerm(true);
            changes = changes && (data.z[gateStep][a] == data.x[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.z[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][1][a],
                    rChanges ^ (data.x[gateStep - 1][a] & data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][1][a], changes);

            data.lb->assertFormula(changes);

            // S
            changes = logicbase::LogicTerm(true);
            changes = changes && (data.z[gateStep][a] ==
                                  (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][2][a],
                    rChanges ^ (data.x[gateStep - 1][a] & data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][2][a], changes);
            data.lb->assertFormula(changes);

            // Sdag
            changes =
                    (data.z[gateStep][a] == (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]));
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a],
                    rChanges ^ (data.x[gateStep - 1][a] & (data.x[gateStep - 1][a] ^ data.z[gateStep - 1][a])), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Sdag)][a], changes);
            data.lb->assertFormula(changes);

            // Z
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a],
                    rChanges ^ (data.x[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Z)][a], changes);
            data.lb->assertFormula(changes);

            // X
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a],
                    rChanges ^ (data.z[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::X)][a], changes);
            data.lb->assertFormula(changes);

            // Y
            changes =
                    (data.z[gateStep][a] == data.z[gateStep - 1][a]);
            changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);

            rChanges = logicbase::LogicTerm::ite(
                    data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a],
                    rChanges ^ (data.z[gateStep - 1][a] ^ data.x[gateStep - 1][a]), rChanges);
            changes = logicbase::LogicTerm::implies(data.gS[gateStep][Gates::toIndex(Gates::GATES::Y)][a], changes);
            data.lb->assertFormula(changes);

            // CNOT
            for (unsigned int b = 0; b < data.nqubits; ++b) {
                if (data.reducedCM.find({data.qubitChoice.at(a), data.qubitChoice.at(b)}) ==
                    data.reducedCM.end()) {
                    data.lb->assertFormula(!data.gC[gateStep][a][b]);
                } else {
                    changes  = logicbase::LogicTerm(true);
                    rChanges = logicbase::LogicTerm::ite(
                            data.gC[gateStep][a][b],
                            (rChanges ^ ((data.x[gateStep - 1][a] & data.z[gateStep - 1][b]) &
                                         ((data.x[gateStep - 1][b] ^ data.z[gateStep - 1][a]) ^
                                          logicbase::LogicTerm((1 << data.nqubits) - 1, data.nqubits)))),
                            rChanges);
                    changes = changes && (data.x[gateStep][b] ==
                                          (data.x[gateStep - 1][b] ^ data.x[gateStep - 1][a]));
                    changes = changes && (data.z[gateStep][a] ==
                                          (data.z[gateStep - 1][a] ^ data.z[gateStep - 1][b]));
                    changes = changes && (data.x[gateStep][a] == data.x[gateStep - 1][a]);
                    changes = changes && (data.z[gateStep][b] == data.z[gateStep - 1][b]);
                    changes = logicbase::LogicTerm::implies(data.gC[gateStep][a][b], changes);
                    data.lb->assertFormula(changes);
                }
            }
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
        case cs::TargetMetric::GATES:
        case cs::TargetMetric::TWO_QUBIT_GATES:
        default:
            makeSingleGateEncoding(data);
            break;
    }
}
