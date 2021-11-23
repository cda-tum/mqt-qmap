/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#ifndef QMAP_MAPPINGSETTINGS_HPP
#define QMAP_MAPPINGSETTINGS_HPP

#include "nlohmann/json.hpp"

/// Identity: q_i -> Q_i
/// Static: first layer is mapped q_c -> Q_c and q_t -> Q_t
/// Dynamic: Layout is generated on demand upon encountering a specific gate
enum class InitialLayoutStrategy {
    Identity,
    Static,
    Dynamic,
    None
};
static std::string toString(const InitialLayoutStrategy strategy) {
    switch (strategy) {
        case InitialLayoutStrategy::Identity:
            return "identity";
        case InitialLayoutStrategy::Static:
            return "static";
        case InitialLayoutStrategy::Dynamic:
            return "dynamic";
        case InitialLayoutStrategy::None:
            return "none";
    }
    return " ";
}
// map InitialLayoutStrategy values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM(InitialLayoutStrategy, {
                                                            {InitialLayoutStrategy::None, "none"},
                                                            {InitialLayoutStrategy::Identity, "identity"},
                                                            {InitialLayoutStrategy::Static, "static"},
                                                            {InitialLayoutStrategy::Dynamic, "dynamic"},
                                                    })

enum class LayeringStrategy {
    IndividualGates,
    DisjointQubits,
    OddGates,
    QubitTriangle,
    None
};

enum class Encodings {
    Naive,
    Commander,
    Bimander
};
enum class CMDRVariableGroupings {
    Fixed2,
    Fixed3,
    Halves,
    Logarithm
};
enum class SwapReductionStrategy {
    None,
    Custom,
    CouplingLimit,
    Increasing
};
static std::string toString(const LayeringStrategy strategy) {
    switch (strategy) {
        case LayeringStrategy::IndividualGates:
            return "individual_gates";
        case LayeringStrategy::DisjointQubits:
            return "disjoint_qubits";
        case LayeringStrategy::OddGates:
            return "odd_gates";
        case LayeringStrategy::QubitTriangle:
            return "qubit_triangle";
        case LayeringStrategy::None:
            return "none";
    }
    return " ";
}
// map LayeringStrategy values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM(LayeringStrategy, {
                                                       {LayeringStrategy::None, "none"},
                                                       {LayeringStrategy::IndividualGates, "individual_gates"},
                                                       {LayeringStrategy::DisjointQubits, "disjoint_qubits"},
                                                       {LayeringStrategy::OddGates, "odd_gates"},
                                                       {LayeringStrategy::QubitTriangle, "qubit_triangle"},
                                               })

static std::string toString(const Encodings encoding) {
    switch (encoding) {
        case Encodings::Naive:
            return "naive";
        case Encodings::Commander:
            return "commander";
        case Encodings::Bimander:
            return "bimander";
    }
    return " ";
}

NLOHMANN_JSON_SERIALIZE_ENUM(Encodings, {
                                                {Encodings::Naive, "naive"},
                                                {Encodings::Commander, "commander"},
                                                {Encodings::Bimander, "bimander"},
                                        })

static std::string toString(const CMDRVariableGroupings grouping) {
    switch (grouping) {
        case CMDRVariableGroupings::Fixed2:
            return "fixed2";
        case CMDRVariableGroupings::Fixed3:
            return "fixed3";
        case CMDRVariableGroupings::Logarithm:
            return "logarithm";
        case CMDRVariableGroupings::Halves:
            return "halves";
    }
    return " ";
}

NLOHMANN_JSON_SERIALIZE_ENUM(CMDRVariableGroupings, {
                                                            {CMDRVariableGroupings::Fixed2, "fixed2"},
                                                            {CMDRVariableGroupings::Fixed3, "fixed3"},
                                                            {CMDRVariableGroupings::Halves, "halves"},
                                                            {CMDRVariableGroupings::Logarithm, "logarithm"},
                                                    })

static std::string toString(const SwapReductionStrategy strategy) {
    switch (strategy) {
        case SwapReductionStrategy::CouplingLimit:
            return "coupling_limit";
        case SwapReductionStrategy::Custom:
            return "custom";
        case SwapReductionStrategy::None:
            return "none";
        case SwapReductionStrategy::Increasing:
            return "increasing";
    }
    return " ";
}

NLOHMANN_JSON_SERIALIZE_ENUM(SwapReductionStrategy, {{SwapReductionStrategy::None, "none"},
                                                     {SwapReductionStrategy::CouplingLimit, "coupling_limit"},
                                                     {SwapReductionStrategy::Custom, "custom"},
                                                     {SwapReductionStrategy::Increasing, "increasing"}})

struct MappingSettings {
    MappingSettings() = default;

    unsigned int timeout = 3600000; // 60min timeout
    void         setTimeout(unsigned int sec) { timeout = sec; }

    LayeringStrategy layeringStrategy = LayeringStrategy::None;

    /// Settings for heuristic approach
    InitialLayoutStrategy initialLayoutStrategy = InitialLayoutStrategy::None;

    bool admissibleHeuristic = true;
    bool verbose             = false;

    bool                  lookahead            = true;
    int                   nrLookaheads         = 15;
    int                   teleportationQubits  = 0;
    unsigned long long    teleportationSeed    = 0;
    bool                  teleportationFake    = false;
    double                firstLookaheadFactor = 0.75;
    double                lookaheadFactor      = 0.5;
    Encodings             encoding             = Encodings::Naive;
    CMDRVariableGroupings grouping             = CMDRVariableGroupings::Halves;
    bool                  enableLimits         = true;
    bool                  useBDD               = false;
    SwapReductionStrategy strategy             = SwapReductionStrategy::CouplingLimit;
    int                   limit                = 0;
    bool                  useQubitSubsets      = true;
};

#endif //QMAP_MAPPINGSETTINGS_HPP
