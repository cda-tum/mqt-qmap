//
// Created by Sarah on 27.10.2022.
//

#ifndef QMAP_GATES_HPP
#define QMAP_GATES_HPP

#include "QuantumComputation.hpp"
namespace cs {

    class Gates {
    public:
        enum GATES {
            NOP,
            H,
            S,
            X,
            Y,
            Z,
            Sdag,
            CX,
        };

        static std::string gateName(GATES gate) {
            switch (gate) {
                case NOP:
                    return "NOP";
                case H:
                    return "H";
                case S:
                    return "S";
                case X:
                    return "X";
                case Y:
                    return "Y";
                case Z:
                    return "Z";
                case Sdag:
                    return "Sdag";
                case CX:
                    return "CX";
                default:
                    return "";
            }
        }

        static int toIndex(GATES gate) {
            switch (gate) {
                case NOP:
                    return 0;
                case H:
                    return 1;
                case S:
                    return 2;
                case X:
                    return 3;
                case Y:
                    return 4;
                case Z:
                    return 5;
                case Sdag:
                    return 6;
                case CX:
                    return 7;
                default:
                    return -1;
            }
        }

        static qc::OpType toOpType(GATES gate) {
            switch (gate) {
                case NOP:
                    return qc::OpType::None;
                case H:
                    return qc::OpType::H;
                case S:
                    return qc::OpType::S;
                case X:
                    return qc::OpType::X;
                case Y:
                    return qc::OpType::Y;
                case Z:
                    return qc::OpType::Z;
                case Sdag:
                    return qc::OpType::Sdag;
                case CX:
                    return qc::OpType::X;
                default:
                    return qc::OpType::None;
            }
        }

        static constexpr auto SINGLE_QUBIT = std::array<GATES, 7>{
                GATES::NOP,
                GATES::H,
                GATES::S,
                GATES::X,
                GATES::Y,
                GATES::Z,
                GATES::Sdag};

        static constexpr auto TWO_QUBIT = std::array<GATES, 3>{
                GATES::CX};

        static constexpr auto SINGLE_QUBIT_WITHOUT_NOP = std::array<GATES, 6>{
                GATES::H,
                GATES::S,
                GATES::X,
                GATES::Y,
                GATES::Z,
                GATES::Sdag};
    };
}

#endif //QMAP_GATES_HPP
