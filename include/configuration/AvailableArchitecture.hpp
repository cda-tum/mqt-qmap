/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_AVAILABLEARCHITECTURE_HPP
#define QMAP_AVAILABLEARCHITECTURE_HPP

#include <iostream>
#include <map>

enum class AvailableArchitecture {
    IBM_QX4,
    IBM_QX5,
    IBMQ_Yorktown,
    IBMQ_London,
    IBMQ_Bogota,
    IBMQ_Tokyo,
    IBMQ_Ehningen
};

[[maybe_unused]] static std::string toString(const AvailableArchitecture architecture) {
    switch (architecture) {
        case AvailableArchitecture::IBM_QX4:
            return "IBM_QX4";
        case AvailableArchitecture::IBM_QX5:
            return "IBM_QX5";
        case AvailableArchitecture::IBMQ_Yorktown:
            return "IBMQ_Yorktown";
        case AvailableArchitecture::IBMQ_London:
            return "IBMQ_London";
        case AvailableArchitecture::IBMQ_Bogota:
            return "IBMQ_Bogota";
        case AvailableArchitecture::IBMQ_Tokyo:
            return "IBMQ_Tokyo";
        case AvailableArchitecture::IBMQ_Ehningen:
            return "IBMQ_Ehningen";
    }
    return " ";
}

[[maybe_unused]] static AvailableArchitecture architectureFromString(const std::string& architecture) {
    if (architecture == "IBM_QX4" || architecture == "0") {
        return AvailableArchitecture::IBM_QX4;
    } else if (architecture == "IBM_QX5" || architecture == "1") {
        return AvailableArchitecture::IBM_QX5;
    } else if (architecture == "IBMQ_Yorktown" || architecture == "2") {
        return AvailableArchitecture::IBMQ_Yorktown;
    } else if (architecture == "IBMQ_London" || architecture == "3") {
        return AvailableArchitecture::IBMQ_London;
    } else if (architecture == "IBMQ_Bogota" || architecture == "4") {
        return AvailableArchitecture::IBMQ_Bogota;
    } else if (architecture == "IBMQ_Tokyo" || architecture == "5") {
        return AvailableArchitecture::IBMQ_Tokyo;
    } else if (architecture == "IBMQ_Ehningen" || architecture == "6") {
        return AvailableArchitecture::IBMQ_Ehningen;
    } else {
        throw std::invalid_argument("Invalid architecture value: " + architecture);
    }
}

[[maybe_unused]] static std::string getCouplingMapSpecification(AvailableArchitecture architecture) {
    static std::map<AvailableArchitecture, std::string> architectureMap{
            {AvailableArchitecture::IBM_QX4, "5\n1 0\n2 0\n2 1\n3 2\n3 4\n2 4"},
            {AvailableArchitecture::IBM_QX5, "16\n1 0\n15 0\n1 2\n2 3\n15 2\n3 4\n3 14\n5 4\n13 4\n6 5\n12 5\n6 7\n6 11\n8 7\n7 10\n9 8\n9 10\n11 10\n12 11\n12 13\n13 14\n15 14"},
            {AvailableArchitecture::IBMQ_Yorktown, "5\n0 1\n1 0\n0 2\n2 0\n1 2\n2 1\n2 3\n3 2\n3 4\n4 3\n2 4\n4 2"},
            {AvailableArchitecture::IBMQ_London, "5\n0 1\n1 0\n1 2\n2 1\n1 3\n3 1\n3 4\n4 3"},
            {AvailableArchitecture::IBMQ_Bogota, "5\n0 1\n1 0\n1 2\n2 1\n2 3\n3 2\n3 4\n4 3"},
            {AvailableArchitecture::IBMQ_Tokyo, "20\n0 1\n1 0\n1 2\n2 1\n2 3\n3 2\n3 4\n4 3\n5 6\n6 5\n6 7\n7 6\n7 8\n8 7\n8 9\n9 8\n10 11\n11 10\n11 12\n12 11\n12 13\n13 12\n13 14\n14 13\n15 16\n16 15\n16 17\n17 16\n17 18\n18 17\n18 19\n19 18\n0 5\n5 0\n5 10\n10 5\n10 15\n15 10\n1 6\n6 1\n6 11\n11 6\n11 16\n16 11\n2 7\n7 2\n7 12\n12 7\n12 17\n17 12\n3 8\n8 3\n8 13\n13 8\n13 18\n18 13\n4 9\n9 4\n9 14\n14 9\n14 19\n19 14\n5 11\n11 5\n11 17\n17 11\n1 7\n7 1\n7 13\n13 7\n13 9\n9 13\n3 9\n9 3\n2 6\n6 2\n6 10\n10 6\n4 8\n8 4\n8 12\n12 8\n12 16\n16 12\n14 18\n18 14"},
            {AvailableArchitecture::IBMQ_Ehningen, "27\n0 1\n1 0\n1 2\n1 4\n2 1\n2 3\n3 2\n3 5\n4 1\n4 7\n5 3\n5 8\n6 7\n7 4\n7 6\n7 10\n8 5\n8 9\n8 11\n9 8\n10 7\n10 12\n11 8\n11 14\n12 10\n12 13\n12 15\n13 12\n13 14\n14 11\n14 13\n14 16\n15 12\n15 18\n16 14\n16 19\n17 18\n18 15\n18 17\n18 21\n19 16\n19 20\n19 22\n20 19\n21 18\n21 23\n22 19\n22 25\n23 21\n23 24\n24 23\n24 25\n25 22\n25 24\n25 26\n26 25"}
    };
    return architectureMap.at(architecture);
}

#endif //QMAP_AVAILABLEARCHITECTURE_HPP
