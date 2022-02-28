/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#ifndef QMAP_METHOD_HPP
#define QMAP_METHOD_HPP

#include <iostream>

enum class Method {
    None,
    Exact,
    Heuristic
};

static std::string toString(const Method method) {
    switch (method) {
        case Method::None:
            return "none";
        case Method::Exact:
            return "exact";
        case Method::Heuristic:
            return "heuristic";
    }
    return " ";
}

[[maybe_unused]] static Method methodFromString(const std::string& method) {
    if (method == "none" || method == "0") {
        return Method::None;
    } else if (method == "exact" || method == "1") {
        return Method::Exact;
    } else if (method == "heuristic" || method == "2") {
        return Method::Heuristic;
    } else {
        throw std::invalid_argument("Invalid method value: " + method);
    }
}

#endif //QMAP_METHOD_HPP
