/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#pragma once

#include <iostream>

enum class Encoding { Naive, Commander, Bimander };

static std::string toString(const Encoding encoding) {
  switch (encoding) {
  case Encoding::Naive:
    return "naive";
  case Encoding::Commander:
    return "commander";
  case Encoding::Bimander:
    return "bimander";
  }
  return " ";
}

[[maybe_unused]] static Encoding
encodingFromString(const std::string& encoding) {
  if (encoding == "naive" || encoding == "0") {
    return Encoding::Naive;
  } else if (encoding == "commander" || encoding == "1") {
    return Encoding::Commander;
  } else if (encoding == "bimander" || encoding == "2") {
    return Encoding::Bimander;
  } else {
    throw std::invalid_argument("Invalid encoding value: " + encoding);
  }
}
