//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "Configuration.hpp"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>

using json = nlohmann::json;

namespace na {

Configuration::Configuration(const std::string& filename) {
  std::ifstream fs(filename);
  if (!fs.good()) {
    std::stringstream ss;
    ss << "Could not open JSON file " << filename << ".";
    throw std::runtime_error(ss.str());
  }
  *this = Configuration(fs);
}

Configuration::Configuration(std::istream& fs) {
  json data;
  try {
    fs >> data;
  } catch (const std::exception& e) {
    std::stringstream ss;
    ss << "Error parsing JSON: " << e.what() << std::endl;
    throw std::runtime_error(ss.str());
  }

  if (data.contains("patch") && data["patch"].is_object()) {
    auto& patch = data["patch"];
    if (patch.contains("rows") && patch["rows"].is_number_unsigned()) {
      patchRows = patch["rows"];
    }
    if (patch.contains("cols") && patch["cols"].is_number_unsigned()) {
      patchCols = patch["cols"];
    }
  }
  if (data.contains("method") && data["method"].is_string()) {
    method = getMethodOfString(data["method"]);
  }
}
} // namespace na
