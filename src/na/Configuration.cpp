//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "Configuration.hpp"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace na {

Configuration::Configuration(const std::string& filename) {
  std::ifstream fs(filename);
  if (!fs.good()) {
    throw std::runtime_error("Could not open JSON file.");
  }
  *this = Configuration(fs);
}

Configuration::Configuration(std::istream& fs) {
  json data;
  try {
    fs >> data;
  } catch (const std::exception& e) {
    std::cerr << "Error parsing JSON: " << e.what() << std::endl;
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
};
} // namespace na
