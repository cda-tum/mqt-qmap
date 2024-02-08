#include "Configuration.hpp"
#include <iostream>
#include <fstream>
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
        if (patch.contains("rows") && patch["rows"].is_number_integer()) {
            patchRows = patch["rows"];
        }
        if (patch.contains("cols") && patch["cols"].is_number_integer()) {
            patchCols = patch["cols"];
        }
    }

    if (data.contains("singleQubitScheduling") && data["singleQubitScheduling"].is_binary()) {
        singleQubitScheduling = data["singleQubitScheduling"];
    }
};
} // namespace na
