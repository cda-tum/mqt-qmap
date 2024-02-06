#include "na/Configuration.hpp"
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// Constructor
Configuration::Configuration() : patchRows(0), patchCols(0) {}

// Parse configuration from a JSON file
bool Configuration::parseFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    json config;
    try {
        file >> config;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        return false;
    }

    if (config.contains("patch") && config["patch"].is_object()) {
        auto& patch = config["patch"];
        if (patch.contains("rows") && patch["rows"].is_number_integer()) {
            patchRows = patch["rows"];
        }
        if (patch.contains("cols") && patch["cols"].is_number_integer()) {
            patchCols = patch["cols"];
        }
    }

    if (config.contains("singleQubitScheduling") && config["singleQubitScheduling"].is_string()) {
        singleQubitScheduling = config["singleQubitScheduling"];
    }

    return true;
}
