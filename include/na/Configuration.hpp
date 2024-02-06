#pragma once

#include <string>

class Configuration {
public:
    int patchRows;
    int patchCols;
    std::string singleQubitScheduling;

    // Constructor
    Configuration();

    // Parse configuration from a JSON file
    bool parseFromFile(const std::string& filename);
};
