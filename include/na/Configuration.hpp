#pragma once

#include <string>

namespace na {
class Configuration {
private:
    int patchRows = 1;
    int patchCols = 1;
    bool singleQubitScheduling = false;

public:
    explicit Configuration(const std::string& filename);
    explicit Configuration(std::istream& fs);
    virtual ~Configuration() = default;

    [[nodiscard]] auto getPatchRows() const -> int { return patchRows; }
    [[nodiscard]] auto getPatchCols() const -> int { return patchCols; }
    [[nodiscard]] auto isSingleQubitSchedulingAllowed() const -> bool { return singleQubitScheduling; }
};
} // namespace na
