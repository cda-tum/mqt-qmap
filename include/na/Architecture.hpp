//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"

#include <map>
#include <string>

namespace na {

enum class Scope { Global, Local };

enum class Type { SLM, AOD };
static std::map<std::string, Type> stringToType = { {"SLM", Type::SLM}, {"AOD", Type::AOD} };
Type getTypeOfString(std::string& s) {
  auto it = stringToType.find(s);
  return it->second;
}

class Point {
public:
  const std::uint16_t x;
  const std::uint16_t y;
  Point(std::uint16_t x, std::uint16_t y) : x(x), y(y){};
  inline Point operator-(Point p) { return Point(x - p.x, y - p.y); }
  inline Point operator+(Point p) { return Point(x + p.x, y + p.y); }
  [[nodiscard]] auto length() { return std::abs(x * x + y * y); }
};

// template <typename Value> using Grid = std::vector<std::vector<Value>>;
using Index  = size_t;
using Zone   = Index;
using Site   = std::tuple<Point, Zone, Type>;
using Value  = qc::fp;
using Number = std::uint16_t;

class Architecture {
public:
  /**
   * @brief Class to store the decoherence times of a neutral atom
   * architecture
   * @details
   * The decoherence times of a neutral atom architecture are:
   * - T1
   * - T2
   * - effective decoherence time
   */
  class DecoherenceTimes {
  public:
    Value const t1;
    Value const t2;
    Value const tEff;
    inline DecoherenceTimes(Value t1, Value t2)
        : t1(t1), t2(t2), tEff(t1 * t2 / (t1 + t2)) {}
    operator const double() const { return tEff; }
  };

protected:
  std::string              name;
  std::vector<std::string> zones;
  std::vector<Site>        sites;
  Number                   nSites;
  DecoherenceTimes         decoherenceTimes;
  Number                   nAods;
  Value                    minAtomDistance;
  Value                    interactionRadius;

public:
  Architecture(std::string& filename);
  virtual ~Architecture() = default;
  [[nodiscard]] inline auto getName() { return name; }
  [[nodiscard]] inline auto getNZones() { return zones.size(); }
  [[nodiscard]] inline auto getZoneLabel(const Index& i) { return zones.at(i); }
  [[nodiscard]] inline auto getNSites() { return nSites; }
  [[nodiscard]] inline auto getType(const Index& i) {
    return std::get<Type>(sites.at(i));
  }
  [[nodiscard]] inline auto getZone(const Index& i) {
    return std::get<Zone>(sites.at(i));
  }
  [[nodiscard]] inline auto getPos(const Index& i) {
    return std::get<Point>(sites.at(i));
  }
  [[nodiscard]] inline auto getDecoherenceTimes() { return decoherenceTimes; }
  [[nodiscard]] inline auto getNAods() { return nAods; }
  [[nodiscard]] inline auto getMinAtomDistance() { return minAtomDistance; }
  [[nodiscard]] inline auto getInteractionRadius() { return interactionRadius; }
  [[nodiscard]] static Architecture fromJSON(const std::string& filename);
  [[nodiscard]] static Architecture fromJSON(std::ifstream& fs);

  [[nodiscard]] auto getDistance(Index i, Index j) {
    return (getPos(j) - getPos(i)).length();
  }
};
} // namespace na
