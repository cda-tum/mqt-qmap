//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Configuration.hpp"
#include "Definitions.hpp"
#include "na/NADefinitions.hpp"

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

/// The scope of an operation (Global or Local)
enum class Scope : uint8_t { Global, Local };
static const std::unordered_map<std::string, Scope> STRING_TO_SCOPE = {
    {"Global", Scope::Global},
    {"Local", Scope::Local},
    {"global", Scope::Global},
    {"local", Scope::Local}};

/**
 * @brief Get the Scope of a gate from a string
 *
 * @param s the name
 * @return Type
 */
inline Scope getScopeOfString(const std::string& s) {
  if (const auto it = STRING_TO_SCOPE.find(s); it != STRING_TO_SCOPE.end()) {
    return it->second;
  }
  std::stringstream ss;
  ss << "The scope " << s << " is not supported.";
  throw std::invalid_argument(ss.str());
}

/// For Indices
using Index = std::size_t;
/// The zones are just stored as int
using Zone = Index;
/// Any double-valued property
using Value = qc::fp;
/// Any information on numbers of something
using Number = std::int64_t;

class Architecture {
public:
  /**
   * @brief Struct to store the decoherence times of a neutral atom
   * architecture
   * @details
   * The decoherence times of a neutral atom architecture are:
   * - T1
   * - T2
   * - effective decoherence time
   */
  struct DecoherenceTimes {
    Value t1                                  = 0;
    Value t2                                  = 0;
    Value tEff                                = 0;
    DecoherenceTimes()                        = default;
    DecoherenceTimes(const DecoherenceTimes&) = default;
    virtual ~DecoherenceTimes()               = default;
    DecoherenceTimes(const Value t1, const Value t2)
        : t1(t1), t2(t2), tEff(t1 * t2 / (t1 + t2)) {}
    explicit operator double() const { return tEff; }
  };
  /**
   * @brief Strcut to store the properties of an operation.
   * @details Times are in µs, fidelities are in [0,1].
   */
  struct OperationProperties {
    Scope                    scope; // local or global
    std::unordered_set<Zone> zones; // the zones where the gate can be applied
    Value time;     // the time the gate takes to be applied in µs
    Value fidelity; // the fidelity of the gate
  };
  /**
   * @brief Struct to store the properties of a Shuttling operation (i.e. of the
   * AOD).
   * @details Times are in µs, fidelities are in [0,1], and velocities are in
   * µm/µs.
   */
  struct ShuttlingProperties {
    Number rows          = 0; // maximum number of rows in one AOD
    Number cols          = 0; // maximum number of columns in one AOD
    Number minX          = 0; // minimum x position of the AOD
    Number maxX          = 0; // maximum x position of the AOD
    Number minY          = 0; // minimum y position of the AOD
    Number maxY          = 0; // maximum y position of the AOD
    Value  speed         = 0; // speed of the AOD in µm/µs
    Value  fidelity      = 1; // fidelity during the shuttling
    Value  loadTime      = 0; // time to activate the AOD in µs
    Value  loadFidelity  = 1; // fidelity of the load
    Value  storeTime     = 0; // time to deactivate the AOD in µs
    Value  storeFidelity = 1; // fidelity of the store
  };
  struct ZoneProperties {
    std::string name;         // the name of the zone
    Number      minX     = 0; // minimum x dimension
    Number      maxX     = 0; // maximum x dimension
    Number      minY     = 0; // minimum y dimension
    Number      maxY     = 0; // maximum y dimension
    Value       fidelity = 1; // fidelity during idling
  };

protected:
  std::string name; // the name of the architecture
  std::vector<ZoneProperties>
      zones; // a mapping from zones (int) to their name from the config
  std::vector<Point> sites; // a vector of sites
  std::unordered_map<OpType, OperationProperties>
      gateSet; // all possible operations by their type, i.e. gate set
  DecoherenceTimes decoherenceTimes;          // the decoherence characteristic
  std::vector<ShuttlingProperties> shuttling; // all properties regarding AODs
  Index minAtomDistance = 0; // minimal distance that must be kept between atoms
  Index interactionRadius = 0; // the Rydberg radius
  Index noInteractionRadius =
      0; // sufficient radius to avoid Rydberg interaction
  std::vector<Zone> initialZones; // the zones where the atoms are initially

public:
  /**
   * @brief Import a new architecture from a file.
   *
   * @param jsonFn The path to the JSON file
   * @param csvFn The path to the CSV file
   */
  Architecture(const std::string& jsonFn, const std::string& csvFn);
  Architecture(std::istream& jsonS, std::istream& csvS);
  Architecture(const Architecture&)            = default;
  Architecture(Architecture&&)                 = default;
  virtual ~Architecture()                      = default;
  Architecture& operator=(const Architecture&) = default;
  Architecture& operator=(Architecture&&)      = default;

  [[nodiscard]] auto getName() const -> const std::string& { return name; }
  [[nodiscard]] auto getNZones() const -> Index { return zones.size(); }
  [[nodiscard]] auto getZoneLabel(const Index& i) const -> const std::string& {
    return zones[i].name;
  }
  [[nodiscard]] auto getInitialZones() const -> const std::vector<Zone>& {
    return initialZones;
  }
  [[nodiscard]] auto getNSites() const -> Index { return sites.size(); }
  [[nodiscard]] auto getPositionOfSite(const Index& i) const -> const Point& {
    return sites[i];
  }
  [[nodiscard]] auto getDecoherenceTimes() const -> const DecoherenceTimes& {
    return decoherenceTimes;
  }
  [[nodiscard]] auto getNShuttlingUnits() const -> Index {
    return shuttling.size();
  }
  [[nodiscard]] auto getPropertiesOfShuttlingUnit(const Index& i) const
      -> const ShuttlingProperties& {
    return shuttling[i];
  }
  [[nodiscard]] auto getMinAtomDistance() const -> Index {
    return minAtomDistance;
  }
  [[nodiscard]] auto getInteractionRadius() const -> Index {
    return interactionRadius;
  }
  [[nodiscard]] auto getNoInteractionRadius() const -> Index {
    return noInteractionRadius;
  }
  [[nodiscard]] auto
  getPropertiesOfZone(const Zone& zone) const -> const ZoneProperties& {
    return zones[zone];
  }
  [[nodiscard]] auto getPropertiesOfOperation(const OpType& t) const
      -> const OperationProperties& {
    if (auto it = gateSet.find(t); it != gateSet.end()) {
      return it->second;
    }
    std::stringstream ss;
    ss << "The operation " << t << " is not supported.";
    throw std::invalid_argument(ss.str());
  }
  /**
   * @brief Returns the distance between two sites.
   *
   * @param i address of first site
   * @param j address of second site
   * @return the distance in µm
   */
  [[nodiscard]] auto getDistance(const Index& i,
                                 const Index& j) const -> Index {
    return (getPositionOfSite(j) - getPositionOfSite(i)).length();
  }
  [[nodiscard]] auto getZoneAt(const Point& p) const -> Zone;
  [[nodiscard]] auto getZoneOfSite(const Index& i) const -> Zone {
    return getZoneAt(getPositionOfSite(i));
  }
  /// Checks whether the gate can be applied at all.
  [[nodiscard]] auto isAllowedLocally(const OpType& t) const -> bool;
  /// Checks whether the gate can be applied (locally) in this zone.
  [[nodiscard]] auto isAllowedLocally(const OpType& t,
                                      const Zone&   zone) const -> bool;
  /// Checks whether the gate can be applied (locally) on this qubit.
  [[nodiscard]] auto isAllowedLocallyAtSite(const OpType& t,
                                            const Index&  qubit) const -> bool;
  /// Checks whether the gate can be applied (locally) on this qubit.
  [[nodiscard]] auto isAllowedLocallyAt(const OpType& t,
                                        const Point&  p) const -> bool;
  /// Checks whether the gate is a global gate for this Zone.
  [[nodiscard]] auto isAllowedGlobally(const OpType& t) const -> bool;
  [[nodiscard]] auto isAllowedGlobally(const OpType& t,
                                       const Zone&   zone) const -> bool;
  [[nodiscard]] auto isInSameRow(const Index& i, const Index& j) const -> bool;
  [[nodiscard]] auto isInSameCol(const Index& i, const Index& j) const -> bool;
  [[nodiscard]] auto getNrowsInZone(const Zone& z) const -> Index;
  [[nodiscard]] auto getNcolsInZone(const Zone& z) const -> Index;
  [[nodiscard]] auto
  getSitesInRow(const Zone& z, const Index& row) const -> std::vector<Index>;
  [[nodiscard]] auto
  getSitesInCol(const Zone& z, const Index& col) const -> std::vector<Index>;
  [[nodiscard]] auto getRowInZoneOf(const Index& i) const -> Index;
  [[nodiscard]] auto getColInZoneOf(const Index& i) const -> Index;
  [[nodiscard]] auto getNearestXLeft(const Number& x, const Zone& z,
                                     bool proper = true) const -> Number;
  [[nodiscard]] auto getNearestXRight(const Number& x, const Zone& z,
                                      bool proper = true) const -> Number;
  [[nodiscard]] auto getNearestYUp(const Number& y, const Zone& z,
                                   bool proper = true) const -> Number;
  [[nodiscard]] auto getNearestYDown(const Number& y, const Zone& z,
                                     bool proper = true) const -> Number;
  [[nodiscard]] auto hasSiteLeft(const Point& p, bool proper = false,
                                        bool sameZone = false) const -> bool;
  [[nodiscard]] auto hasSiteRight(const Point& p, bool proper = false,
                                         bool sameZone = false) const -> bool;
  [[nodiscard]] auto hasSiteUp(const Point& p, bool proper = false,
                                      bool sameZone = false) const -> bool;
  [[nodiscard]] auto hasSiteDown(const Point& p, bool proper = false,
                                        bool sameZone = false) const -> bool;
  [[nodiscard]] auto getNearestSiteLeft(const Point& p, bool proper = false,
                                        bool sameZone = false) const -> Index;
  [[nodiscard]] auto getNearestSiteRight(const Point& p, bool proper = false,
                                         bool sameZone = false) const -> Index;
  [[nodiscard]] auto getNearestSiteUp(const Point& p, bool proper = false,
                                      bool sameZone = false) const -> Index;
  [[nodiscard]] auto getNearestSiteDown(const Point& p, bool proper = false,
                                        bool sameZone = false) const -> Index;
  [[nodiscard]] auto
                     getNearestSiteUpRight(const Point& p, bool proper = false,
                                           bool sameZone = false) const -> Index;
  [[nodiscard]] auto getNearestSiteUpLeft(const Point& p, bool proper = false,
                                          bool sameZone = false) const -> Index;
  [[nodiscard]] auto
  getNearestSiteDownLeft(const Point& p, bool proper = false,
                         bool sameZone = false) const -> Index;
  [[nodiscard]] auto
  getNearestSiteDownRight(const Point& p, bool proper = false,
                          bool sameZone = false) const -> Index;
  [[nodiscard]] auto getSiteAt(const Point& p) const -> Index;
  [[nodiscard]] auto getSitesInZone(const Zone& z) const -> std::vector<Index>;
  [[nodiscard]] auto
  withConfig(const Configuration& config) const -> Architecture;
  [[nodiscard]] auto getPositionOffsetBy(const Point& p, const Number& rows,
                                         const Number& cols) const -> Point;

private:
  [[nodiscard]] auto getRowsInZone(const Zone& z) const -> std::vector<Number>;
  [[nodiscard]] auto getColsInZone(const Zone& z) const -> std::vector<Number>;
  [[nodiscard]] auto getRows() const -> std::vector<Number>;
  [[nodiscard]] auto getCols() const -> std::vector<Number>;
};
} // namespace na
