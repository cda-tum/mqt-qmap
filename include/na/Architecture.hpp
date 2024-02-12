//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include "Definitions.hpp"
#include "na/Definitions.hpp"
#include "operations/OpType.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {

/// The scope of an operation (Global or Local)
enum class Scope : uint8_t { Global, Local };
static const std::unordered_map<std::string, Scope> STRING_TO_SCOPE = {
    {"Global", Scope::Global}, {"Local", Scope::Local},
    {"global", Scope::Global}, {"local", Scope::Local}};
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

/**
 * @brief The type of a site (SLM or AOD)
 * @details SLM comprises both AOD and SLM; AOD denotes AOD only.
 */
enum class Type { SLM, AOD };
static std::map<std::string, Type> const STRING_TO_TYPE = {{"SLM", Type::SLM},
                                                           {"AOD", Type::AOD}};
/**
 * @brief Get the Type of a site (SLM or AOD) from a string
 *
 * @param s the name
 * @return Type
 */
inline Type getTypeOfString(const std::string& s) {
  if (auto it = STRING_TO_TYPE.find(s); it != STRING_TO_TYPE.end()) {
    return it->second;
  }
  std::stringstream ss;
  ss << "The type " << s << " is not supported.";
  throw std::invalid_argument(ss.str());
}

/// For Indices
using Index = std::size_t;
/// The zones are just stored as int
using Zone = Index;
/// A site is defined by a position, a zone, and a type
using Site = std::tuple<Point, Zone, Type>;
/// Any double-valued property
using Value = qc::fp;
/// Any information on numbers of something
using Number = std::size_t;

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
    inline DecoherenceTimes(Value t1, Value t2)
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
    Number rows                 = 0; // maximum number of rows in one AOD
    Number cols                 = 0; // maximum number of columns in one AOD
    Value  speed                = 0; // speed of the AOD in µm/µs
    Value  fidelity             = 1; // fidelity during the shuttling
    Value  activationTime       = 0; // time to activate the AOD in µs
    Value  activationFidelity   = 1; // fidelity of the activation
    Value  deactivationTime     = 0; // time to deactivate the AOD in µs
    Value  deactivationFidelity = 1; // fidelity of the deactivation
  };

protected:
  std::string name{}; // the name of the architecture
  std::vector<std::string>
      zones{}; // a mapping from zones (int) to their name from the config
  std::vector<Site> sites{}; // a vector of sites (Position, Zone, Type)
  std::map<qc::OpType, OperationProperties>
      gateSet{}; // all possible operations by their type, i.e. gate set
  DecoherenceTimes    decoherenceTimes{};  // the decoherence characteristic
  Number              nShuttlingUnits = 0; // number of AODs for atom movement
  ShuttlingProperties shuttling{};         // all properties regarding AODs
  Value minAtomDistance = 0; // minimal distance that must be kept between atoms
  Value interactionRadius = 0; // the Rydberg radius

public:
  /**
   * @brief Import a new architecture from a file.
   *
   * @param filename The path to the JSON file
   */
  Architecture(const std::string& jsonFn, const std::string& csvFn);
  Architecture(std::istream& jsonS, std::istream& csvS);
  virtual ~Architecture() = default;

  [[nodiscard]] auto getName() const { return name; }
  [[nodiscard]] auto getNZones() const { return zones.size(); }
  [[nodiscard]] auto getZoneLabel(const Index& i) const { return zones.at(i); }
  [[nodiscard]] auto getNSites() const { return sites.size(); }
  [[nodiscard]] auto getType(const Index& i) const {
    return std::get<Type>(sites.at(i));
  }
  [[nodiscard]] auto getZone(const Index& i) const {
    return std::get<Zone>(sites.at(i));
  }
  [[nodiscard]] auto getPos(const Index& i) const {
    return std::get<Point>(sites.at(i));
  }
  [[nodiscard]] auto getDecoherenceTimes() const -> DecoherenceTimes {
    return decoherenceTimes;
  }
  [[nodiscard]] auto getNShuttlingUnits() const { return nShuttlingUnits; }
  [[nodiscard]] auto getShuttling() const -> ShuttlingProperties {
    return shuttling;
  }
  [[nodiscard]] auto getMinAtomDistance() const { return minAtomDistance; }
  [[nodiscard]] auto getInteractionRadius() const { return interactionRadius; }
  [[nodiscard]] auto getOpPropsByOpType(const qc::OpType& t) const
      -> OperationProperties {
    if (auto it = gateSet.find(t); it != gateSet.end()) {
      return it->second;
    }
    std::stringstream ss;
    ss << "The operation " << qc::toString(t) << " is not supported.";
    throw std::invalid_argument(ss.str());
  }
  /**
   * @brief Returns the distance between two sites.
   *
   * @param i address of first site
   * @param j address of second site
   * @return the distance in µm
   */
  [[nodiscard]] auto getDistance(Index i, Index j) const {
    return (getPos(j) - getPos(i)).length();
  }
  /**
   * @brief Checks whether the gate can be applied (locally) on this qubit.
   *
   * @param gate the gate
   * @param qubit the qubit
   * @return true if the gate is a local operation and available in the zone of
   * the qubit,
   * @return false otherwise
   */
  [[nodiscard]] auto isAllowedLocallyAt(qc::OpType gate, Index qubit) const -> bool;
  /**
   * @brief Checks whether the gate can be applied (locally) in this zone.
   *
   * @param gate the gate
   * @param qubit the qubit
   * @return true if the gate is a local operation and available in the zone of
   * the qubit,
   * @return false otherwise
   */
  [[nodiscard]] auto isAllowedLocallyIn(qc::OpType gate, Zone zone) const -> bool;
  /**
   * @brief Checks whether the gate is a global gate for this Zone.
   *
   * @param gate the gate
   * @param zone the zone
   * @return true if the gate is global and applicable in this zone,
   * @return false otherwise
   */
  [[nodiscard]] auto isAllowedGlobally(qc::OpType gate, Zone zone) const -> bool;
};
} // namespace na
