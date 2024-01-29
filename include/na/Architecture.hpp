//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "operations/OpType.hpp"
#include "na/Definitions.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace na {

/// The scope of an operation (Global or Local)
enum class Scope: uint8_t { Global, Local };
static constexpr std::map<std::string, Scope> STRING_TO_SCOPE = {
    {"Global", Scope::Global}, {"Local", Scope::Local}};
/**
 * @brief Get the Scope of a gate from a string
 *
 * @param s the name
 * @return Type
 */
inline Scope getScopeOfString(const std::string& s) {
  auto it = STRING_TO_SCOPE.find(s);
  return it->second;
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
  auto it = STRING_TO_TYPE.find(s);
  return it->second;
}

/**
 * @brief For Indices
 */
using Index = std::size_t;
/**
 * @brief The zones are just stored as int
 */
using Zone = Index;
/**
 * @brief A site is defined by a position, a zone, and a type
 */
using Site = std::tuple<Point, Zone, Type>;
/**
 * @brief Any double-valued property
 */
using Value = qc::fp;
/**
 * @brief Any information on numbers of something
 */
using Number = std::size_t

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
  struct DecoherenceTimes {
    Value t1                              = 0;
    Value t2                              = 0;
    Value tEff                            = 0;
    DecoherenceTimes()                    = default;
    DecoherenceTimes(DecoherenceTimes& t) = default;
    inline DecoherenceTimes(Value t1, Value t2)
        : t1(t1), t2(t2), tEff(t1 * t2 / (t1 + t2)) {}
    virtual ~DecoherenceTimes() = default;
    explicit          operator double() const { return tEff; }
    DecoherenceTimes& operator=(DecoherenceTimes&& t) noexcept {
      if (this != &t) {
        t1   = t.t1;
        t2   = t.t2;
        tEff = t.tEff;
      }
      return *this;
    }
  };
  class Operation {
  public:
    qc::OpType type;  // the type of the gate, use also RY for global ones here
    Scope      scope; // local or global
    std::set<Zone> zones;    // the zones where the gate can be applied
    Value          time;     // the time the gate takes to be applied
    Value          fidelity; // the fidelity of the gate
    Operation(const Operation& op) = default;
    Operation(qc::OpType type, Scope scope, std::set<Zone> zones, Value time,
              Value fidelity)
        : type(type), scope(scope), zones(std::move(zones)), time(time),
          fidelity(fidelity) {}
  };
  class Shuttling {
  public:
    Number rows;
    Number cols;
    Value  speed;
    Value  fidelity;
    Value  activationTime;
    Value  activationFidelity;
    Value  deactivationTime;
    Value  deactivationFidelity;
    Shuttling()               = default;
    Shuttling(Shuttling& sh) = default;
    Shuttling(Number rs, Number cs, Value sp, Value fi, Value ta, Value fa,
               Value td, Value fd)
        : rows(rs), cols(cs), speed(sp), fidelity(fi), activationTime(ta),
          activationFidelity(fa), deactivationTime(td),
          deactivationFidelity(fd) {}
    Shuttling& operator=(Shuttling&& s) noexcept {
      if (this != &s) {
        rows                 = s.rows;
        cols                 = s.cols;
        speed                = s.speed;
        fidelity             = s.fidelity;
        activationTime       = s.activationTime;
        activationFidelity   = s.activationFidelity;
        deactivationTime     = s.deactivationTime;
        deactivationFidelity = s.deactivationFidelity;
      }
      return *this;
    }
  };

protected:
  std::string name; // the name of the architecture
  std::vector<std::string>
      zones; // a mapping from zones (int) to their name from the config
  std::vector<Site> sites;  // a vector of sites (Position, Zone, Type)
  Number            nSites; // number of sites
  std::map<qc::OpType, Operation>
      operations; // all possible operations by their type, i.e. gate set
  DecoherenceTimes decoherenceTimes; // the decoherence characteristic
  Number           nAods;            // number of AODs for atom movement
  Shuttling       shuttling;       // all properties regarding AODs
  Value minAtomDistance;   // minimal distance that must be kept between atoms
  Value interactionRadius; // the Rydberg radius

public:
  /**
   * @brief Import a new architecture from a file.
   *
   * @param filename The path to the JSON file
   */
  explicit Architecture(std::string& filename);
  virtual ~Architecture() = default;
  [[nodiscard]] inline auto getName() { return name; }
  [[nodiscard]] inline auto getNZones() { return zones.size(); }
  [[nodiscard]] inline auto getZoneLabel(const Index& i) { return zones.at(i); }
  [[nodiscard]] inline auto getNSites() const { return nSites; }
  [[nodiscard]] auto        getType(const Index& i) {
    return std::get<Type>(sites.at(i));
  }
  [[nodiscard]] auto getZone(const Index& i) {
    return std::get<Zone>(sites.at(i));
  }
  [[nodiscard]] inline auto getPos(const Index& i) {
    return std::get<Point>(sites.at(i));
  }
  [[nodiscard]] inline auto getDecoherenceTimes() { return decoherenceTimes; }
  [[nodiscard]] inline auto getNAods() const { return nAods; }
  [[nodiscard]] inline auto getShuttling() { return shuttling; }
  [[nodiscard]] inline auto getMinAtomDistance() const {
    return minAtomDistance;
  }
  [[nodiscard]] inline auto getInteractionRadius() const {
    return interactionRadius;
  }
  [[nodiscard]] auto getOperationByOpType(const qc::OpType& t) {
    auto it = operations.find(t);
    if (it == operations.end()) {
      throw std::invalid_argument(
          "This operation is not supported by this architecture.");
    }
    return it->second;
  }
  /**
   * @brief Returns the distance between two sites.
   *
   * @param i address of first site
   * @param j address of second site
   * @return the distance in µm
   */
  [[nodiscard]] auto getDistance(Index i, Index j) {
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
  [[nodiscard]] bool isAllowedLocally(qc::OpType gate, Index qubit);
  /**
   * @brief Checks whether the gate is a global gate for this Zone.
   *
   * @param gate the gate
   * @param zone the zone
   * @return true if the gate is global and applicable in this zone,
   * @return false otherwise
   */
  [[nodiscard]] bool isAllowedGlobally(qc::OpType gate, Zone zone);
};
} // namespace na
