//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#include "Architecture.hpp"

#include "nlohmann/json.hpp"
#include "operations/OpType.hpp"

#include <algorithm>
#include <climits>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace na {

Architecture::Architecture(const std::string& jsonFn,
                           const std::string& csvFn) {
  std::ifstream jsonS(jsonFn);
  if (!jsonS.good()) {
    throw std::runtime_error("Could not open JSON file.");
  }
  std::ifstream csvS(csvFn);
  if (!csvS.good()) {
    throw std::runtime_error("Could not open CSV file.");
  }
  *this = Architecture(jsonS, csvS);
}

Architecture::Architecture(std::istream& jsonS, std::istream& csvS) {
  nlohmann::json data;
  // load CSV
  // auxiliary variables
  std::string line;
  try {
    csvS >> line;                         // skip first line, i.e. header
    while (csvS >> line) {                // read one line
      std::stringstream lineStream(line); // make a stream of this line
      std::string       sX;
      std::string       sY;

      std::getline(lineStream, sX, ',');
      std::getline(lineStream, sY);

      sites.emplace_back(
          std::stoi(sX),
          std::stoi(sY)); // make point from x and y (convert to int)
    }
  } catch (std::exception& e) {
    throw std::runtime_error(
        "While parsing the CSV file, the following error occurred: " +
        std::string(e.what()));
  }
  // load JSON
  try {
    jsonS >> data;
  } catch (std::exception& e) {
    throw std::runtime_error(
        "While parsing the JSON file, the following error occurred: " +
        std::string(e.what()));
  }
  try {
    // load rest of JSON
    name = data["name"];
    std::map<std::string, Zone> nameToZone;
    for (const auto& zone : data["zones"]) {
      nameToZone[zone["name"]] = zones.size();
      ZoneProperties zp{zone["name"], zone["xmin"], zone["xmax"],
                        zone["ymin"], zone["ymax"], zone["fidelity"]};
      zones.emplace_back(zp);
    }
    for (auto const& zone : data["initialZones"]) {
      initialZones.emplace_back(nameToZone.find(zone)->second);
    }
    for (auto const& op : data["operations"]) {
      const std::string        opName = op["name"];
      const OpType             ty     = {qc::opTypeFromString(opName),
                                         opName.find_first_not_of('c')};
      const Scope              sc     = getScopeOfString(op["type"]);
      std::unordered_set<Zone> zo     = {};
      for (auto const& zs : op["zones"]) {
        zo.emplace(nameToZone.find(zs)->second);
      }
      const Value               fi = op["fidelity"];
      const Value               ti = op["time"];
      OperationProperties const o  = {sc, zo, ti, fi};
      gateSet.emplace(ty, o);
    }
    decoherenceTimes = Architecture::DecoherenceTimes(
        data["decoherence"]["t1"], data["decoherence"]["t2"]);
    for (const auto& sh : data["shuttling"]) {
      ShuttlingProperties sp{sh["rows"],          sh["columns"],
                             sh["xmin"],          sh["xmax"],
                             sh["ymin"],          sh["ymax"],
                             sh["move"]["speed"], sh["move"]["fidelity"],
                             sh["load"]["time"],  sh["load"]["fidelity"],
                             sh["store"]["time"], sh["store"]["fidelity"]};
      shuttling.emplace_back(sp);
    }
    minAtomDistance   = data["minAtomDistance"];
    interactionRadius = data["interactionRadius"];
  } catch (std::exception& e) {
    throw std::runtime_error(
        "While reading the JSON data, the following error occurred: " +
        std::string(e.what()));
  }
}

[[nodiscard]] auto Architecture::getZoneAt(const Point& p) const -> Zone {
  const auto& it =
      std::find_if(zones.cbegin(), zones.cend(), [&](const auto& zProp) {
        return p.x >= zProp.minX && p.x <= zProp.maxX && p.y >= zProp.minY &&
               p.y <= zProp.maxY;
      });
  if (it == zones.cend()) {
    std::stringstream ss;
    ss << "The point " << p << " is not in any zone.";
    throw std::invalid_argument(ss.str());
  }
  return static_cast<Zone>(std::distance(zones.cbegin(), it));
}

auto Architecture::isAllowedLocally(const OpType& t) const -> bool {
  const auto it = gateSet.find(t);
  return it != gateSet.end() && it->second.scope == Scope::Local;
}
auto Architecture::isAllowedLocally(const OpType& t, const Zone& zone) const
    -> bool {
  if (!isAllowedLocally(t)) {
    return false; // gate not supported at all
  }
  const auto  it        = gateSet.find(t);
  const auto& gateZones = it->second.zones;
  // zone exists in gateZones
  return gateZones.find(zone) != gateZones.end();
}

auto Architecture::isAllowedLocallyAtSite(const OpType& t,
                                          const Index&  qubit) const -> bool {
  const auto zone = getZoneOfSite(qubit);
  return isAllowedLocally(t, zone);
}

auto Architecture::isAllowedLocallyAt(const OpType& t, const Point& p) const
    -> bool {
  const auto& it =
      std::find_if(zones.cbegin(), zones.cend(), [&](const auto& zProp) {
        return p.x >= zProp.minX && p.x <= zProp.maxX && p.y >= zProp.minY &&
               p.y <= zProp.maxY;
      });
  if (it == zones.cend()) {
    std::stringstream ss;
    ss << "The point " << p << " is not in any zone.";
    throw std::invalid_argument(ss.str());
  }
  return isAllowedLocally(t,
                          static_cast<Zone>(std::distance(zones.cbegin(), it)));
}

auto Architecture::isAllowedGlobally(const OpType& t) const -> bool {
  const auto it = gateSet.find(t);
  return it != gateSet.end() && it->second.scope == Scope::Global;
}

auto Architecture::isAllowedGlobally(const OpType& t, const Zone& zone) const
    -> bool {
  if (!isAllowedGlobally(t)) {
    return false; // gate not supported at all
  }
  const auto  it        = gateSet.find(t);
  const auto& gateZones = it->second.zones;
  // zone exists in gateZones
  return gateZones.find(zone) != gateZones.end();
}
[[nodiscard]] auto Architecture::getRows() const -> std::vector<Number> {
  std::unordered_set<Number> rows;
  std::for_each(sites.cbegin(), sites.cend(),
                [&](const auto& s) { rows.insert(s.y); });
  return {rows.cbegin(), rows.cend()};
}
[[nodiscard]] auto Architecture::getCols() const -> std::vector<Number> {
  std::unordered_set<Number> cols;
  std::for_each(sites.cbegin(), sites.cend(),
                [&](const auto& s) { cols.insert(s.x); });
  return {cols.cbegin(), cols.cend()};
}
[[nodiscard]] auto Architecture::getRowsInZone(const Zone& z) const
    -> std::vector<Number> {
  std::unordered_set<Number> rows;
  std::for_each(sites.cbegin(), sites.cend(), [&](const auto& s) {
    if (s.y >= zones[z].minY && s.y <= zones[z].maxY) {
      rows.insert(s.y);
    }
  });
  return {rows.cbegin(), rows.cend()};
}
[[nodiscard]] auto Architecture::getColsInZone(const Zone& z) const
    -> std::vector<Number> {
  std::unordered_set<Number> cols;
  std::for_each(sites.cbegin(), sites.cend(), [&](const auto& s) {
    if (s.x >= zones[z].minX && s.x <= zones[z].maxX) {
      cols.insert(s.x);
    }
  });
  return {cols.cbegin(), cols.cend()};
}
[[nodiscard]] auto Architecture::isInSameRow(const Index &i, const Index& j) const -> bool {
  return getPositionOfSite(i).y == getPositionOfSite(j).y;
}
[[nodiscard]] auto Architecture::isInSameCol(const Index &i, const Index& j) const -> bool {
  return getPositionOfSite(i).x == getPositionOfSite(j).x;
}
[[nodiscard]] auto Architecture::getNrowsInZone(const Zone& z) const -> Index {
  return Architecture::getRowsInZone(z).size();
}
[[nodiscard]] auto Architecture::getNcolsInZone(const Zone& z) const -> Index {
  return Architecture::getColsInZone(z).size();
}
[[nodiscard]] auto Architecture::getSitesInRow(const Zone&  z,
                                               const Index& row) const
    -> std::vector<Index> {
  const auto         y = Architecture::getRowsInZone(z)[row];
  std::vector<Index> atoms;
  std::vector<Index> enumerate(sites.size());
  std::iota(enumerate.begin(), enumerate.end(), 0);
  std::copy_if(enumerate.cbegin(), enumerate.cend(), std::back_inserter(atoms),
               [&](const auto& i) {
                 const auto& s = sites[i];
                 return s.y == y && s.x >= zones[z].minX &&
                        s.x <= zones[z].maxX;
               });
  return atoms;
}
[[nodiscard]] auto Architecture::getSitesInCol(const Zone&  z,
                                               const Index& col) const
    -> std::vector<Index> {
  const auto         x = Architecture::getColsInZone(z)[col];
  std::vector<Index> atoms;
  std::vector<Index> enumerate(sites.size());
  std::iota(enumerate.begin(), enumerate.end(), 0);
  std::copy_if(enumerate.cbegin(), enumerate.cend(), std::back_inserter(atoms),
               [&](const auto& i) {
                 const auto& s = sites[i];
                 return s.x == x && s.y >= zones[z].minY &&
                        s.y <= zones[z].maxY;
               });
  return atoms;
}
[[nodiscard]] auto Architecture::getSitesInZone(const Zone& z) const
    -> std::vector<Index> {
  std::vector<Index> atoms;
  std::vector<Index> enumerate(sites.size());
  std::iota(enumerate.begin(), enumerate.end(), 0);
  std::copy_if(enumerate.cbegin(), enumerate.cend(), std::back_inserter(atoms),
               [&](const auto& i) {
                 const auto& s = sites[i];
                 return s.x >= zones[z].minX && s.x <= zones[z].maxX &&
                        s.y >= zones[z].minY && s.y <= zones[z].maxY;
               });
  return atoms;
}
[[nodiscard]] auto Architecture::getRowInZoneOf(const Index& i) const -> Index {
  const auto& rows = Architecture::getRowsInZone(getZoneOfSite(i));
  const auto& it   = std::find(rows.cbegin(), rows.cend(), sites[i].y);
  assert(it != rows.cend());
  return static_cast<Index>(std::distance(rows.cbegin(), it));
}
[[nodiscard]] auto Architecture::getColInZoneOf(const Index& i) const -> Index {
  const auto& cols = Architecture::getColsInZone(getZoneOfSite(i));
  const auto& it   = std::find(cols.cbegin(), cols.cend(), sites[i].x);
  assert(it != rows.cend());
  return static_cast<Index>(std::distance(cols.cbegin(), it));
}
[[nodiscard]] auto Architecture::getNearestXLeft(const Number& x) const
    -> Number {
  const auto& cols = getCols();
  if (!std::any_of(cols.cbegin(), cols.cend(),
                   [&](const auto& c) { return c <= x; })) {
    return x;
  }
  return std::accumulate(cols.cbegin(), cols.cend(), LONG_LONG_MIN,
                         [&](const auto& acc, const auto& c) {
                           return acc < c and c < x ? c : acc;
                         });
}
[[nodiscard]] auto Architecture::getNearestXRight(const Number& x) const
    -> Number {
  const auto& cols = getCols();
  if (!std::any_of(cols.cbegin(), cols.cend(),
                   [&](const auto& c) { return c >= x; })) {
    return x;
  }
  return std::accumulate(cols.cbegin(), cols.cend(), LONG_LONG_MAX,
                         [&](const auto& acc, const auto& c) {
                           return acc > c and c > x ? c : acc;
                         });
}
[[nodiscard]] auto Architecture::getNearestYUp(const Number& y) const
    -> Number {
  const auto& rows = getRows();
  if (!std::any_of(rows.cbegin(), rows.cend(),
                   [&](const auto& c) { return c <= y; })) {
    return y;
  }
  return std::accumulate(rows.cbegin(), rows.cend(), LONG_LONG_MIN,
                         [&](const auto& acc, const auto& r) {
                           return acc < r and r < y ? r : acc;
                         });
}
[[nodiscard]] auto Architecture::getNearestYDown(const Number& y) const
    -> Number {
  const auto& rows = getRows();
  if (!std::any_of(rows.cbegin(), rows.cend(),
                   [&](const auto& c) { return c >= y; })) {
    return y;
  }
  return std::accumulate(rows.cbegin(), rows.cend(), LONG_LONG_MAX,
                         [&](const auto& acc, const auto& r) {
                           return acc > r and r > y ? r : acc;
                         });
}
[[nodiscard]] auto Architecture::getNearestSiteLeft(const Point& p) const
    -> Index {
  const auto& it =
      std::find_if(sites.cbegin(), sites.cend(),
                   [&](const auto& s) { return s.y == p.y and s.x <= p.x; });
  if (it == sites.cend()) {
    throw std::invalid_argument("No site found.");
  }
  return static_cast<std::size_t>(std::distance(sites.cbegin(), it));
}
[[nodiscard]] auto Architecture::getNearestSiteRight(const Point& p) const
    -> Index {
  const auto& it =
      std::find_if(sites.cbegin(), sites.cend(),
                   [&](const auto& s) { return s.y == p.y and s.x >= p.x; });
  if (it == sites.cend()) {
    throw std::invalid_argument("No site found.");
  }
  return static_cast<std::size_t>(std::distance(sites.cbegin(), it));
}
[[nodiscard]] auto Architecture::getNearestSiteUp(const Point& p) const
    -> Index {
  const auto& it =
      std::find_if(sites.cbegin(), sites.cend(),
                   [&](const auto& s) { return s.x == p.x and s.y <= p.y; });
  if (it == sites.cend()) {
    throw std::invalid_argument("No site found.");
  }
  return static_cast<std::size_t>(std::distance(sites.cbegin(), it));
}
[[nodiscard]] auto Architecture::getNearestSiteDown(const Point& p) const
    -> Index {
  const auto& it =
      std::find_if(sites.cbegin(), sites.cend(),
                   [&](const auto& s) { return s.x == p.x and s.y >= p.y; });
  if (it == sites.cend()) {
    throw std::invalid_argument("No site found.");
  }
  return static_cast<std::size_t>(std::distance(sites.cbegin(), it));
}
[[nodiscard]] auto Architecture::getSiteAt(const Point& p) const -> Index {
  std::vector<std::size_t> enumerate(sites.size());
  std::iota(enumerate.begin(), enumerate.end(), 0ULL);
  const auto& it = std::find_if(enumerate.cbegin(), enumerate.cend(),
                                [&](const auto i) { return sites[i] == p; });
  if (it == enumerate.end()) {
    throw std::invalid_argument("No site at this position.");
  }
  return *it;
}

} // namespace na
