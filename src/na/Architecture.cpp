//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"

#include "nlohmann/json.hpp"
#include "operations/OpType.hpp"

#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

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
  std::string                 line;
  std::map<std::string, Zone> nameToZone;
  try {
    csvS >> line;                         // skip first line, i.e. header
    while (csvS >> line) {                // read one line
      std::stringstream lineStream(line); // make a stream of this line
      std::string       sAddr;
      std::string       sX;
      std::string       sY;
      std::string       sZone;
      std::string       sType;

      std::getline(lineStream, sAddr, ','); // read until the next separator
      std::getline(lineStream, sX, ',');
      std::getline(lineStream, sY, ',');
      std::getline(lineStream, sZone, ',');
      std::getline(lineStream, sType);

      const auto& [it, success] =
          nameToZone.try_emplace(sZone, nameToZone.size());
      if (success) {
        zones.emplace_back(sZone);
      }
      sites.emplace_back(
          Point(std::stoi(sX),
                std::stoi(sY)),    // make point from x and y (convert to int)
          it->second,              // find the id of the zone
          getTypeOfString(sType)); // convert to type
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
    name            = data["name"];
    nShuttlingUnits = sites.size();
    // load rest of JSON
    for (auto const& op : data["operations"]) {
      const std::string                    opName = op["name"];
      const std::pair<qc::OpType, Number> ty     = {
          qc::opTypeFromString(opName), opName.find_first_not_of('c')};
      const Scope              sc = getScopeOfString(op["type"]);
      std::unordered_set<Zone> zo = {};
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
    const auto& aodData = data["AOD"];
    nShuttlingUnits     = aodData["number"];
    shuttling           = {aodData["rows"],
                           aodData["columns"],
                           aodData["move"]["speed"],
                           aodData["move"]["fidelity"],
                           aodData["activate"]["time"],
                           aodData["activate"]["fidelity"],
                           aodData["deactivate"]["time"],
                           aodData["deactivate"]["fidelity"]};
    minAtomDistance     = data["minAtomDistance"];
    interactionRadius   = data["interactionRadius"];
  } catch (std::exception& e) {
    throw std::runtime_error(
        "While reading the JSON data, the following error occurred: " +
        std::string(e.what()));
  }
}

auto Architecture::isAllowedLocally(const qc::OpType gate, Number nctrl) const
    -> bool {
  const auto it = gateSet.find({gate, nctrl});
  return it != gateSet.end() && it->second.scope == Scope::Local;
}
auto Architecture::isAllowedLocally(const qc::OpType gate, const Zone zone,
                                    Number nctrl) const -> bool {
  if (!isAllowedLocally(gate, nctrl)) {
    return false; // gate not supported at all
  }
  const auto  it        = gateSet.find({gate, nctrl});
  const auto& gateZones = it->second.zones;
  // zone exists in gateZones
  return gateZones.find(zone) != gateZones.end();
}

auto Architecture::isAllowedLocallyAt(const qc::OpType gate, const Index qubit,
                                      Number nctrl) const -> bool {
  const auto zone = getZone(qubit);
  return isAllowedLocally(gate, zone, nctrl);
}

auto Architecture::isAllowedGlobally(const qc::OpType gate, Number nctrl) const
    -> bool {
  const auto it = gateSet.find({gate, nctrl});
  return it != gateSet.end() && it->second.scope == Scope::Global;
}

auto Architecture::isAllowedGlobally(const qc::OpType gate, const Zone zone,
                                     Number nctrl) const -> bool {
  if (!isAllowedGlobally(gate, nctrl)) {
    return false; // gate not supported at all
  }
  const auto  it        = gateSet.find({gate, nctrl});
  const auto& gateZones = it->second.zones;
  // zone exists in gateZones
  return gateZones.find(zone) != gateZones.end();
}
} // namespace na
