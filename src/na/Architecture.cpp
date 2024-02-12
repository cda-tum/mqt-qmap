//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"

#include "nlohmann/json.hpp"
#include "operations/OpType.hpp"

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
      qc::OpType const         ty = qc::opTypeFromString(op["name"]);
      Scope const              sc = getScopeOfString(op["type"]);
      std::unordered_set<Zone> zo = {};
      for (auto const& zs : op["zones"]) {
        zo.emplace(nameToZone.find(zs)->second);
      }
      Value const               fi = op["fidelity"];
      Value const               ti = op["time"];
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

auto Architecture::isAllowedLocallyIn(const qc::OpType gate,
                                      const Zone       zone) const -> bool {
  const auto it = gateSet.find(gate);
  if (it == gateSet.end()) {
    return false; // gate not supported at all
  }
  if (it->second.scope != Scope::Local) {
    return false; // gate cannot be applied individually
  }
  const auto& gateZones = it->second.zones;
  // zone exists in gateZones
  return gateZones.find(zone) != gateZones.end();
}

auto Architecture::isAllowedLocallyAt(const qc::OpType gate,
                                      const Index      qubit) const -> bool {
  const auto zone = getZone(qubit);
  return isAllowedLocallyIn(gate, zone);
}

auto Architecture::isAllowedGlobally(const qc::OpType gate,
                                     const Zone       zone) const -> bool {
  const auto it = gateSet.find(gate);
  if (it == gateSet.end()) {
    return false; // gate not supported at all
  }
  if (it->second.scope != Scope::Global) {
    return false; // gate cannot be applied globally
  }
  const auto& gateZones = it->second.zones;
  // zone exists in gateZones
  return gateZones.find(zone) != gateZones.end();
}
} // namespace na
