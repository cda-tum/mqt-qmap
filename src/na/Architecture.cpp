//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "na/Architecture.hpp"

#include "nlohmann/json.hpp"
#include "nlohmann/json_fwd.hpp"
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
std::string getDirnameOf(std::string& filename) {
  return std::filesystem::path(filename).parent_path().u8string();
}

Architecture::Architecture(std::string& filename) {
  std::ifstream  fs(filename);
  nlohmann::json data;
  if (!fs.is_open()) {
    throw std::runtime_error("Could not open file.");
  }
  try {
    fs >> data;
    name = data["name"];
    // load CSV
    std::stringstream ss;
    ss << getDirnameOf(filename) << "/" << data["grid"].get<std::string>();
    std::ifstream gridFs(ss.str());
    if (!gridFs.is_open()) {
      std::stringstream msg;
      msg << "Could not open file with file name: " << ss.str() << ".";
      throw std::runtime_error(msg.str());
    }
    // auxiliary variables
    std::string                 line;
    std::map<std::string, Zone> nameToZone;

    gridFs >> line;                       // skip first line, i.e. header
    while (gridFs >> line) {              // read one line
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
      std::getline(lineStream, sType, ',');

      if (nameToZone.try_emplace(sZone, nameToZone.size()).second) {
        zones.emplace_back(sZone);
      }
      sites.emplace_back(
          Point(static_cast<std::uint32_t>(std::stoul(sX)),
                static_cast<std::uint32_t>(std::stoul(
                    sY))), // make point from x and y (convert to int)
          nameToZone.find(sZone)->second, // find the id of the zone
          getTypeOfString(sType));        // convert to type
    }
    nSites = sites.size();
    // load rest of JSON
    for (auto const& op : data["operations"]) {
      qc::OpType const ty = qc::opTypeFromString(op["name"]);
      Scope const      sc = getScopeOfString(op["type"]);
      std::set<Zone>   zo = {};
      for (auto const& zs : op["zones"]) {
        zo.emplace(nameToZone.find(zs)->second);
      }
      Value const     fi = op["fidelity"];
      Value const     ti = op["time"];
      Operation const o  = {ty, sc, zo, ti, fi};
      operations.emplace(ty, o);
    }
    decoherenceTimes = Architecture::DecoherenceTimes(
        data["decoherence"]["t1"], data["decoherence"]["t2"]);
    nAods      = data["AOD"]["number"];
    shutteling = Architecture::Shutteling(
        data["AOD"]["rows"], data["AOD"]["columns"],
        data["AOD"]["move"]["speed"], data["AOD"]["move"]["fidelity"],
        data["AOD"]["activate"]["time"], data["AOD"]["activate"]["fidelity"],
        data["AOD"]["deactivate"]["time"],
        data["AOD"]["deactivate"]["fidelity"]);
    minAtomDistance   = data["minAtomDistance"];
    interactionRadius = data["interactionRadius"];
  } catch (std::exception& e) {
    throw std::runtime_error(
        "While parsing the JSON file, the following error occurred: " +
        std::string(e.what()));
  }
}

bool Architecture::isAllowedLocally(qc::OpType gate, Index qubit) {
  auto it = operations.find(gate);
  if (it == operations.end()) {
    return false; // gate not supported at all
  }
  if (it->second.scope != Scope::Local) {
    return false; // gate cannot be applied individually
  }
  auto zone = getZone(qubit);
  auto zit  = it->second.zones.find(zone);
  return zit == it->second.zones.end();
}

bool Architecture::isAllowedGlobally(qc::OpType gate, Zone zone) {
  auto it = operations.find(gate);
  if (it == operations.end()) {
    return false; // gate not supported at all
  }
  auto zit = it->second.zones.find(zone);
  return zit == it->second.zones.end();
}
} // namespace na
