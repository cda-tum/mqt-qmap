//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "na/Architecture.hpp"

#include "OpType.hpp"
#include "nlohmann/json.hpp"
#include "nlohmann/json_fwd.hpp"

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
    ss << getDirnameOf(filename) << data["grid"];
    std::ifstream gridFs(ss.str());
    // auxiliary variables
    std::string                 line;
    std::map<std::string, Zone> nameToZone;

    while (std::getline(gridFs, line)) {  // read one line
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
          Point(stoi(sX), stoi(sY)), // make point from x and y (convert to int)
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
      Value const fi   = op["fidelity"];
      Value const ti   = op["time"];
      Operation const o = {ty, sc, zo, ti, fi};
      operations.emplace(ty, o);
    }
    decoherenceTimes = Architecture::DecoherenceTimes(
        data["decoherence"]["t1"], data["decoherence"]["t2"]);
    nAods             = data["AOD"]["number"];
    shutteling.rows          = data["AOD"]["rows"];
    shutteling.cols          = data["AOD"]["cols"];
    shutteling.speed = data["AOD"]["move"]["speed"];
    shutteling.fidelity = data["AOD"]["move"]["fidelity"];
    shutteling.activationTime = data["AOD"]["activate"]["time"];
    shutteling.activationFidelity = data["AOD"]["activate"]["fidelity"];
    shutteling.deactivationTime = data["AOD"]["deactivate"]["time"];
    shutteling.deactivationFidelity = data["AOD"]["deactivate"]["fidelity"];
    minAtomDistance   = data["minAtomDistance"];
    interactionRadius = data["interactionRadius"];
  } catch (std::exception& e) {
    throw std::runtime_error("Could not parse JSON file." +
                             std::string(e.what()));
  }
}
} // namespace na
