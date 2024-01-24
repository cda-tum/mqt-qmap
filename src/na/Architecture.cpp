//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "na/Architecture.hpp"

#include "nlohmann/json.hpp"

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace na {
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
    std::string                 gridFname = data["grid"];
    std::ifstream               gridFs(gridFname);
    std::string                 line;
    std::map<std::string, Zone> nameToZone();

    while (std::getline(gridFs, line)) {  // read one line
      std::stringstream lineStream(line); // make a stream of this line
      std::string       sAddr, sX, sY, sZone, sType;

      std::getline(lineStream, sAddr, ','); // read until the next separator
      std::getline(lineStream, sX, ',');
      std::getline(lineStream, sY, ',');
      std::getline(lineStream, sZone, ',');
      std::getline(lineStream, sType, ',');

      if( nameToZone.try_emplace(sZone, nameToZone.count())) { // FIXME
        zones.emplace_back(sZone);
      }
      sites.emplace_back(Point(stoi(sX), stoi(sY)), nameToZone.find(sZone), getTypeOfString(sType));
    }
    nSites = sites.size();
    // load rest of JSON
    decoherenceTimes = Architecture::DecoherenceTimes(
        data["decoherence"]["t1"], data["decoherence"]["t2"]); // FIXME
    nAods             = data["AOD"]["number"];
    minAtomDistance   = data["AOD"]["minAtomDistance"];
    interactionRadius = data["interactionRadius"];
  } catch (std::exception& e) {
    throw std::runtime_error("Could not parse JSON file." +
                             std::string(e.what()));
  }
}
} // namespace na
