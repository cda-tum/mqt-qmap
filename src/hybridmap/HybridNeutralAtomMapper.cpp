//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/HybridNeutralAtomMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/MoveToAodConverter.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomLayer.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "operations/OpType.hpp"
#include "operations/Operation.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <iomanip>

namespace na {

qc::QuantumComputation NeutralAtomMapper::map(qc::QuantumComputation& qc,
                                              InitialMapping initialMapping,
                                              InitialCoordinateMapping initialCoordinateMapping) {
  mappedQc = qc::QuantumComputation(arch.getNpositions());
  nMoves     = 0;
  nSwaps     = 0;
  nBridges   = 0;
  nFAncillas = 0;
  nFQubits   = 0;
  qc::CircuitOptimizer::replaceMCXWithMCZ(qc);
  qc::CircuitOptimizer::singleQubitGateFusion(qc);
  qc::CircuitOptimizer::flattenOperations(qc);
  qc::CircuitOptimizer::removeFinalMeasurements(qc);

  auto dag = qc::CircuitOptimizer::constructDAG(qc);
  qc::DAG dagUpdate = dag;

  // init mapping
  this->mapping = Mapping(arch.getNqubits(), initialMapping);
  
  // init coord mapping
  auto archNqubits    = arch.getNqubits();
  auto archNpositions = arch.getNpositions();
  std::vector<CoordIndex> qubitIndices(archNqubits, std::numeric_limits<unsigned int>::max());
  std::vector<CoordIndex> hwIndices(archNpositions, std::numeric_limits<unsigned int>::max());

  if(initialCoordinateMapping == Graph){
    graphMatching(qubitIndices, hwIndices, dag);
    this->hardwareQubits = HardwareQubits(arch, initialCoordinateMapping, qubitIndices, hwIndices, parameters.seed);
  }
  
  /*
  if(this->parameters.verbose){
    std::cout << "* Init Coord Mapping w/ [row:" << arch.getNrows() << " X col:" << arch.getNcolumns() << "] hardware"  << std::endl;
    for(uint32_t q=0; q<qc.getNqubits(); q++){
      std::cout << "q " << std::setw(3) << q;
      std::cout << " -> h " << std::setw(3) << this->mapping.getHwQubit(q);
      std::cout << " -> c " << std::setw(3) << this->hardwareQubits.getCoordIndex(q) << std::endl;
    }
    std::cout << std::endl;           
  }
  */

  // init layers
  NeutralAtomLayer frontLayer(dag);
  frontLayer.initLayerOffset();
  mapAllPossibleGates(frontLayer, dagUpdate);
  NeutralAtomLayer lookaheadLayer(dag);
  lookaheadLayer.initLayerOffset(frontLayer.getIteratorOffset());

  // Checks
  if (dag.size() > arch.getNqubits()) {
    throw std::runtime_error("More qubits in circuit than in architecture");
  }

  //   precompute exponential decay weights
  this->decayWeights.reserve(this->arch.getNcolumns());
  for (uint32_t i = this->arch.getNcolumns(); i > 0; --i) {
    this->decayWeights.emplace_back(std::exp(-this->parameters.decay * i));
  }

  auto i = 0;
  while (!frontLayer.getGates().empty()) {
    // assign gates to layers
    reassignGatesToLayers(frontLayer.getGates(), lookaheadLayer.getGates());
    if (this->parameters.verbose) {
      printLayers();
    }

    // save last swap to prevent immediate swap back
    Swap lastSwap = {0, 0};

    // first do all gate based mapping gates
    while (!this->frontLayerGate.empty()) {
      GateList gatesToExecute;
      while (gatesToExecute.empty() && !this->frontLayerGate.empty()) {
        ++i;
        if (this->parameters.verbose) {
          std::cout << "\niteration " << i << '\n';
        }

        // find swap
        auto bestSwap = findBestSwap(lastSwap);
        
        // find bridge
        auto allBridges = findAllBridges(qc, bestSwap);
        std::vector<std::pair<const qc::Operation*, Bridge>> ExecutableBridges;
        if(!allBridges.empty()){
          ExecutableBridges = compareCostBridgeSwap(allBridges, bestSwap, dagUpdate, frontLayer, qc);
        }

        // execute bridge gate
        if(!ExecutableBridges.empty()){
          updateMappingBridge(ExecutableBridges, frontLayer, lookaheadLayer, dagUpdate);
          mapAllPossibleGates(frontLayer, dagUpdate);
          lookaheadLayer.initLayerOffset(frontLayer.getIteratorOffset());
          reassignGatesToLayers(frontLayer.getGates(), lookaheadLayer.getGates());
          if (this->parameters.verbose) {
            printLayers();
          }
        }
        // execute swap gate
        else{
          lastSwap      = bestSwap;
          updateMappingSwap(bestSwap);
        }
        gatesToExecute = getExecutableGates(frontLayer.getGates());
      }
      mapAllPossibleGates(frontLayer, dagUpdate);
      lookaheadLayer.initLayerOffset(frontLayer.getIteratorOffset());
      reassignGatesToLayers(frontLayer.getGates(), lookaheadLayer.getGates());
      if (this->parameters.verbose) {
        printLayers();
      }
    }
    // then do all shuttling based mapping gates
    while (!this->frontLayerShuttling.empty()) {
      GateList gatesToExecute;
      while (gatesToExecute.empty() && !this->frontLayerShuttling.empty()) {
        ++i;
        if (this->parameters.verbose) {
          std::cout << "\niteration " << i << '\n';
        }

        // find best move
        auto [bestMove, bestMoveComb, opForMove] = findBestAtomMoveWithOp();

        // find flying ancilla
        auto bestFA = findBestFlyingAncilla(qc, opForMove);

        // find best flying qubit
        auto bestFQ = findBestFlyingQubit(qc, opForMove);

        // compare shuttling vs f.a.
        auto [useShuttling, useFlyingAncilla] = compareCostMoveFAandFQ(bestMoveComb, bestFA, bestFQ, opForMove, dagUpdate, frontLayer, qc);

        // execute shuttling
        if(useShuttling){
          updateMappingMove(bestMove);
        } 
        // execute flying ancilla / qubit
        else{
          if(useFlyingAncilla){
            updateMappingFlyingAncilla(bestFA, opForMove, frontLayer, lookaheadLayer);
          } else{
            updateMappingFlyingQubit(bestFQ, opForMove, frontLayer, lookaheadLayer);
          }
          mapAllPossibleGates(frontLayer, dagUpdate);
          lookaheadLayer.initLayerOffset(frontLayer.getIteratorOffset());
          reassignGatesToLayers(frontLayer.getGates(), lookaheadLayer.getGates());
          if(this->parameters.verbose){
            printLayers();
          }
        }

        gatesToExecute = getExecutableGates(frontLayer.getGates());
      }
      mapAllPossibleGates(frontLayer, dagUpdate);
      lookaheadLayer.initLayerOffset(frontLayer.getIteratorOffset());
      reassignGatesToLayers(frontLayer.getGates(), lookaheadLayer.getGates());
      if (this->parameters.verbose) {
        printLayers();
      }
    }
  }
  if (this->parameters.verbose) {
    std::cout << "nSwaps: "     << nSwaps     << '\n';
    std::cout << "nBridges: "   << nBridges   << '\n';
    std::cout << "nFAncillas: " << nFAncillas << '\n';
    std::cout << "nFQubits: "   << nFQubits   << '\n';
    std::cout << "nMoves: "     << nMoves     << '\n';
  }
  return mappedQc;
}
  
void NeutralAtomMapper::graphMatching(std::vector<CoordIndex>& qubitIndices, std::vector<CoordIndex>& hwIndices, const qc::DAG& dag){
  auto archNqubits    = arch.getNqubits();
  auto archNpositions = arch.getNpositions();
  auto archNrows      = arch.getNrows();
  auto archNcolumns   = arch.getNcolumns();
  // interaction graph
  std::vector<std::vector<double>> circGraph(dag.size(), std::vector<double> (dag.size(), 0.0) );
  std::vector<std::pair<int, std::pair<int, double>>> circGraph_degree(dag.size());
  std::vector<std::vector<std::pair<int,double>>> circGraph_neighbor(dag.size());
  for(uint32_t qubit=0; qubit<dag.size(); ++qubit){
    for(auto opPtr : dag[qubit]){
      auto* op = opPtr->get();
      auto usedQubits = op->getUsedQubits();
      if(usedQubits.size() > 1){
        for(auto i : usedQubits){
          if(i!=qubit){
            circGraph[qubit][i] += 1;
          }
        }
      }
    }
  }
  
  // generate graph matching queue
  for(uint32_t qubit=0; qubit<dag.size(); qubit++){
    int cnt = 0;
    double sum = 0;
    for(uint32_t i=0; i<dag.size(); i++){
      double weight = circGraph[qubit][i];
      if(weight > 0){
        cnt++;
        sum += weight;
        circGraph_neighbor[qubit].emplace_back( i, weight );
      }
    }
    circGraph_degree[qubit] = std::make_pair( qubit, std::make_pair(cnt, sum));
  }
  sort(circGraph_degree.begin(), circGraph_degree.end(), 
    [](std::pair<int, std::pair<int,double>> a, std::pair<int, std::pair<int, double>> b){
      if(a.second.first == b.second.first)
        return a.second.second > b.second.second;
      else
        return a.second.first > b.second.first;
    });
  std::queue<uint32_t> circGraph_queue;
  for(const auto i : circGraph_degree){
    circGraph_queue.push(i.first);
  }
  for(auto& innerVec : circGraph_neighbor){
    sort(innerVec.begin(), innerVec.end(), 
      [](std::pair<int,double>& a, std::pair<int,double>& b){
        return a.second > b.second;
      });
  }
  
  // graph matching
  bool firstCenter = true;
  uint32_t nMapped = 0;
  uint32_t archCenterX = (archNcolumns % 2 == 0) ? (archNcolumns/2 - 1) : (archNcolumns - 1)/2;
  uint32_t archCenterY = (archNrows % 2 == 0) ? (archNrows/2 - 1) : (archNrows - 1)/2;
  uint32_t archCenter = archCenterY * archNcolumns + archCenterX;
  
  while(!circGraph_queue.empty() && nMapped!=dag.size()){
    uint32_t qc = circGraph_queue.front();
    uint32_t hc = std::numeric_limits<unsigned int>::max(); //hardwarCenter
    // center mapping
    if(firstCenter){
      hc = archCenter;
      qubitIndices[qc] = hc;
      hwIndices[hc] = qc;
      firstCenter = false;
      nMapped++;
    }
    else if(qubitIndices[qc]==std::numeric_limits<unsigned int>::max()){
      
      //ref loc
      std::vector<int> refLoc;
      for(auto i : circGraph_neighbor[qc]){
        if(qubitIndices[i.first] != std::numeric_limits<unsigned int>::max()){
          refLoc.push_back( qubitIndices[i.first] );
        }
      }
      
      //candidate loc
      std::vector< std::pair<int,int> > distCandiLoc;
      std::vector<int> initCandiLoc;
      for(int i=0; i < archNpositions; i++){
        if(hwIndices[i] == std::numeric_limits<unsigned int>::max()){
          initCandiLoc.push_back(i);
        }
      }
      for(auto v : initCandiLoc){
        int distSum = 0;
        for(auto r : refLoc){
          int dist = std::abs(v%archNcolumns - r%archNcolumns) + std::abs(v/archNcolumns-r/archNcolumns);
          distSum += dist;
        }
        distCandiLoc.emplace_back(v, distSum);
      }
      sort(distCandiLoc.begin(), distCandiLoc.end(), 
        [](std::pair<int,int> a, std::pair<int,int> b){
          return a.second < b.second;
      });
      
      //find position
      hc = distCandiLoc[0].first;
      qubitIndices[qc] = hc;
      hwIndices[hc] = qc;
      nMapped++;
    }
    else{
      hc = qubitIndices[qc];
    }

    // neighbor mapping
    if(!circGraph_neighbor[qc].empty()){
      int idx_qc_n = 0;
      std::vector<int> qc_n;
      for(auto i : circGraph_neighbor[qc]){
        if( qubitIndices[i.first] != std::numeric_limits<unsigned int>::max() ) continue;
        if( idx_qc_n >= 4 ) continue;
        else{
          qc_n.push_back( i.first );
          idx_qc_n++;
        }
      }
      std::vector<int> hw_n;
      if(hc!=std::numeric_limits<unsigned int>::max() && qc_n.size()>0){
        if((hc+1)%archNcolumns != 0 && hwIndices[hc+1]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc+1); //right
        if(hc/archNcolumns < (archNrows-1) && hwIndices[hc+archNcolumns]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc+archNcolumns);  //down
        if(hc%archNcolumns != 0 && hwIndices[hc-1]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc-1); //left
        if(hc>archNcolumns && hwIndices[hc-archNcolumns]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc-archNcolumns); //up
      }

      int minSize = std::min( qc_n.size(), hw_n.size() );
      for(int i=0; i<minSize; i++){
        int qc_i = qc_n[i];
        int hw_i = hw_n[i];
        qubitIndices[qc_i] = hw_i;
        hwIndices[hw_i] = qc_i;
        nMapped++;
      }

    }
    circGraph_queue.pop();
  }
}


void NeutralAtomMapper::mapAllPossibleGates(NeutralAtomLayer& layer, qc::DAG& dag) {
  // map single qubit gates
  for (const auto* opPointer : layer.getMappedSingleQubitGates()) {
    mapGate(opPointer);
    // update dag
    auto qubits = opPointer->getUsedQubits();
    auto q = *qubits.begin();
    auto it = std::find_if( dag[q].begin(), dag[q].end(),
                            [opPointer](const std::unique_ptr<qc::Operation>* op){
                                return op->get() == opPointer;
                                });
    if(it!=dag[q].end()){
      dag[q].erase(it);
    }
  }
  layer.removeGatesAndUpdate({});
  // check and map multi qubit gates
  auto executableGates = getExecutableGates(layer.getGates());
  while (!executableGates.empty()) {
    for (const auto* opPointer : layer.getMappedSingleQubitGates()) {
      mapGate(opPointer);
      // update dag
      auto qubits = opPointer->getUsedQubits();
      auto q = *qubits.begin();
      auto it = std::find_if( dag[q].begin(), dag[q].end(),
                              [opPointer](const std::unique_ptr<qc::Operation>* op){
                                  return op->get() == opPointer;
                                  });
      if(it!=dag[q].end()){
        dag[q].erase(it);
      }
    }
    for (const auto* opPointer : executableGates) {
      mapGate(opPointer);
      // update dag
      auto qubits = opPointer->getUsedQubits();
      for(auto q : qubits){
        auto it = std::find_if( dag[q].begin(), dag[q].end(),
                                [opPointer](const std::unique_ptr<qc::Operation>* op){
                                    return op->get() == opPointer;
                                    });
        if(it!=dag[q].end()){
          dag[q].erase(it);
        }
      }
    }
    layer.removeGatesAndUpdate(executableGates);
    executableGates = getExecutableGates(layer.getGates());
  }
}

qc::QuantumComputation
NeutralAtomMapper::convertToAod(qc::QuantumComputation& qc) {
  // decompose SWAP gates
  qc::CircuitOptimizer::decomposeSWAP(qc, false);
  qc::CircuitOptimizer::replaceMCXWithMCZ(qc);
  qc::CircuitOptimizer::singleQubitGateFusion(qc);
  qc::CircuitOptimizer::flattenOperations(qc);
  // decompose AOD moves
  //MoveToAodConverter aodScheduler(arch);
  MoveToAodConverter aodScheduler(arch, this->hardwareQubits);
  mappedQcAOD = aodScheduler.schedule(qc); 
  if (this->parameters.verbose) {
    std::cout << "nMoveGroups: " << aodScheduler.getNMoveGroups() << '\n';
  }
  return mappedQcAOD;
}

void NeutralAtomMapper::reassignGatesToLayers(const GateList& frontGates,
                                              const GateList& lookaheadGates) {
  // assign gates to gates or shuttling
  this->frontLayerGate.clear();
  this->frontLayerShuttling.clear();
  for (const auto& gate : frontGates) {
    if (swapGateBetter(gate)) {
      this->frontLayerGate.emplace_back(gate);
    } else {
      this->frontLayerShuttling.emplace_back(gate);
    }
  }

  this->lookaheadLayerGate.clear();
  this->lookaheadLayerShuttling.clear();
  for (const auto& gate : lookaheadGates) {
    if (swapGateBetter(gate)) {
      this->lookaheadLayerGate.emplace_back(gate);
    } else {
      this->lookaheadLayerShuttling.emplace_back(gate);
    }
  }
}

void NeutralAtomMapper::mapGate(const qc::Operation* op) {
  if (op->getType() == qc::OpType::I) {
    return;
  }
  // Safety check
  if (std::find(this->executedCommutingGates.begin(),
                this->executedCommutingGates.end(),
                op) != this->executedCommutingGates.end()) {
    return;
  }
  this->executedCommutingGates.emplace_back(op);
  if (this->parameters.verbose) {
    std::cout << "mapped " << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << "\n";
  }
  // convert circuit qubits to CoordIndex and append to mappedQc
  auto  opCopyUnique = op->clone();
  auto* opCopy       = opCopyUnique.get();
  this->mapping.mapToHwQubits(opCopy);
  this->hardwareQubits.mapToCoordIdx(opCopy);
  this->mappedQc.emplace_back(opCopy->clone());
}

bool NeutralAtomMapper::isExecutable(const qc::Operation* opPointer) {
  auto usedQubits  = opPointer->getUsedQubits();
  auto nUsedQubits = usedQubits.size();
  if (nUsedQubits == 1) {
    return true;
  }
  std::set<qc::Qubit> usedHwQubits;
  for (auto qubit : usedQubits) {
    usedHwQubits.emplace(this->mapping.getHwQubit(qubit));
  }
  return this->hardwareQubits.getAllToAllSwapDistance(usedHwQubits) == 0;
}

void NeutralAtomMapper::printLayers() {
  std::cout << "f,g: ";
  for (const auto* op : this->frontLayerGate) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << "f,s: ";
  for (const auto* op : this->frontLayerShuttling) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << "l,g: ";
  for (const auto* op : this->lookaheadLayerGate) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
  std::cout << "l,s: ";
  for (const auto* op : this->lookaheadLayerShuttling) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

GateList NeutralAtomMapper::getExecutableGates(const GateList& gates) {
  GateList executableGates;
  for (const auto* opPointer : gates) {
    if (isExecutable(opPointer)) {
      executableGates.emplace_back(opPointer);
    }
  }
  return executableGates;
}

void NeutralAtomMapper::updateMappingSwap(Swap swap) {
  nSwaps++;
  // save to lastSwaps
  this->lastBlockedQubits.emplace_back(
      this->hardwareQubits.getBlockedQubits({swap.first, swap.second}));
  if (this->lastBlockedQubits.size() > this->arch.getNcolumns()) {
    this->lastBlockedQubits.pop_front();
  }
  this->mapping.applySwap(swap);
  // convert circuit qubits to CoordIndex and append to mappedQc
  auto idxFirst  = this->hardwareQubits.getCoordIndex(swap.first);
  auto idxSecond = this->hardwareQubits.getCoordIndex(swap.second);
  this->mappedQc.swap(idxFirst, idxSecond);
  if (this->parameters.verbose) {
    std::cout << "swapped " << swap.first << " " << swap.second;
    std::cout << "  logical qubits: ";
    if (this->mapping.isMapped(swap.first)) {
      std::cout << this->mapping.getCircQubit(swap.first);
    } else {
      std::cout << "not mapped";
    }
    if (this->mapping.isMapped(swap.second)) {
      std::cout << " " << this->mapping.getCircQubit(swap.second);
    } else {
      std::cout << " not mapped";
    }
    std::cout << '\n';
  }
}

void NeutralAtomMapper::updateMappingMove(AtomMove move) {
  this->lastMoves.emplace_back(move);
  if (this->lastMoves.size() > 4) {
    this->lastMoves.pop_front();
  }
  mappedQc.move(move.first, move.second);
  auto toMoveHwQubit = this->hardwareQubits.getHwQubit(move.first);
  this->hardwareQubits.move(toMoveHwQubit, move.second);
  if (this->parameters.verbose) {
    std::cout << "moved " << move.first << " to " << move.second;
    if (this->mapping.isMapped(toMoveHwQubit)) {
      std::cout << "  logical qubit: "
                << this->mapping.getCircQubit(toMoveHwQubit) << '\n';
    } else {
      std::cout << "  not mapped" << '\n';
    }
  }
  nMoves++;
}

Swap NeutralAtomMapper::findBestSwap(const Swap& lastSwap) {
  // compute necessary movements
  auto swapsFront     = initSwaps(this->frontLayerGate);
  auto swapsLookahead = initSwaps(this->lookaheadLayerGate);
  setTwoQubitSwapWeight(swapsFront.second);

  // evaluate swaps based on cost function
  auto swaps = getAllPossibleSwaps(swapsFront);
  // remove last swap to prevent immediate swap back
  swaps.erase(lastSwap);
  swaps.erase({lastSwap.second, lastSwap.first});

  // no swap possible
  if (swaps.empty()) {
    return {std::numeric_limits<qc::Qubit>::max(),
            std::numeric_limits<qc::Qubit>::max()};
  }
  std::vector<std::pair<Swap, qc::fp>> swapCosts;
  swapCosts.reserve(swaps.size());
  for (const auto& swap : swaps) {
    swapCosts.emplace_back(swap, swapCost(swap, swapsFront, swapsLookahead));
  }
  std::sort(swapCosts.begin(), swapCosts.end(),
            [](const auto& swap1, const auto& swap2) {
              return swap1.second < swap2.second;
            });
  // get swap of minimal cost
  auto bestSwap = std::min_element(swapCosts.begin(), swapCosts.end(),
                                   [](const auto& swap1, const auto& swap2) {
                                     return swap1.second < swap2.second;
                                   });
  return bestSwap->first;
}

void NeutralAtomMapper::setTwoQubitSwapWeight(const WeightedSwaps& swapExact) {
  for (const auto& [swap, weight] : swapExact) {
    this->twoQubitSwapWeight = std::min(weight, this->twoQubitSwapWeight);
  }
}

std::set<Swap> NeutralAtomMapper::getAllPossibleSwaps(
    const std::pair<Swaps, WeightedSwaps>& swapsFront) const {
  auto [swapCloseByFront, swapExactFront] = swapsFront;
  std::set<Swap> swaps;
  for (const auto& swapNearby : swapCloseByFront) {
    const auto nearbySwapsFirst =
        this->hardwareQubits.getNearbySwaps(swapNearby.first);
    for (const auto& swapFirst : nearbySwapsFirst) {
      swaps.emplace(swapFirst);
    }
    const auto nearbySwapsSecond =
        this->hardwareQubits.getNearbySwaps(swapNearby.second);
    for (const auto& swapSecond : nearbySwapsSecond) {
      swaps.emplace(swapSecond);
    }
  }
  for (const auto& [swap, weight] : swapExactFront) {
    const auto nearbySwapsFirst =
        this->hardwareQubits.getNearbySwaps(swap.first);
    for (const auto& swapFirst : nearbySwapsFirst) {
      swaps.emplace(swapFirst);
    }
  }
  return swaps;
}

qc::fp NeutralAtomMapper::swapCost(
    const Swap& swap, const std::pair<Swaps, WeightedSwaps>& swapsFront,
    const std::pair<Swaps, WeightedSwaps>& swapsLookahead) {
  auto [swapCloseByFront, swapExactFront]         = swapsFront;
  auto [swapCloseByLookahead, swapExactLookahead] = swapsLookahead;
  // compute the change in total distance
  auto distanceChangeFront =
      swapCostPerLayer(swap, swapCloseByFront, swapExactFront) /
      static_cast<qc::fp>(this->frontLayerGate.size());
  qc::fp distanceChangeLookahead = 0;
  if (!this->lookaheadLayerGate.empty()) {
    distanceChangeLookahead =
        swapCostPerLayer(swap, swapCloseByLookahead, swapExactLookahead) /
        static_cast<qc::fp>(this->lookaheadLayerGate.size());
  }
  auto cost = parameters.lookaheadWeightSwaps * distanceChangeLookahead +
              distanceChangeFront;
  //  compute the last time one of the swap qubits was used
  if (this->parameters.decay != 0) {
    uint32_t idxLastUsed = 0;
    for (uint32_t i = 0; i < this->lastBlockedQubits.size(); ++i) {
      if (this->lastBlockedQubits[i].find(swap.first) !=
              this->lastBlockedQubits[i].end() ||
          this->lastBlockedQubits[i].find(swap.second) !=
              this->lastBlockedQubits[i].end()) {
        idxLastUsed = i;
        break;
      }
    }
    cost *= this->decayWeights[idxLastUsed];
  }
  return cost;
}

std::pair<Swaps, WeightedSwaps>
NeutralAtomMapper::initSwaps(const GateList& layer) {
  Swaps         swapCloseBy = {};
  WeightedSwaps swapExact   = {};
  // computes for each gate the necessary moves to execute it
  for (const auto& gate : layer) {
    auto usedQubits   = gate->getUsedQubits();
    auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
    if (usedQubits.size() == 2) {
      // swap close by for two qubit gates
      swapCloseBy.emplace_back(*usedHwQubits.begin(), *usedHwQubits.rbegin());
    } else {
      // for multi-qubit gates, find the best position around the gate qubits
      auto bestPos = getBestMultiQubitPosition(gate);
      if (this->parameters.verbose) {
        std::cout << "bestPos: ";
        for (auto qubit : bestPos) {
          std::cout << qubit << " ";
        }
        std::cout << '\n';
      }
      // then compute the exact moves to get to the best position
      auto exactSwapsToPos = getExactSwapsToPosition(gate, bestPos);
      swapExact.insert(swapExact.end(), exactSwapsToPos.begin(),
                       exactSwapsToPos.end());
    }
  }
  // sort and remove duplicates from moveExact
  swapExact.erase(std::unique(swapExact.begin(), swapExact.end(),
                              [](const auto& swap1, const auto& swap2) {
                                return swap1.first == swap2.first;
                              }),
                  swapExact.end());
  return {swapCloseBy, swapExact};
}

qc::fp NeutralAtomMapper::swapCostPerLayer(const Swap&          swap,
                                           const Swaps&         swapCloseBy,
                                           const WeightedSwaps& swapExact) {
  SwapDistance distBefore = 0;
  SwapDistance distAfter  = 0;
  qc::fp       distChange = 0;
  // bring close only until swap distance =0, bring exact to the exact position
  // bring qubits together to execute gate
  for (const auto& [q1, q2] : swapCloseBy) {
    // distance before
    distBefore = this->hardwareQubits.getSwapDistance(q1, q2);
    if (distBefore == std::numeric_limits<SwapDistance>::max()) {
      continue;
    }
    // do swap
    if (q1 == swap.first) {
      distAfter = this->hardwareQubits.getSwapDistance(swap.second, q2);
    } else if (q2 == swap.second) {
      distAfter = this->hardwareQubits.getSwapDistance(q1, swap.first);
    } else if (q1 == swap.second) {
      distAfter = this->hardwareQubits.getSwapDistance(swap.first, q2);
    } else if (q2 == swap.first) {
      distAfter = this->hardwareQubits.getSwapDistance(q1, swap.second);
    } else {
      continue;
    }
    distChange +=
        static_cast<qc::fp>(distAfter - distBefore) * this->twoQubitSwapWeight;
  }

  // move qubits to the exact position for multi-qubit gates
  for (const auto& [exactSwap, weight] : swapExact) {
    auto origin      = exactSwap.first;
    auto destination = exactSwap.second;
    distBefore =
        this->hardwareQubits.getSwapDistance(origin, destination, false);
    if (distBefore == std::numeric_limits<SwapDistance>::max()) {
      continue;
    }
    if (origin == swap.first) {
      if (destination == swap.second) {
        distAfter = 0;
      } else {
        distAfter = this->hardwareQubits.getSwapDistance(swap.second,
                                                         destination, false);
      }
    } else if (origin == swap.second) {
      if (destination == swap.first) {
        distAfter = 0;
      } else {
        distAfter = this->hardwareQubits.getSwapDistance(swap.first,
                                                         destination, false);
      }
    } else {
      continue;
    }
    // multiply by multi-qubit weight
    // is larger for more qubits and if the qubits are closer together
    distChange += static_cast<qc::fp>(distAfter - distBefore) * weight;
  }

  return distChange;
}

HwQubits NeutralAtomMapper::getBestMultiQubitPosition(const qc::Operation* op) {
  // try to find position around gate Qubits recursively
  // if not, search through coupling graph until found according to a
  // priority queue based on the distance to the other qubits

  std::priority_queue<std::pair<qc::fp, HwQubit>,
                      std::vector<std::pair<qc::fp, HwQubit>>, std::greater<>>
      qubitQueue;
  // add the gate qubits to the priority queue
  auto gateQubits   = op->getUsedQubits();
  auto gateHwQubits = this->mapping.getHwQubits(gateQubits);
  // add the gate qubits to the priority queue
  for (const auto& gateQubit : gateHwQubits) {
    qc::fp totalDist = 0;
    for (const auto& otherGateQubit : gateHwQubits) {
      if (gateQubit == otherGateQubit) {
        continue;
      }
      totalDist +=
          this->hardwareQubits.getSwapDistance(gateQubit, otherGateQubit, true);
    }
    qubitQueue.emplace(totalDist, gateQubit);
  }

  // run through the priority queue until a position is found
  std::set<HwQubit> visitedQubits;
  while (!qubitQueue.empty()) {
    auto qubit = qubitQueue.top().second;
    visitedQubits.emplace(qubit);
    qubitQueue.pop();

    // remove selected qubit from the gate qubits
    auto tempGateHwQubits = gateHwQubits;
    if (tempGateHwQubits.find(qubit) != tempGateHwQubits.end()) {
      tempGateHwQubits.erase(tempGateHwQubits.find(qubit));
    }
    auto bestPos = getBestMultiQubitPositionRec(
        tempGateHwQubits, {qubit}, this->hardwareQubits.getNearbyQubits(qubit));
    if (!bestPos.empty()) {
      return bestPos;
    }
    // add nearby qubits to the priority queue
    for (const auto& nearbyQubit :
         this->hardwareQubits.getNearbyQubits(qubit)) {
      if (visitedQubits.find(nearbyQubit) != visitedQubits.end()) {
        continue;
      }
      // compute total distance to all other gate qubits
      qc::fp totalDist = 0;
      for (const auto& otherGateQubit : gateHwQubits) {
        if (nearbyQubit == otherGateQubit) {
          continue;
        }
        totalDist +=
            this->hardwareQubits.getSwapDistance(nearbyQubit, otherGateQubit);
      }
      qubitQueue.emplace(totalDist, nearbyQubit);
    }
  }
  // find gate and move it to the shuttling layer
  auto idxFrontGate =
      std::find(this->frontLayerGate.begin(), this->frontLayerGate.end(), op);
  if (idxFrontGate != this->frontLayerGate.end()) {
    this->frontLayerGate.erase(idxFrontGate);
    this->frontLayerShuttling.emplace_back(op);
  }
  // remove from lookahead layer if there
  auto idxLookaheadGate = std::find(this->lookaheadLayerGate.begin(),
                                    this->lookaheadLayerGate.end(), op);
  if (idxLookaheadGate != this->lookaheadLayerGate.end()) {
    this->lookaheadLayerGate.erase(idxLookaheadGate);
    this->lookaheadLayerShuttling.emplace_back(op);
  }
  return {};
}

HwQubits NeutralAtomMapper::getBestMultiQubitPositionRec(
    HwQubits remainingGateQubits, std::vector<HwQubit> selectedQubits,
    HwQubits remainingNearbyQubits) {
  // check if done
  if (remainingGateQubits.empty()) {
    HwQubits bestPos(selectedQubits.begin(), selectedQubits.end());
    return bestPos;
  }
  // update remainingNearbyQubits
  auto newQubit        = *selectedQubits.rbegin();
  auto nearbyNextQubit = this->hardwareQubits.getNearbyQubits(newQubit);
  // compute remaining qubits as the intersection with current
  Qubits newRemainingQubits;
  std::set_intersection(
      remainingNearbyQubits.begin(), remainingNearbyQubits.end(),
      nearbyNextQubit.begin(), nearbyNextQubit.end(),
      std::inserter(newRemainingQubits, newRemainingQubits.begin()));
  for (const auto& qubit : selectedQubits) {
    if (newRemainingQubits.find(qubit) != newRemainingQubits.end()) {
      newRemainingQubits.erase(newRemainingQubits.find(qubit));
    }
  }
  remainingNearbyQubits = newRemainingQubits;

  // if not enough space
  if (remainingNearbyQubits.size() < remainingGateQubits.size()) {
    return {};
  }

  std::vector<std::pair<HwQubit, qc::fp>> summedDistances;
  for (const auto& hwQubit : remainingNearbyQubits) {
    qc::fp distance = 0;
    for (const auto& gateHwQubit : remainingGateQubits) {
      if (hwQubit == gateHwQubit) {
        // gate qubit is already at one of the positions -> assign it
        selectedQubits.emplace_back(hwQubit);
        remainingGateQubits.erase(remainingGateQubits.find(gateHwQubit));
        return getBestMultiQubitPositionRec(remainingGateQubits, selectedQubits,
                                            remainingNearbyQubits);
      }
      distance +=
          this->hardwareQubits.getSwapDistance(hwQubit, gateHwQubit, true);
    }
    summedDistances.emplace_back(hwQubit, distance);
  }
  // select next qubit as the one with minimal distance
  auto nextQubitDist =
      std::min_element(summedDistances.begin(), summedDistances.end(),
                       [](const auto& qubit1, const auto& qubit2) {
                         return qubit1.second < qubit2.second;
                       });
  auto nextQubit = nextQubitDist->first;
  selectedQubits.emplace_back(nextQubit);
  // remove from remaining gate qubits the one that is closest to the next
  auto closesGateQubits = *remainingGateQubits.begin();
  auto closesDistance =
      this->hardwareQubits.getSwapDistance(closesGateQubits, nextQubit, true);
  for (const auto& gateQubit : remainingGateQubits) {
    auto distance =
        this->hardwareQubits.getSwapDistance(gateQubit, nextQubit, true);
    if (distance < closesDistance) {
      closesGateQubits = gateQubit;
      closesDistance   = distance;
    }
  }
  remainingGateQubits.erase(remainingGateQubits.find(closesGateQubits));

  return getBestMultiQubitPositionRec(remainingGateQubits, selectedQubits,
                                      remainingNearbyQubits);
}

WeightedSwaps
NeutralAtomMapper::getExactSwapsToPosition(const qc::Operation* op,
                                           HwQubits             position) {
  if (position.empty()) {
    return {};
  }
  auto          gateQubits   = op->getUsedQubits();
  auto          gateHwQubits = this->mapping.getHwQubits(gateQubits);
  WeightedSwaps swapsExact;
  while (!position.empty() && !gateHwQubits.empty()) {
    std::vector<std::tuple<HwQubit, std::set<HwQubit>, SwapDistance>>
                      minimalDistances;
    std::set<HwQubit> minimalDistancePosQubit;
    for (const auto& gateQubit : gateHwQubits) {
      SwapDistance minimalDistance = std::numeric_limits<SwapDistance>::max();
      for (const auto& posQubit : position) {
        auto distance =
            this->hardwareQubits.getSwapDistance(gateQubit, posQubit, false);
        if (distance < minimalDistance) {
          minimalDistance = distance;
          minimalDistancePosQubit.clear();
          minimalDistancePosQubit.emplace(posQubit);
        } else if (distance == minimalDistance) {
          minimalDistancePosQubit.emplace(posQubit);
        }
      }
      if (minimalDistance == std::numeric_limits<SwapDistance>::max()) {
        // not possible to move to position
        // move gate to shuttling layer
        auto idxFrontGate = std::find(this->frontLayerGate.begin(),
                                      this->frontLayerGate.end(), op);
        if (idxFrontGate != this->frontLayerGate.end()) {
          this->frontLayerGate.erase(idxFrontGate);
          this->frontLayerShuttling.emplace_back(op);
        }
        // remove from lookahead layer if there
        auto idxLookaheadGate = std::find(this->lookaheadLayerGate.begin(),
                                          this->lookaheadLayerGate.end(), op);
        if (idxLookaheadGate != this->lookaheadLayerGate.end()) {
          this->lookaheadLayerGate.erase(idxLookaheadGate);
          this->lookaheadLayerShuttling.emplace_back(op);
        }
        return {};
      }
      minimalDistances.emplace_back(gateQubit, minimalDistancePosQubit,
                                    minimalDistance);
    }
    // find gate qubit with maximal minimal distance to assign first to a
    // position
    auto assignFirst =
        std::max_element(minimalDistances.begin(), minimalDistances.end(),
                         [](const auto& qubit1, const auto& qubit2) {
                           return std::get<2>(qubit1) < std::get<2>(qubit2);
                         });

    auto assignedGateQubit = std::get<0>(*assignFirst);
    auto assignedPosQubits = std::get<1>(*assignFirst);
    // for multiple equal good positions, choose the one that
    // is not assigned to one of the other ones
    HwQubit assignedPosQubit = *assignedPosQubits.begin();
    if (assignedPosQubits.size() > 1) {
      for (const auto& posQubit : assignedPosQubits) {
        // as all places within the position can reach each other, it is
        // sufficient to check for a single unoccupied position
        // check if posQubit is assigned at its current position
        if (std::none_of(minimalDistances.begin(), minimalDistances.end(),
                         [&posQubit](const auto& qubit) {
                           return std::get<0>(qubit) == posQubit &&
                                  *(std::get<1>(qubit).begin()) == posQubit;
                         })) {
          assignedPosQubit = posQubit;
          break;
        }
      }
    }

    // assign gateQubit to position by removing both from gateHwQubits and
    // position
    gateHwQubits.erase(gateHwQubits.find(assignedGateQubit));
    position.erase(position.find(assignedPosQubit));
    // and add to exactMove if not swap with one of the other qubits
    // only problem if their exact swap distance is 1
    if (std::none_of(gateHwQubits.begin(), gateHwQubits.end(),
                     [&assignedGateQubit, this](const auto& qubit) {
                       return assignedGateQubit == qubit &&
                              this->hardwareQubits.getSwapDistance(
                                  assignedGateQubit, qubit, false) == 1;
                     }) &&
        assignedGateQubit != assignedPosQubit) {
      swapsExact.emplace_back(
          std::make_pair(assignedGateQubit, assignedPosQubit), 0);
    }
  }

  // compute total distance of all moves
  SwapDistance totalDistance = 0;
  for (const auto& [swap, weight] : swapsExact) {
    auto [q1, q2] = swap;
    totalDistance += this->hardwareQubits.getSwapDistance(q1, q2, false);
  }
  // add cost to the moves -> move first qubit corresponding to almost finished
  // positions
  auto nQubits = op->getUsedQubits().size();
  auto multiQubitFactor =
      (static_cast<qc::fp>(nQubits) * static_cast<qc::fp>(nQubits - 1)) / 2;
  for (auto& move : swapsExact) {
    move.second = multiQubitFactor / static_cast<qc::fp>(totalDistance);
  }

  return swapsExact;
}

AtomMove NeutralAtomMapper::findBestAtomMove() {
  auto moveCombs = getAllMoveCombinations();

  // compute cost for each move combination
  std::vector<std::pair<MoveComb, qc::fp>> moveCosts;
  moveCosts.reserve(moveCombs.size());
  for (const auto& moveComb : moveCombs) {
    moveCosts.emplace_back(moveComb, moveCostComb(moveComb));
  }

  std::sort(moveCosts.begin(), moveCosts.end(),
            [](const auto& move1, const auto& move2) {
              return move1.second < move2.second;
            });

  // get move of minimal cost
  auto bestMove = std::min_element(moveCosts.begin(), moveCosts.end(),
                                   [](const auto& move1, const auto& move2) {
                                     return move1.second < move2.second;
                                   });
  return bestMove->first.getFirstMove();
}

std::tuple<AtomMove, MoveComb, const qc::Operation*> NeutralAtomMapper::findBestAtomMoveWithOp() {
  auto [moveCombs, moveCombsWithOp] = getAllMoveCombinationsWithOp();

  // compute cost for each move combination
  std::vector<std::pair<MoveComb, qc::fp>> moveCosts;
  moveCosts.reserve(moveCombs.size());
  for (const auto& moveComb : moveCombs) {
    moveCosts.emplace_back(moveComb, moveCostComb(moveComb));
  }

  std::sort(moveCosts.begin(), moveCosts.end(),
            [](const auto& move1, const auto& move2) {
              return move1.second < move2.second;
            });

  // get move of minimal cost
  auto bestMove = std::min_element(moveCosts.begin(), moveCosts.end(),
                                   [](const auto& move1, const auto& move2) {
                                     return move1.second < move2.second;
                                   });
  MoveComb bestAtomMove = bestMove->first;
  const qc::Operation* corresOp = nullptr;
  for(const auto& pair : moveCombsWithOp){
    if(pair.first == bestAtomMove){
      corresOp = pair.second;
      break;
    }
  }
  return make_tuple(bestAtomMove.getFirstMove(), bestAtomMove, corresOp);
}

qc::fp NeutralAtomMapper::moveCostComb(const MoveComb& moveComb) {
  qc::fp costComb = 0;
  for (const auto& move : moveComb.moves) {
    costComb += moveCost(move);
  }
  return costComb;
}

qc::fp NeutralAtomMapper::moveCost(const AtomMove& move) {
  qc::fp cost      = 0;
  auto   frontCost = moveCostPerLayer(move, this->frontLayerShuttling) /
                   static_cast<qc::fp>(this->frontLayerShuttling.size());
  cost += frontCost;
  if (!lookaheadLayerShuttling.empty()) {
    auto lookaheadCost =
        moveCostPerLayer(move, this->lookaheadLayerShuttling) /
        static_cast<qc::fp>(this->lookaheadLayerShuttling.size());
    cost += parameters.lookaheadWeightMoves * lookaheadCost;
  }
  if (!this->lastMoves.empty()) {
    auto parallelCost = parameters.shuttlingTimeWeight *
                        parallelMoveCost(move) /
                        static_cast<qc::fp>(this->lastMoves.size()) /
                        static_cast<qc::fp>(this->frontLayerShuttling.size());
    cost += parallelCost;
  }

  return cost;
}

qc::fp NeutralAtomMapper::moveCostPerLayer(const AtomMove& move,
                                           GateList&       layer) {
  // compute cost assuming the move was applied
  qc::fp distChange    = 0;
  auto   toMoveHwQubit = this->hardwareQubits.getHwQubit(move.first);
  if (this->mapping.isMapped(toMoveHwQubit)) {
    auto toMoveCircuitQubit = this->mapping.getCircQubit(toMoveHwQubit);
    for (const auto& gate : layer) {
      auto usedQubits = gate->getUsedQubits();
      if (usedQubits.find(toMoveCircuitQubit) != usedQubits.end()) {
        // check distance reduction
        qc::fp distanceBefore = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          auto dist    = this->arch.getEuclideanDistance(
              this->hardwareQubits.getCoordIndex(hwQubit),
              this->hardwareQubits.getCoordIndex(toMoveHwQubit));
          distanceBefore += dist;
        }
        qc::fp distanceAfter = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          auto dist    = this->arch.getEuclideanDistance(
              this->hardwareQubits.getCoordIndex(hwQubit), move.second);
          distanceAfter += dist;
        }
        distChange += distanceAfter - distanceBefore;
      }
    }
  }
  return distChange;
}

qc::fp NeutralAtomMapper::parallelMoveCost(const AtomMove& move) {
  qc::fp parallelCost = 0;
  auto   moveVector   = this->arch.getVector(move.first, move.second);
  std::vector<CoordIndex> lastEndingCoords;
  if (this->lastMoves.empty()) {
    parallelCost += arch.getVectorShuttlingTime(moveVector);
  }
  for (const auto& lastMove : this->lastMoves) {
    lastEndingCoords.emplace_back(lastMove.second);
    // decide of shuttling can be done in parallel
    auto lastMoveVector = this->arch.getVector(lastMove.first, lastMove.second);
    if (moveVector.overlap(lastMoveVector)) {
      if (moveVector.direction != lastMoveVector.direction) {
        parallelCost += arch.getVectorShuttlingTime(moveVector);
      } else {
        // check if move can be done in parallel
        if (moveVector.include(lastMoveVector)) {
          parallelCost += arch.getVectorShuttlingTime(moveVector);
        }
      }
    }
  }
  // check if in same row/column like last moves
  // then can may be loaded in parallel
  auto moveCoordInit = this->arch.getCoordinate(move.first);
  auto moveCoordEnd  = this->arch.getCoordinate(move.second);
  parallelCost += arch.getShuttlingTime(qc::OpType::AodActivate) +
                  arch.getShuttlingTime(qc::OpType::AodDeactivate);
  for (const auto& lastMove : this->lastMoves) {
    auto lastMoveCoordInit = this->arch.getCoordinate(lastMove.first);
    auto lastMoveCoordEnd  = this->arch.getCoordinate(lastMove.second);
    if (moveCoordInit.x == lastMoveCoordInit.x ||
        moveCoordInit.y == lastMoveCoordInit.y) {
      parallelCost -= arch.getShuttlingTime(qc::OpType::AodActivate);
    }
    if (moveCoordEnd.x == lastMoveCoordEnd.x ||
        moveCoordEnd.y == lastMoveCoordEnd.y) {
      parallelCost -= arch.getShuttlingTime(qc::OpType::AodDeactivate);
    }
  }
  // check if move can use AOD atom from last moves
  //  if (std::find(lastEndingCoords.begin(), lastEndingCoords.end(),
  //  move.first) ==
  //      lastEndingCoords.end()) {
  //    parallelCost += arch.getShuttlingTime(qc::OpType::AodActivate) +
  //                    arch.getShuttlingTime(qc::OpType::AodDeactivate);
  //  }
  return parallelCost;
}

MultiQubitMovePos
NeutralAtomMapper::getMovePositionRec(MultiQubitMovePos   currentPos,
                                      const CoordIndices& gateCoords,
                                      const size_t&       maxNMoves) {
  if (currentPos.coords.size() == gateCoords.size()) {
    return currentPos;
  }
  if (currentPos.nMoves > maxNMoves) {
    return {};
  }

  auto nearbyCoords = this->arch.getNearbyCoordinates(currentPos.coords.back());
  // filter out coords that have a SWAP distance unequal to 0 to any of the
  // current qubits. Also sort out coords that are already in the vector
  std::vector<CoordIndex> filteredNearbyCoords;
  for (const auto& coord : nearbyCoords) {
    bool valid = true;
    for (const auto& qubit : currentPos.coords) {
      if (this->arch.getSwapDistance(qubit, coord) != 0 || coord == qubit) {
        valid = false;
        break;
      }
    }
    if (valid) {
      filteredNearbyCoords.emplace_back(coord);
    }
  }

  // differentiate between free and occupied coords
  CoordIndices freeNearbyCoords;
  CoordIndices occupiedNearbyCoords;
  CoordIndices occupiedGateCoords;
  for (const auto& coord : filteredNearbyCoords) {
    if (this->hardwareQubits.isMapped(coord)) {
      if (std::find(gateCoords.begin(), gateCoords.end(), coord) !=
          gateCoords.end()) {
        occupiedGateCoords.emplace_back(coord);
      } else {
        occupiedNearbyCoords.emplace_back(coord);
      }
    } else {
      freeNearbyCoords.emplace_back(coord);
    }
  }

  // compute minimal possible moves
  size_t       minPossibleMoves = currentPos.nMoves;
  size_t const nMissingQubits   = gateCoords.size() - currentPos.coords.size();
  auto         itGate           = occupiedGateCoords.begin();
  auto         itFree           = freeNearbyCoords.begin();
  auto         itOcc            = occupiedNearbyCoords.begin();
  for (size_t i = 0; i < nMissingQubits; ++i) {
    if (itGate != occupiedGateCoords.end()) {
      ++itGate;
    } else if (itFree != freeNearbyCoords.end()) {
      ++itFree;
      minPossibleMoves += 1;
    } else if (itOcc != occupiedNearbyCoords.end()) {
      ++itOcc;
      minPossibleMoves += 2;
    }
  }
  if (minPossibleMoves > maxNMoves) {
    return {};
  }

  for (const auto& gateCoord : occupiedGateCoords) {
    MultiQubitMovePos nextPos = MultiQubitMovePos(currentPos);
    nextPos.coords.emplace_back(gateCoord);
    auto bestPos = getMovePositionRec(nextPos, gateCoords, maxNMoves);
    if (bestPos.coords.size() == gateCoords.size()) {
      return bestPos;
    }
  }

  for (const auto& freeCoord : freeNearbyCoords) {
    MultiQubitMovePos nextPos = MultiQubitMovePos(currentPos);
    nextPos.coords.emplace_back(freeCoord);
    nextPos.nMoves += 1;
    auto bestPos = getMovePositionRec(nextPos, gateCoords, maxNMoves);
    if (bestPos.coords.size() == gateCoords.size()) {
      return bestPos;
    }
  }

  for (const auto& occCoord : occupiedNearbyCoords) {
    MultiQubitMovePos nextPos = MultiQubitMovePos(currentPos);
    nextPos.coords.emplace_back(occCoord);
    nextPos.nMoves += 2;
    auto bestPos = getMovePositionRec(nextPos, gateCoords, maxNMoves);
    if (bestPos.coords.size() == gateCoords.size()) {
      return bestPos;
    }
  }

  // if no position found, return empty
  return {};
}

MoveCombs NeutralAtomMapper::getAllMoveCombinations() {
  MoveCombs allMoves;
  for (const auto& op : this->frontLayerShuttling) {
    auto usedQubits    = op->getUsedQubits();
    auto usedHwQubits  = this->mapping.getHwQubits(usedQubits);
    auto usedCoordsSet = this->hardwareQubits.getCoordIndices(usedHwQubits);
    auto usedCoords =
        std::vector<CoordIndex>(usedCoordsSet.begin(), usedCoordsSet.end());
    auto bestPos = getBestMovePos(usedCoords);
    auto moves   = getMoveCombinationsToPosition(usedHwQubits, bestPos);
    allMoves.addMoveCombs(moves);
  }
  allMoves.removeLongerMoveCombs();
  return allMoves;
}

std::pair<MoveCombs, std::vector<std::pair<MoveComb, const qc::Operation*>> > NeutralAtomMapper::getAllMoveCombinationsWithOp() {
  MoveCombs allMoves;
  std::vector< std::pair< MoveComb, const qc::Operation*>> allMovesWithOp;
  int i=1;
  for (const auto& op : this->frontLayerShuttling) {
    auto usedQubits    = op->getUsedQubits();
    auto usedHwQubits  = this->mapping.getHwQubits(usedQubits);
    auto usedCoordsSet = this->hardwareQubits.getCoordIndices(usedHwQubits);
    auto usedCoords =
        std::vector<CoordIndex>(usedCoordsSet.begin(), usedCoordsSet.end());
    auto bestPos = getBestMovePos(usedCoords);
    auto moves   = getMoveCombinationsToPosition(usedHwQubits, bestPos);
    allMoves.addMoveCombs(moves);
    for(auto move : moves){
      allMovesWithOp.push_back( std::make_pair( move, op ) );
    }
  }
  allMoves.removeLongerMoveCombs();
  return make_pair(allMoves, allMovesWithOp);
}

CoordIndices NeutralAtomMapper::getBestMovePos(const CoordIndices& gateCoords) {
  size_t const maxMoves   = gateCoords.size() * 2;
  size_t const minMoves   = gateCoords.size();
  size_t       nMovesGate = maxMoves;
  // do a breadth first search for the best position
  // start with the used coords
  std::queue<CoordIndex> q;
  for (const auto& coord : gateCoords) {
    q.push(coord);
  }
  std::vector<CoordIndex> visited;

  auto finalBestPos = MultiQubitMovePos();
  while (!q.empty()) {
    auto coord = q.front();
    q.pop();
    if (std::find(visited.begin(), visited.end(), coord) != visited.end()) {
      continue;
    }
    visited.emplace_back(coord);
    MultiQubitMovePos currentPos;
    currentPos.coords.emplace_back(coord);
    if (this->hardwareQubits.isMapped(coord)) {
      if (std::find(gateCoords.begin(), gateCoords.end(), coord) !=
          gateCoords.end()) {
        currentPos.nMoves = 0;
      } else {
        currentPos.nMoves = 2;
      }
    } else {
      currentPos.nMoves = 1;
    }
    auto bestPos = getMovePositionRec(currentPos, gateCoords, nMovesGate);
    if (!bestPos.coords.empty() && bestPos.nMoves <= minMoves) {
      return bestPos.coords;
    }

    // min not yet reached, check nearby
    if (!bestPos.coords.empty()) {
      nMovesGate = std::min(nMovesGate, bestPos.nMoves);
    }
    for (const auto& nearbyCoord : this->arch.getNearbyCoordinates(coord)) {
      if (std::find(visited.begin(), visited.end(), nearbyCoord) ==
          visited.end()) {
        q.push(nearbyCoord);
      }
    }
  }
  throw qc::QFRException(
      "No move position found (check if enough free coords are available)");
}

MoveCombs
NeutralAtomMapper::getMoveCombinationsToPosition(HwQubits&     gateQubits,
                                                 CoordIndices& position) {
  if (position.empty()) {
    throw qc::QFRException("No position given");
  }
  // compute for each qubit the best position around it based on the cost of
  // the single move choose best one
  MoveCombs const      moveCombinations;
  std::set<CoordIndex> gateQubitCoords;
  for (const auto& gateQubit : gateQubits) {
    gateQubitCoords.emplace(this->hardwareQubits.getCoordIndex(gateQubit));
  }

  auto     remainingCoords = position;
  MoveComb moveComb;
  // compute cost for each candidate and each gateQubit
  auto remainingGateCoords = gateQubitCoords;
  // pre-filter away all gateQubitCoords which are already in the position
  for (auto it = remainingGateCoords.begin();
       it != remainingGateCoords.end();) {
    if (std::find(remainingCoords.begin(), remainingCoords.end(), *it) !=
        remainingCoords.end()) {
      remainingCoords.erase(
          std::find(remainingCoords.begin(), remainingCoords.end(), *it));
      it = remainingGateCoords.erase(it);
    } else {
      ++it;
    }
  }

  while (!remainingGateCoords.empty()) {
    auto currentGateQubit = *remainingGateCoords.begin();
    // compute costs and find best coord
    std::vector<std::pair<CoordIndex, qc::fp>> costs;
    for (const auto& remainingCoord : remainingCoords) {
      if (this->hardwareQubits.isMapped(remainingCoord)) {
        const auto moveAwayComb = getMoveAwayCombinations(
            currentGateQubit, remainingCoord, remainingCoords);
        for (const auto& moveAway : moveAwayComb) {
          auto cost = moveCostComb(moveAway);
          costs.emplace_back(remainingCoord, cost);
        }
      } else {
        auto cost = moveCost({currentGateQubit, remainingCoord});
        costs.emplace_back(remainingCoord, cost);
      }
    }
    // find minimal cost
    auto bestCost  = std::min_element(costs.begin(), costs.end(),
                                      [](const auto& cost1, const auto& cost2) {
                                       return cost1.second < cost2.second;
                                     });
    auto bestCoord = bestCost->first;
    if (this->hardwareQubits.isMapped(bestCoord)) {
      auto moveAwayComb =
          getMoveAwayCombinations(currentGateQubit, bestCoord, remainingCoords);
      for (const auto& moveAway : moveAwayComb) {
        moveComb.append(moveAway);
      }
    } else {
      moveComb.append(AtomMove{currentGateQubit, bestCoord});
    }
    remainingGateCoords.erase(currentGateQubit);
    remainingCoords.erase(
        std::find(remainingCoords.begin(), remainingCoords.end(), bestCoord));
  }
  return MoveCombs({moveComb});
}

MoveCombs
NeutralAtomMapper::getMoveAwayCombinations(CoordIndex          startCoord,
                                           CoordIndex          targetCoord,
                                           const CoordIndices& excludedCoords) {
  MoveCombs  moveCombinations;
  auto const originalVector    = this->arch.getVector(startCoord, targetCoord);
  auto const originalDirection = originalVector.direction;
  // Find move away target in the same direction as the original move
  auto moveAwayTargets = this->hardwareQubits.findClosestFreeCoord(
      targetCoord, originalDirection, excludedCoords);
  for (const auto& moveAwayTarget : moveAwayTargets) {
    const AtomMove move     = {startCoord, targetCoord};
    const AtomMove moveAway = {targetCoord, moveAwayTarget};
    moveCombinations.addMoveComb(MoveComb({moveAway, move}));
  }
  if (moveCombinations.empty()) {
    throw QMAPException("No move away target found");
  }
  return moveCombinations;
}

std::pair<uint32_t, qc::fp>
NeutralAtomMapper::estimateNumSwapGates(const qc::Operation* opPointer) {
  auto   usedQubits   = opPointer->getUsedQubits();
  auto   usedHwQubits = this->mapping.getHwQubits(usedQubits);
  qc::fp minNumSwaps  = 0;
  if (usedHwQubits.size() == 2) {
    SwapDistance minDistance = std::numeric_limits<SwapDistance>::max();
    for (const auto& hwQubit : usedHwQubits) {
      for (const auto& otherHwQubit : usedHwQubits) {
        if (hwQubit == otherHwQubit) {
          continue;
        }
        auto distance =
            this->hardwareQubits.getSwapDistance(hwQubit, otherHwQubit);
        minDistance = std::min(distance, minDistance);
      }
    }
    minNumSwaps = minDistance;
  } else { // multi-qubit gates
    auto bestPos = getBestMultiQubitPosition(opPointer);
    if (bestPos.empty()) {
      return {std::numeric_limits<SwapDistance>::max(),
              std::numeric_limits<qc::fp>::max()};
    }
    auto exactSwaps = getExactSwapsToPosition(opPointer, bestPos);
    if (exactSwaps.empty()) {
      return {std::numeric_limits<SwapDistance>::max(),
              std::numeric_limits<qc::fp>::max()};
    }
    for (const auto& [swap, weight] : exactSwaps) {
      auto [q1, q2] = swap;
      minNumSwaps += this->hardwareQubits.getSwapDistance(q1, q2, false);
    }
  }
  const qc::fp minTime = minNumSwaps * this->arch.getGateTime("swap");
  return {minNumSwaps, minTime};
}

std::pair<uint32_t, qc::fp>
NeutralAtomMapper::estimateNumMove(const qc::Operation* opPointer) {
  auto usedQubits   = opPointer->getUsedQubits();
  auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
  auto usedCoords   = this->hardwareQubits.getCoordIndices(usedHwQubits);
  // estimate the number of moves as:
  // compute distance between qubits
  // 1. for each free coord in the vicinity = 1 move with corresponding
  // distance
  // 2. for each occupied coord in the vicinity = 2 moves with corresponding
  // distance

  uint32_t minMoves = std::numeric_limits<uint32_t>::max();
  qc::fp   minTime  = std::numeric_limits<qc::fp>::max();
  for (const auto& coord : usedCoords) {
    qc::fp   totalTime  = 0;
    uint32_t totalMoves = 0;
    auto     nearbyFreeCoords =
        this->hardwareQubits.getNearbyFreeCoordinatesByCoord(coord);
    auto nearbyOccupiedCoords =
        this->hardwareQubits.getNearbyOccupiedCoordinatesByCoord(coord);
    auto otherQubitsIt = usedCoords.begin();
    auto nearbyFreeIt  = nearbyFreeCoords.begin();
    auto nearbyOccIt   = nearbyOccupiedCoords.begin();
    while (otherQubitsIt != usedCoords.end()) {
      auto otherCoord = *otherQubitsIt;
      if (otherCoord == coord) {
        otherQubitsIt++;
        continue;
      }
      if (nearbyFreeIt != nearbyFreeCoords.end()) {
        totalTime += this->arch.getVectorShuttlingTime(
            this->arch.getVector(otherCoord, *nearbyFreeIt));
        totalTime += this->arch.getShuttlingTime(qc::OpType::AodActivate) +
                     this->arch.getShuttlingTime(qc::OpType::AodDeactivate);
        nearbyFreeIt++;
        totalMoves++;
      } else if (nearbyOccIt != nearbyOccupiedCoords.end()) {
        totalTime += 2 * this->arch.getVectorShuttlingTime(
                             this->arch.getVector(otherCoord, *nearbyOccIt));
        totalTime +=
            2 * (this->arch.getShuttlingTime(qc::OpType::AodActivate) +
                 this->arch.getShuttlingTime(qc::OpType::AodDeactivate));
        nearbyOccIt++;
        totalMoves += 2;
      } else {
        throw std::runtime_error("No space to "
                                 "execute a multi-qubit gate. "
                                 "Check int radius. Op:" +
                                 opPointer->getName() + " nQubit: " +
                                 std::to_string(usedQubits.size()));
      }

      otherQubitsIt++;
    }

    if (totalTime < minTime) {
      minTime  = totalTime;
      minMoves = totalMoves;
    }
  }

  return {minMoves, minTime};
}
      

bool NeutralAtomMapper::swapGateBetter(const qc::Operation* opPointer) {
  auto [minNumSwaps, minTimeSwaps] = estimateNumSwapGates(opPointer);
  if (minNumSwaps == 0) {
    return true;
  }
  auto [minMoves, minTimeMoves] = estimateNumMove(opPointer);
  auto fidSwaps =
      std::exp(-minTimeSwaps * this->arch.getNqubits() /
               this->arch.getDecoherenceTime()) *
      std::pow(this->arch.getGateAverageFidelity("swap"), minNumSwaps);
  auto fidMoves =
      std::exp(-minTimeMoves * this->arch.getNqubits() /
               this->arch.getDecoherenceTime()) *
      std::pow(
          this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove) *
              this->arch.getShuttlingAverageFidelity(qc::OpType::AodActivate) *
              this->arch.getShuttlingAverageFidelity(qc::OpType::AodDeactivate),
          minMoves);

  return fidSwaps * parameters.gateWeight >
         fidMoves * parameters.shuttlingWeight;
}

std::vector<std::pair<const qc::Operation*, Bridge>> NeutralAtomMapper::findAllBridges(qc::QuantumComputation& qc, Swap bestSwap){
  std::vector<std::pair<const qc::Operation*, Bridge>> allBridges;
  for (const auto* op : this->frontLayerGate) {
    auto usedQubits = op->getUsedQubits();
    size_t usedQubitSize = usedQubits.size();
    if(usedQubitSize == 2){
      // logical qubits
      qc::Qubit q1 = *(usedQubits.begin());
      qc::Qubit q2 = *(std::next( usedQubits.begin(), 1));
      // hardware qubits
      HwQubit h1 = this->mapping.getHwQubit(q1);
      HwQubit h2 = this->mapping.getHwQubit(q2);
      qc::fp  dist = this->hardwareQubits.getSwapDistance(h1, h2);
      qc::fp  distswap;
      if(dist == 1){
        //only consider the case that bestSwap can reduce the swap distance of targetOp
        if (h1 == bestSwap.first) {
          distswap = this->hardwareQubits.getSwapDistance(bestSwap.second, h2);
        } else if (h1 == bestSwap.second) {
          distswap = this->hardwareQubits.getSwapDistance(bestSwap.first, h2);
        } else if (h2 == bestSwap.first) {
          distswap = this->hardwareQubits.getSwapDistance(h1, bestSwap.second);
        } else if (h2 == bestSwap.second) {
          distswap = this->hardwareQubits.getSwapDistance(h1, bestSwap.first);
        } 

        if(distswap==0){
          // get nearby
          HwQubits h1Near = this->hardwareQubits.getNearbyQubits(h1);
          HwQubits h2Near = this->hardwareQubits.getNearbyQubits(h2);
          Qubits Qbtw;
          for(const auto& h : h1Near){
            if(h2Near.find(h) != h2Near.end()){
              qc::Qubit qBtw = this->mapping.getCircQubit(h);
              if(qBtw < qc.getNqubits() && qBtw!=-1){
                Qbtw.insert( qBtw );
              }
            }
          }
          allBridges.emplace_back(op, Bridge(q1, q2, Qbtw));
        }
      }
    }
  }
  return allBridges;
}

std::vector<std::pair<const qc::Operation*, Bridge>> NeutralAtomMapper::compareCostBridgeSwap(
  std::vector<std::pair<const qc::Operation*, Bridge>> allBridges, Swap bestSwap,
  const qc::DAG& dag, NeutralAtomLayer& frontLayer, qc::QuantumComputation& qc){

  std::vector<std::pair<const qc::Operation*, Bridge>> ExecutableBridges;
  std::vector< std::tuple<qc::fp, qc::fp, const qc::Operation*, Bridge> > bridgePairCost;
  std::vector< std::tuple<qc::fp, qc::fp, const qc::Operation*, Bridge> > finalBridgePairCost;
  qc::Qubit qs1 = this->mapping.getCircQubit(bestSwap.first);
  qc::Qubit qs2 = this->mapping.getCircQubit(bestSwap.second);

  int i=0;
  for(auto bridgePair : allBridges){
    auto op = bridgePair.first;
    auto bridge = bridgePair.second;
    qc::Qubit q1 = std::get<0>(bridge);
    qc::Qubit q2 = std::get<1>(bridge);
    Qubits    Qb = std::get<2>(bridge);

    qc::fp requiredCNOTbridge = 0;
    qc::fp requiredCNOTswap   = 0;
    qc::fp bridgeCost = 0;
    qc::fp swapCost = 0;
    std::set<qc::Qubit> qGroup = {qs1, qs2};
    for(auto q : qGroup){
      if(q >= qc.getNqubits()) continue; 
      auto tempIter = dag[q].begin();
      qc::fp discountFactor = 1;
      bool isFirstBridge = true;
      while(tempIter < dag[q].end() && discountFactor > 0.5){
        auto* dagOp = (*tempIter)->get();
        auto dagOpUsedQubits = dagOp->getUsedQubits();
        if(dagOpUsedQubits.size() != 1){
          qc::fp distbefore = 0;
          qc::fp distswap   = 0;
          // if bridged gate -> ignore in calculating the cost
          if(isFirstBridge && dagOpUsedQubits.size()==2 && *dagOpUsedQubits.begin()==q1 && *(std::next(dagOpUsedQubits.begin(), 1))==q2){
            isFirstBridge = false;
          }
          else{
            for(auto it1 = dagOpUsedQubits.begin(); it1 != dagOpUsedQubits.end(); ++it1){
              qc::Qubit qi = *it1;
              HwQubit pi = this->mapping.getHwQubit(qi);
              for(auto it2 = std::next(it1); it2 != dagOpUsedQubits.end(); ++it2){
                qc::Qubit qj = *it2;
                HwQubit pj = this->mapping.getHwQubit(qj);

                qc::fp dist = this->hardwareQubits.getSwapDistance(pi, pj);
                distbefore += dist;
                if (pi == bestSwap.first) {
                  distswap += this->hardwareQubits.getSwapDistance(bestSwap.second, pj);
                } else if (pi == bestSwap.second) {
                  distswap += this->hardwareQubits.getSwapDistance(bestSwap.first, pj);
                } else if (pj == bestSwap.first) {
                  distswap += this->hardwareQubits.getSwapDistance(pi, bestSwap.second);
                } else if (pj == bestSwap.second) {
                  distswap += this->hardwareQubits.getSwapDistance(pi, bestSwap.first);
                } else {
                  distswap += dist;
                }
              }
            }
            requiredCNOTbridge += distbefore * discountFactor;
            requiredCNOTswap   += distswap   * discountFactor;
          }
        }
        else{
          discountFactor *= (1-this->parameters.decay);
        }
        tempIter++;
      }
    }
    // cost calculation
    qc::fp bridgeFidelOp = std::pow(this->arch.getGateAverageFidelity("cz"), 3*allBridges.size());
                              //std::pow(this->arch.getGateAverageFidelity("h"), 6*allBridges.size());
    qc::fp swapFidelOp   = std::pow(this->arch.getGateAverageFidelity("cz"), 3);
                              //std::pow(this->arch.getGateAverageFidelity("h"), 8);
    qc::fp bridgeFidelRequired = std::pow(this->arch.getGateAverageFidelity("cz"), requiredCNOTbridge);
    qc::fp swapFidelRequired   = std::pow(this->arch.getGateAverageFidelity("cz"), requiredCNOTswap);
    qc::fp bridgeTime = std::exp(-requiredCNOTbridge * this->arch.getNqubits() / this->arch.getDecoherenceTime());
    qc::fp swapTime   = std::exp(-requiredCNOTswap   * this->arch.getNqubits() / this->arch.getDecoherenceTime());

    bridgeCost = bridgeFidelOp * bridgeFidelRequired * bridgeTime;
    swapCost   = swapFidelOp   * swapFidelRequired   * swapTime;
    bridgePairCost.emplace_back( bridgeCost, swapCost, op, Bridge(q1, q2, Qb) );
    //}
  }

  // if other bridgePair shares the same control or same target -> remove larger cost 
  std::sort( bridgePairCost.begin(), bridgePairCost.end(),
            [](const auto& lhs, const auto& rhs){
              return std::get<0>(lhs) < std::get<0>(rhs);
            });
  std::vector<qc::Qubit> usedQubitsInBridge;
  for(auto bridgePair : bridgePairCost){
    auto bridge = std::get<3>(bridgePair);
    auto q1 = std::get<0>(bridge);
    auto q2 = std::get<1>(bridge);
    auto Qb = std::get<2>(bridge);

    auto itQ1 = std::find( usedQubitsInBridge.begin(), usedQubitsInBridge.end(), q1 );
    auto itQ2 = std::find( usedQubitsInBridge.begin(), usedQubitsInBridge.end(), q2 );
    for(auto qb : Qb){
      auto itQb = std::find( usedQubitsInBridge.begin(), usedQubitsInBridge.end(), qb );
      if( itQ1==usedQubitsInBridge.end() && itQ2==usedQubitsInBridge.end() && itQb==usedQubitsInBridge.end() ){
        auto bridgeCost = std::get<0>(bridgePair);
        auto swapCost   = std::get<1>(bridgePair);
        auto op     = std::get<2>(bridgePair);
        finalBridgePairCost.emplace_back( bridgeCost, swapCost, op, Bridge(q1, q2, {qb}) );
        usedQubitsInBridge.push_back( q1 );
        usedQubitsInBridge.push_back( q2 );
        usedQubitsInBridge.push_back( qb );
      }
    }
  }

  // calculate the swap cost and bridge cost
  qc::fp finalBridgeCost = 1;
  qc::fp finalSwapCost   = 1;
  for(auto bridgePair : finalBridgePairCost){
    auto bridgeCost = std::get<0>(bridgePair);
    auto swapCost   = std::get<1>(bridgePair);
    auto op     = std::get<2>(bridgePair);
    auto bridge = std::get<3>(bridgePair);
     
    auto q1 = std::get<0>(bridge);
    auto q2 = std::get<1>(bridge);
    auto Qb = std::get<2>(bridge);
    finalBridgeCost *= bridgeCost;
    finalSwapCost   *= swapCost;
  }

  // return ExecutableBridges
  finalSwapCost = std::pow(finalSwapCost, 1.0/bridgePairCost.size());
  if( finalBridgeCost >= finalSwapCost ){
    for(auto bridgePair : bridgePairCost){
      auto op     = std::get<2>(bridgePair);
      auto bridge = std::get<3>(bridgePair);
      ExecutableBridges.emplace_back(op, bridge);
    }
  }

  return ExecutableBridges;
}

void NeutralAtomMapper::updateMappingBridge(std::vector<std::pair<const qc::Operation*, Bridge>> ExecutableBridges,
    NeutralAtomLayer& frontLayer, NeutralAtomLayer& lookaheadLayer, qc::DAG& dag){
  // CX to Bridge
  //[q1] ---c--- = ---------c-----------c---
  //[qb]    |    = ---c---H-Z-H---c---H-Z-H-
  //[q2] -H-Z-H- = -H-Z-H-------H-Z-H-------

  //-> CZ to Bridge w/ QCO
  //[q1] -c- = -------c-----------c---
  //[qb]  |  = -c---H-Z-H---c---H-Z-H-
  //[q2] -Z- = -Z-----------Z--------

  GateList removeGates;
  for(auto bridgePair : ExecutableBridges){
    nBridges++;
    auto op = bridgePair.first;
    auto bridge = bridgePair.second;
    qc::Qubit  q1 = std::get<0>(bridge);
    qc::Qubit  q2 = std::get<1>(bridge);
    Qubits Qb = std::get<2>(bridge);
    qc::Qubit qb = *Qb.begin();
    if (this->parameters.verbose) {
      std::cout << "bridged " << q1 << " " << q2 << " by using " << qb;
      std::cout << "  physical qubits: ";
      std::cout << this->mapping.getHwQubit(q1);
      std::cout << " ";
      std::cout << this->mapping.getHwQubit(q2);
      std::cout << " ";
      std::cout << this->mapping.getHwQubit(qb);
      std::cout << '\n';
    }
    // add BR to mappedQc
    mappedQc.cz(qb, q2);
    mappedQc.h( qb);
    mappedQc.cz(q1, qb);
    mappedQc.h( qb);
    mappedQc.cz(qb, q2);
    mappedQc.h( qb);
    mappedQc.cz(q1, qb);
    mappedQc.h( qb);

    // remove original gate
    removeGates.push_back(op);
  }
  // remove original gate
  frontLayer.removeGatesAndUpdate(removeGates);
  lookaheadLayer.removeGatesAndUpdate(removeGates);
  for (const auto& gate : removeGates) {
    if (std::find(frontLayerGate.begin(), frontLayerGate.end(), gate) != frontLayerGate.end()) {
      frontLayerGate.erase(std::find(frontLayerGate.begin(), frontLayerGate.end(), gate));
    }
    if (std::find(lookaheadLayerGate.begin(), lookaheadLayerGate.end(), gate) != lookaheadLayerGate.end()) {
      lookaheadLayerGate.erase(std::find(lookaheadLayerGate.begin(), lookaheadLayerGate.end(), gate));
    }
    // update dag
    auto qubits = gate->getUsedQubits();
    for(auto q : qubits){
      auto it = std::find_if( dag[q].begin(), dag[q].end(),
                              [gate](const std::unique_ptr<qc::Operation>* op){
                                  return op->get() == gate;
                                  });
      if(it!=dag[q].end()){
        dag[q].erase(it);
      }
    }
  }
}

qc::QuantumComputation NeutralAtomMapper::findBestFlyingAncilla(qc::QuantumComputation& qc, const qc::Operation* targetOp){
  // information of operation
  qc::QuantumComputation bestAddedQc;
  uint32_t bestNumPassby = 0;
  auto usedQubits    = targetOp->getUsedQubits();
  auto QtargetSet = findQtargetSet( usedQubits );
  if(!QtargetSet.empty()){
    int idx = 0;
    int bestNumMoves = std::numeric_limits<int>::max();

    int bestIdx;
    std::set<qc::Qubit> bestQtarget;
    std::vector<qc::Qubit> bestQsource;
    // #F.A. => (Q_source, Q_target) iteration
    for(auto Qtarget : QtargetSet){
      int numFA = usedQubits.size() - Qtarget.size();
      uint32_t NumPassby = 0;
      std::set<qc::Qubit> Qsource;
      std::set_difference(
        usedQubits.begin(), usedQubits.end(),
        Qtarget.begin(), Qtarget.end(),
        std::inserter(Qsource, Qsource.end())
      );

      // permutate Qtarget & Qsource
      std::vector<qc::Qubit> QsourceVec( Qsource.begin(), Qsource.end() );
      do{
        qc::QuantumComputation addedQc;
        addedQc = qc::QuantumComputation(arch.getNpositions());
        qc::fp r_int = arch.getInteractionRadius();
        
        // hardware qubits & coord inices of Qtarget
        auto Htarget = this->mapping.getHwQubits( Qtarget );
        auto Ctarget = this->hardwareQubits.getCoordIndices( Htarget );

        std::vector<qc::Qubit> Qancilla;
        std::vector<HwQubit> Hsource, Hancilla;
        std::vector<CoordIndex> Csource, Cancilla;
        for(auto qs : QsourceVec){
          // hardware qubits & coord inices of Qsource
          auto hs = this->mapping.getHwQubit(qs);
          Hsource.push_back(hs);
          auto cs = this->hardwareQubits.getCoordIndex(hs);
          Csource.push_back(cs);
        }
        std::vector<CoordIndex> excludeCoords;
        std::copy(Ctarget.begin(), Ctarget.end(), std::back_inserter(excludeCoords));
        std::copy(Csource.begin(), Csource.end(), std::back_inserter(excludeCoords));

        std::vector<CoordIndex> passbyCtarget;
        std::copy(Ctarget.begin(), Ctarget.end(), std::back_inserter(passbyCtarget));

        std::vector<bool> needPassby(QsourceVec.size(), false);
        for(uint32_t i=0; i<QsourceVec.size(); i++){
          // find ancillaQubit of Qsource
          auto qi = QsourceVec[i];
          auto ci = Csource[i];
          auto cA = returnClosestAncillaCoord( ci, excludeCoords, qc); //excludeCoord: Ctarget, Csource, Cancilla
          auto hA = this->hardwareQubits.getHwQubit( cA );
          auto qA = this->mapping.getCircQubit( hA );
          Cancilla.push_back( cA );
          Hancilla.push_back( hA );
          Qancilla.push_back( qA );
          excludeCoords.push_back(cA);

          // 1. compare qs, qA 
          if(arch.getEuclideanDistance( this->arch.getCoordinate(ci), this->arch.getCoordinate(cA) ) > r_int ){
            // -> 1-1. passby (qA -> qi) 
            addedQc.passby( cA, {ci} );
            needPassby[i] = true;
            NumPassby++;
          }
          //-> cx (qi, qA)
          addedQc.cx(qi, qA);

          // 2. compare qA, {Qtarget, previous qAs}
          // -> 2-1. passby (cA -> Ct)
          if(needPassby[i]){
            addedQc.passby(ci, passbyCtarget);
          }
          else{
            addedQc.passby(cA, passbyCtarget);
          }
          NumPassby++;
          passbyCtarget.push_back(cA);
        }  
        // 2-2. mcz(Qtarget, Qancilla)
        qc::Controls mczControl;
        mczControl.insert( Qtarget.begin(), Qtarget.end() );
        mczControl.insert( Qancilla.begin(), Qancilla.end()-1 );
        qc::Qubit mczTarget = Qancilla.back();
        addedQc.mcz(mczControl, mczTarget);

        // 3. passby -> cx
        for(uint32_t i=0; i<QsourceVec.size(); i++){
          auto qA = Qancilla[i];
          auto qi = QsourceVec[i];
          auto cA = Cancilla[i];
          auto ci = Csource[i];
          if(needPassby[i]){
            addedQc.passby( passbyCtarget[0], {ci} );
            NumPassby++;
          }
          else{
            addedQc.move(passbyCtarget[0], cA);
          }
          addedQc.cx(QsourceVec[i], Qancilla[i]);

          // 4. move to original coordIdx
          if(needPassby[i]){
            addedQc.move(ci, cA);
          }
        }

        // find bestAddedQc
        qc::CircuitOptimizer::replaceMCXWithMCZ(addedQc);
        if(addedQc.size() < bestNumMoves){
          bestNumMoves = addedQc.size();
          bestAddedQc = addedQc;
          // for debugging
          bestIdx = idx;
          bestQtarget = Qtarget;
          bestQsource = QsourceVec;
          bestNumPassby = NumPassby;
        }
        //TODO: how to find the best addedQc?
        //else if(addedQc.size() == bestNumMoves){
        //}
      } while( std::next_permutation(QsourceVec.begin(), QsourceVec.end()) );

    }
    return bestAddedQc;
  }
}

std::set<std::set<qc::Qubit>> NeutralAtomMapper::findQtargetSet( std::set<qc::Qubit>& usedQubits ){
  std::set<std::set<qc::Qubit>> QtargetSet;
  auto numUsedQubits = usedQubits.size();
  SymmetricMatrix<qc::fp> gateQubitDistances(numUsedQubits);
  for(uint32_t i=0; i<numUsedQubits; ++i){
    for(uint32_t j=0; j<=i; ++j){
      if(i==j) gateQubitDistances(i, j) = 0;
      qc::Qubit qi = *(std::next( usedQubits.begin(), i));
      qc::Qubit qj = *(std::next( usedQubits.begin(), j));
      gateQubitDistances(i, j) = this->hardwareQubits.getSwapDistance( 
        this->mapping.getHwQubit(qi), this->mapping.getHwQubit(qj)); 
    }
  }

  size_t maxSize = 0;
  for(int i=0; i<numUsedQubits; ++i){
    std::vector<qc::Qubit> currentVec;
    qc::Qubit qi = *(std::next( usedQubits.begin(), i));
    currentVec.push_back( qi );
    for(int j=0; j<numUsedQubits; ++j){
      if(i!=j){
        qc::Qubit qj = *(std::next( usedQubits.begin(), j));
        bool isInteractable = true;
        for(auto& q : currentVec){
          auto it = std::find( usedQubits.begin(), usedQubits.end(), q );
          uint32_t idx;
          if(it!=usedQubits.end()){
            idx = std::distance( usedQubits.begin(), it );
          }
          if(gateQubitDistances(idx, j) != 0){
            isInteractable = false;
            break;
          }
        }
        if(isInteractable){
          currentVec.push_back( qj );
        }
      }
    }
    std::set<qc::Qubit> currentSet( currentVec.begin(), currentVec.end() );
    if(currentSet.size() > maxSize){
      maxSize = currentSet.size();
      QtargetSet.clear();
      QtargetSet.insert(currentSet);
    }
    else if(currentSet.size() == maxSize){
      QtargetSet.insert(currentSet);
    }
  }
  return QtargetSet;
}

CoordIndex NeutralAtomMapper::returnClosestAncillaCoord(const CoordIndex& c_target, const CoordIndices& excludeCoords, qc::QuantumComputation& qc){
  auto const originalVector    = this->arch.getVector(c_target+arch.getNcolumns(), c_target); //startCoord, targetCoord
  auto const originalDirection = originalVector.direction;
  auto AncillaTargets = this->hardwareQubits.findClosestAncillaCoord(c_target, originalDirection, qc.getNqubits(), excludeCoords); 
  return AncillaTargets[0];
}

std::pair<bool, bool> NeutralAtomMapper::compareCostMoveFAandFQ(MoveComb bestMoveComb, qc::QuantumComputation& bestFA, qc::QuantumComputation& bestFQ,
  const qc::Operation* targetOp, const qc::DAG& dag, NeutralAtomLayer& frontLayer, qc::QuantumComputation& qc){
  
  //make circuit using the bestMoveComb
  qc::QuantumComputation bestMove;
  bestMove = qc::QuantumComputation(arch.getNpositions());
  std::vector<HwQubit> usedCoordsInMoveComb1, usedCoordsInMoveComb2;
  // layout save
  qc::Permutation layoutBefore = this->hardwareQubits.getHwToCoordIdx();
  qc::Permutation layoutAfter  = layoutBefore;
  std::vector<qc::Qubit> Qmoved;
  for(auto move : bestMoveComb.moves){
    auto c1 = move.first;
    auto c2 = move.second;
    auto itC1 = std::find( usedCoordsInMoveComb1.begin(), usedCoordsInMoveComb1.end(), c1 );
    auto itC2 = std::find( usedCoordsInMoveComb2.begin(), usedCoordsInMoveComb2.end(), c2 );
    if(itC1==usedCoordsInMoveComb1.end() && itC2==usedCoordsInMoveComb2.end()){
      bestMove.move(c1, c2);
      usedCoordsInMoveComb1.push_back( c1 );
      usedCoordsInMoveComb2.push_back( c2 );
      auto h1 = this->hardwareQubits.getHwQubit(c1);
      auto h2 = this->hardwareQubits.getHwQubit(c2);
      auto q1 = this->mapping.getCircQubit(h1);
      //layout update
      layoutAfter[h1] = c2;
      Qmoved.push_back(q1);
    }
  }
  bestMove.mcz( targetOp->getControls(), *targetOp->getTargets().begin() );
  
  // (1) calculate the Fidelity_requiredOp
  qc::fp beforeFidelOp = 1;
  qc::fp afterFidelOp  = 1;
  // check the fidelity of required op (SWAP? or MOVE?)
  int i=0;
  for(auto q : Qmoved){
    auto h = this->mapping.getHwQubit(q);
    qc::fp discountFactor = 1;
    if(q>=qc.getNqubits()) continue;
    auto tempIter = dag[q].begin();
    bool isTargetOp = true;
    while(tempIter < dag[q].end() && discountFactor > 0.5){
      auto* dagOp = (*tempIter)->get();
      auto dagOpUsedQubits = dagOp->getUsedQubits();
      if(dagOpUsedQubits.size() != 1){
        qc::fp distbefore = 0;
        qc::fp distafter  = 0;
        qc::fp numMoveRequired = 0;
        // if target op -> ignore in calculating the cost
        if(isTargetOp && dagOp == targetOp){
          isTargetOp = false;
        }
        else{
          for(auto it1 = dagOpUsedQubits.begin(); it1 != dagOpUsedQubits.end(); ++it1){
            qc::Qubit qi = *it1;
            HwQubit hi = this->mapping.getHwQubit(qi);
            for(auto it2 = std::next(it1); it2 != dagOpUsedQubits.end(); ++it2){
              qc::Qubit qj = *it2;
              HwQubit hj = this->mapping.getHwQubit(qj);

              qc::fp dist = this->hardwareQubits.getSwapDistance(hi, hj);
              distbefore += dist;

              auto itQi = std::find(Qmoved.begin(), Qmoved.end(), qi);
              auto itQj = std::find(Qmoved.begin(), Qmoved.end(), qj);
              if(itQi == Qmoved.end()){
                if(itQj == Qmoved.end()){
                  //qi, qj is not moved
                  distafter += dist;
                }
                else{
                  //qj is moved to layoutAfter[hj]
                  auto hMoved = this->hardwareQubits.getHwQubit(layoutAfter[hj]);
                  if( hMoved < this->arch.getNqubits()){
                    distafter += this->hardwareQubits.getSwapDistance( hMoved, hi );
                  }
                  else{
                    // shuttling is required
                    numMoveRequired++;
                  }
                }
              }
              else{
                if(itQj == Qmoved.end()){
                  //qi is moved to layoutAfter[hi]
                  auto hMoved = this->hardwareQubits.getHwQubit(layoutAfter[hi]);
                  if( hMoved < this->arch.getNqubits()){
                    distafter += this->hardwareQubits.getSwapDistance( hMoved, hj );
                  }
                  else{
                    // shuttling is required
                    numMoveRequired++;
                  }
                }
                else{
                  //qi, qj is moved to layoutAfter[hi] and layoutAfter[hj] each.
                  auto hiMoved = this->hardwareQubits.getHwQubit(layoutAfter[hi]);
                  auto hjMoved = this->hardwareQubits.getHwQubit(layoutAfter[hj]);
                  if( hiMoved < this->arch.getNqubits() ){
                    if( hjMoved < this->arch.getNqubits() ){
                      //hi and hj are both in the coordinate
                      distafter += this->hardwareQubits.getSwapDistance( hiMoved, hjMoved );
                    }
                    else{
                      //hj is far away
                      numMoveRequired++;
                    }
                  }
                  else{
                    if( hjMoved < this->arch.getNqubits() ){
                      //hi is far away
                      numMoveRequired++;
                    }
                    else{
                      //hi & hj are far away
                      distafter += this->hardwareQubits.getSwapDistance( hiMoved, hjMoved );
                    }
                  }
                }
              }
            }
          }
          beforeFidelOp *= std::pow( this->arch.getGateAverageFidelity("swap"), distbefore*discountFactor )
                         * std::exp(-distbefore*discountFactor*this->arch.getNqubits() / this->arch.getDecoherenceTime());
          qc::fp afterFidelOpSwap = 1;
          qc::fp afterFidelOpMove = 1;
          if(distafter!=0){
            afterFidelOpSwap *= std::pow( this->arch.getGateAverageFidelity("swap"), distafter*discountFactor )
                             * std::exp(-distafter*discountFactor*this->arch.getNqubits() / this->arch.getDecoherenceTime());
          }
          if(numMoveRequired!=0){
            afterFidelOpMove *= std::pow( this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove) *
                                          this->arch.getShuttlingAverageFidelity(qc::OpType::AodActivate) *
                                          this->arch.getShuttlingAverageFidelity(qc::OpType::AodDeactivate),
                                          numMoveRequired*discountFactor)
                                    * std::exp(-numMoveRequired*discountFactor * this->arch.getNqubits() / this->arch.getDecoherenceTime());
          }
          afterFidelOp *= afterFidelOpSwap * afterFidelOpMove;
        }
      }
      else{
        discountFactor *= (1-this->parameters.decay);
      }
      tempIter++;
    }
  }
 
  // (2) calculate the Fidelity_op @addedQc
  qc::fp moveFidelOp = 1;
  qc::fp FaFidelOp = 1;
  qc::fp FqFidelOp = 1;
  //-----------------------------------------------------
  for(const auto& opPtr : bestMove){
    const auto* op = opPtr.get();
    if(op->getType() == qc::OpType::Move){
      moveFidelOp *= this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove) *
                     this->arch.getShuttlingAverageFidelity(qc::OpType::AodActivate) *
                     this->arch.getShuttlingAverageFidelity(qc::OpType::AodDeactivate)
                  * std::exp(-0 * this->arch.getNqubits() / this->arch.getDecoherenceTime());
    }
    else{
      auto usedQubits = op->getUsedQubits();
      if(op->getType() == qc::OpType::Z){
        std::string gateName;
        if(usedQubits.size()==2) gateName = "cz";
        else if(usedQubits.size()==3) gateName = "ccz";
        else if(usedQubits.size()==4) gateName = "cccz";
        else gateName = "ccccz";
        moveFidelOp *= this->arch.getGateAverageFidelity(gateName)
                    * std::exp(-0*this->arch.getNqubits() / this->arch.getDecoherenceTime());
      }
      else{
        moveFidelOp *= this->arch.getGateAverageFidelity(toString(op->getType()))
                    * std::exp(-0*this->arch.getNqubits() / this->arch.getDecoherenceTime());
      }
    }
  }
  //-----------------------------------------------------
  for(const auto& opPtr : bestFA){
    const auto* op = opPtr.get();
    if(op->getType() == qc::OpType::PassBy){
      FaFidelOp *= this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove) *
                   this->arch.getShuttlingAverageFidelity(qc::OpType::AodActivate) *
                   this->arch.getShuttlingAverageFidelity(qc::OpType::AodDeactivate)
                * std::exp(-0 * this->arch.getNqubits() / this->arch.getDecoherenceTime());
    }  
    else if(op->getType() == qc::OpType::Move){
      FaFidelOp *= this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove)
                * std::exp(-0 * this->arch.getNqubits() / this->arch.getDecoherenceTime());
    }
    else{
      auto usedQubits = op->getUsedQubits();
      if(op->getType() == qc::OpType::Z){
        std::string gateName;
        if(usedQubits.size()==2) gateName = "cz";
        else if(usedQubits.size()==3) gateName = "ccz";
        else if(usedQubits.size()==4) gateName = "cccz";
        else gateName = "ccccz";
        FaFidelOp *= this->arch.getGateAverageFidelity(gateName)
                  * std::exp(-0*this->arch.getNqubits() / this->arch.getDecoherenceTime());
      }
      else{
        FaFidelOp *= this->arch.getGateAverageFidelity(toString(op->getType()))
                  * std::exp(-0*this->arch.getNqubits() / this->arch.getDecoherenceTime());
      }
    }
  }
  //-----------------------------------------------------
  for(const auto& opPtr : bestFQ){
    const auto* op = opPtr.get();
    if(op->getType() == qc::OpType::PassBy){
      FqFidelOp *= this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove) *
                   this->arch.getShuttlingAverageFidelity(qc::OpType::AodActivate) *
                   this->arch.getShuttlingAverageFidelity(qc::OpType::AodDeactivate)
                * std::exp(-0 * this->arch.getNqubits() / this->arch.getDecoherenceTime());
    }  
    else if(op->getType() == qc::OpType::Move){
      FqFidelOp *= this->arch.getShuttlingAverageFidelity(qc::OpType::AodMove)
                * std::exp(-0 * this->arch.getNqubits() / this->arch.getDecoherenceTime());
    }
    else{
      auto usedQubits = op->getUsedQubits();
      if(op->getType() == qc::OpType::Z){
        std::string gateName;
        if(usedQubits.size()==2) gateName = "cz";
        else if(usedQubits.size()==3) gateName = "ccz";
        else if(usedQubits.size()==4) gateName = "cccz";
        else gateName = "ccccz";
        FqFidelOp *= this->arch.getGateAverageFidelity(gateName)
                  * std::exp(-0*this->arch.getNqubits() / this->arch.getDecoherenceTime());
      }
      else{
        FqFidelOp *= this->arch.getGateAverageFidelity(toString(op->getType()))
                  * std::exp(-0*this->arch.getNqubits() / this->arch.getDecoherenceTime());
      }
    }
  }

  // compare final fidelity -> update Mapping
  moveFidelOp *= afterFidelOp;
  FaFidelOp *= beforeFidelOp;
  FqFidelOp *= beforeFidelOp;

  if(FaFidelOp > moveFidelOp || FqFidelOp > moveFidelOp){
    if(FaFidelOp > FqFidelOp){
      nFAncillas++;
      return std::make_pair(false, true); //useShuttling, useFlyingAncilla;
    }
    else{
      nFQubits++;
      return std::make_pair(false, false); //useShuttling, useFlyingAncilla;
    }
  }
  else{
    return std::make_pair(true, false); //useShuttling, useFlyingAncilla;
  }
}

void NeutralAtomMapper::updateMappingFlyingAncilla(qc::QuantumComputation& bestFA, const qc::Operation* targetOp, 
  NeutralAtomLayer& frontLayer, NeutralAtomLayer& lookaheadLayer){
  // add bestFA to mappedQc
  // TODO: solve the error (gate type pass_by could not be converted to OpenQASM)

  for(const auto& opPtr : bestFA){
    const auto* op = opPtr.get();
    auto usedQubits = op->getUsedQubits();
    if(op->getType() == qc::OpType::H){
      mappedQc.h(*usedQubits.begin());
    }
    if(op->getType() == qc::OpType::Z){
      if(usedQubits.size() > 1){
        mappedQc.mcz(op->getControls(), op->getTargets()[0]);
      }
      else{
        mappedQc.z(*usedQubits.begin());
      }
    }
    if(op->getType() == qc::OpType::PassBy){
      mappedQc.passby(*op->getControls().begin(), op->getTargets());
    }
    if(op->getType() == qc::OpType::Move){
      mappedQc.move(op->getTargets()[0], op->getTargets()[1]);
    }
  }

  // remove original gate
  GateList removeGates;
  removeGates.push_back(targetOp);
  frontLayer.removeGatesAndUpdate(removeGates);
  lookaheadLayer.removeGatesAndUpdate(removeGates);
  if(std::find( frontLayerShuttling.begin(), frontLayerShuttling.end(), targetOp ) != frontLayerShuttling.end()){
    frontLayerShuttling.erase( std::find( frontLayerShuttling.begin(), frontLayerShuttling.end(), targetOp ) );
  }
  if(std::find( lookaheadLayerShuttling.begin(), lookaheadLayerShuttling.end(), targetOp ) != lookaheadLayerShuttling.end()){
    lookaheadLayerShuttling.erase( std::find( lookaheadLayerShuttling.begin(), lookaheadLayerShuttling.end(), targetOp ) );
  }
}

void NeutralAtomMapper::updateMappingFlyingQubit(qc::QuantumComputation& bestFQ, const qc::Operation* targetOp, 
  NeutralAtomLayer& frontLayer, NeutralAtomLayer& lookaheadLayer){
  // add bestFQ to mappedQc
  // TODO: solve the error (gate type pass_by could not be converted to OpenQASM)

  for(const auto& opPtr : bestFQ){
    const auto* op = opPtr.get();
    auto usedQubits = op->getUsedQubits();
    if(op->getType() == qc::OpType::H){
      mappedQc.h(*usedQubits.begin());
    }
    if(op->getType() == qc::OpType::Z){
      if(usedQubits.size() > 1){
        mappedQc.mcz(op->getControls(), op->getTargets()[0]);
      }
      else{
        mappedQc.z(*usedQubits.begin());
      }
    }
    if(op->getType() == qc::OpType::PassBy){
      mappedQc.passby(*op->getControls().begin(), op->getTargets());
    }
    if(op->getType() == qc::OpType::Move){
      mappedQc.move(op->getTargets()[0], op->getTargets()[1]);
    }
  }

  // remove original gate
  GateList removeGates;
  removeGates.push_back(targetOp);
  frontLayer.removeGatesAndUpdate(removeGates);
  lookaheadLayer.removeGatesAndUpdate(removeGates);
  if(std::find( frontLayerShuttling.begin(), frontLayerShuttling.end(), targetOp ) != frontLayerShuttling.end()){
    frontLayerShuttling.erase( std::find( frontLayerShuttling.begin(), frontLayerShuttling.end(), targetOp ) );
  }
  if(std::find( lookaheadLayerShuttling.begin(), lookaheadLayerShuttling.end(), targetOp ) != lookaheadLayerShuttling.end()){
    lookaheadLayerShuttling.erase( std::find( lookaheadLayerShuttling.begin(), lookaheadLayerShuttling.end(), targetOp ) );
  }
}

qc::QuantumComputation NeutralAtomMapper::findBestFlyingQubit(qc::QuantumComputation& qc, const qc::Operation* targetOp){
  qc::QuantumComputation bestAddedQc;
  uint32_t bestNumPassby = std::numeric_limits<int>::max();
  int32_t bestNumParallel = std::numeric_limits<int>::min();
 
  auto usedQubits    = targetOp->getUsedQubits();
  auto QtargetSet = findQtargetSet( usedQubits );
  if(!QtargetSet.empty()){
    int idx = 0;
    int bestNumMoves = std::numeric_limits<int>::max();

    int bestIdx;
    std::set<qc::Qubit> bestQtarget;
    std::vector<qc::Qubit> bestQsource;
    // #F.A. => (Q_source, Q_target) iteration
    for(auto Qtarget : QtargetSet){
      int numFQ = usedQubits.size() - Qtarget.size();
      uint32_t NumPassby = 0;
      std::set<qc::Qubit> Qsource;
      std::set_difference(
        usedQubits.begin(), usedQubits.end(),
        Qtarget.begin(), Qtarget.end(),
        std::inserter(Qsource, Qsource.end())
      );

      // permutate Qtarget & Qsource
      std::vector<qc::Qubit> QsourceVec( Qsource.begin(), Qsource.end() );
      do{
        qc::QuantumComputation addedQc;
        addedQc = qc::QuantumComputation(arch.getNpositions());
        int32_t NumParallel = 0;
        // generate MoveVector Group
        std::vector<MoveVector> MoveVectorGroup;

        // hardware qubits & coord inices
        auto Htarget = this->mapping.getHwQubits( Qtarget );
        auto Ctarget = this->hardwareQubits.getCoordIndices( Htarget );
        std::vector<HwQubit> Hsource;
        std::vector<CoordIndex> Csource;
        for(auto qs : QsourceVec){
          auto hs = this->mapping.getHwQubit(qs);
          Hsource.push_back(hs);
          auto cs = this->hardwareQubits.getCoordIndex(hs);
          Csource.push_back(cs);
        }
        // 1. passby
        std::vector<qc::Qubit> passbyCtarget;
        std::copy(Ctarget.begin(), Ctarget.end(), std::back_inserter(passbyCtarget));
        for(uint32_t i=0; i<QsourceVec.size(); i++){
          auto qi = QsourceVec[i];
          auto ci = Csource[i];

          addedQc.passby(ci, passbyCtarget);
          auto const Vector_i = this->arch.getVector(ci, *(Ctarget.begin()) );
          MoveVectorGroup.push_back(Vector_i);
          NumPassby++;
          passbyCtarget.push_back(ci);
        }  
        // 2. mcz
        addedQc.mcz(targetOp->getControls(), *targetOp->getTargets().begin()); 

        // 3 move to original coordIdx
        for(uint32_t i=0; i<QsourceVec.size(); i++){
          auto qi = QsourceVec[i];
          auto ci = Csource[i];
          addedQc.move(*(Ctarget.begin()), ci);
        }  

        // find bestAddedQc based on the paralellism
        qc::CircuitOptimizer::replaceMCXWithMCZ(addedQc);
        if(MoveVectorGroup.size() >=2){
          for(int i=1; i<MoveVectorGroup.size()-1 ; i++){
            bool notOverlap = !MoveVectorGroup[0].overlap( MoveVectorGroup[i] );
            bool sameDirection = (MoveVectorGroup[0].direction == MoveVectorGroup[i].direction);
            bool include = (MoveVectorGroup[0].include(MoveVectorGroup[i]) || MoveVectorGroup[i].include(MoveVectorGroup[0]));
            bool parallel = notOverlap || sameDirection || include;
            if(parallel) NumParallel++;
          }
          if(bestNumParallel < NumParallel){
            bestNumPassby = NumPassby;
            bestNumParallel = NumParallel;
            bestAddedQc = addedQc;
          }
        }
        else{
          if(bestNumPassby > NumPassby){
            bestNumPassby = NumPassby;
            bestNumParallel = NumParallel;
            bestAddedQc = addedQc;
          }
        }

      } while( std::next_permutation(QsourceVec.begin(), QsourceVec.end()) );

    }
    // return the best result
    return bestAddedQc;
  }
}

} // namespace na
