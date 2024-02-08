#include "LogicBlock.hpp"

#include "Model.hpp"

namespace logicbase {

void LogicBlock::assertFormula(const LogicTerm& a) {
  if (a.getOpType() == OpType::AND) {
    for (const auto& clause : a.getNodes()) {
      clauses.insert(clause);
    }
  } else {
    clauses.insert(a);
  }
}

LogicTerm LogicBlock::makeVariable(const std::string& name, CType type,
                                   uint16_t bvSize) {
  if (type == CType::BITVECTOR && bvSize == 0) {
    throw std::invalid_argument("bv_size must be > 0");
  }
  return LogicTerm(name, type, this, bvSize);
}

void LogicBlock::reset() {
  delete model;
  model = nullptr;
  clauses.clear();
  internalReset();
  gid = 0U;
}

void LogicBlock::dump(const LogicTerm& a, std::ostream& stream) {
  a.prettyPrint(stream);
}

void LogicBlockOptimizer::reset() {
  model = nullptr;
  clauses.clear();
  weightedTerms.clear();
  internalReset();
  gid = 0U;
}

void LogicBlockOptimizer::weightedTerm(const LogicTerm& a, double weight) {
  weightedTerms.emplace_back(a, weight);
}
} // namespace logicbase
