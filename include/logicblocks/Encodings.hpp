#pragma once

#include "LogicBlock.hpp"

#include <cmath>
#include <set>
#include <utility>
#include <vector>

namespace encodings {

using namespace logicbase;

struct NestedVar {
  explicit NestedVar(LogicTerm v) : var(std::move(v)), list(){};
  NestedVar(LogicTerm v, std::vector<NestedVar> l)
      : var(std::move(v)), list(std::move(l)) {}
  LogicTerm              var = LogicTerm::noneTerm();
  std::vector<NestedVar> list;
};

struct WeightedVar {
  WeightedVar(LogicTerm v, const int w) : var(std::move(v)), weight(w) {}
  LogicTerm var    = LogicTerm::noneTerm();
  int       weight = 0;
};
inline bool operator<(const WeightedVar& rhs, const WeightedVar& lhs) {
  return rhs.weight < lhs.weight;
}
inline bool operator==(const WeightedVar& rhs, const WeightedVar& lhs) {
  return rhs.weight == lhs.weight && rhs.var.getID() == lhs.var.getID();
}

enum class Type { Uninitialized, AuxVar, ProgramVar };
struct SavedLit {
  SavedLit() : var(LogicTerm::noneTerm()) {}
  SavedLit(Type t, LogicTerm v) : type(t), var(std::move(v)) {}
  Type      type = Type::Uninitialized;
  LogicTerm var  = LogicTerm::noneTerm();
};

LogicTerm atMostOneCmdr(const std::vector<NestedVar>& subords,
                        const LogicTerm& cmdrVar, LogicBlock* logic);

LogicTerm exactlyOneCmdr(const std::vector<NestedVar>& subords,
                         const LogicTerm& cmdrVar, LogicBlock* logic);

LogicTerm naiveExactlyOne(const std::vector<LogicTerm>& clauseVars);

LogicTerm naiveAtMostOne(const std::vector<LogicTerm>& clauseVars);

LogicTerm naiveAtLeastOne(const std::vector<LogicTerm>& clauseVars);

LogicTerm atMostOneBiMander(const std::vector<LogicTerm>& vars,
                            LogicBlock*                   logic);

std::vector<NestedVar> groupVars(const std::vector<LogicTerm>& vars,
                                 std::size_t                   maxSize);
std::vector<NestedVar> groupVarsAux(const std::vector<NestedVar>& vars,
                                    std::size_t                   maxSize);

std::vector<std::vector<LogicTerm>>
groupVarsBimander(const std::vector<LogicTerm>& vars, std::size_t groupCount);

[[maybe_unused]] LogicTerm BuildBDD(const std::set<WeightedVar>&  inputLiterals,
                                    const std::vector<LogicTerm>& vars, int leq,
                                    LogicBlock* lb);
LogicTerm BuildBDD(uint64_t index, int64_t curSum, int64_t maxSum, int64_t k,
                   const std::vector<WeightedVar>& inputLiterals,
                   const std::vector<LogicTerm>& vars, LogicTerm& formula,
                   LogicTerm& trueLit, LogicBlock* lb);
} // namespace encodings
