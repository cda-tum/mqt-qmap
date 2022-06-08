/*
* This file is part of the MQT QMAP library which is released under the MIT license.
* See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
*/

#include "LBEncodings.hpp"

std::map<std::pair<unsigned long, long>, SavedLit> history;

LogicTerm NaiveExactlyOne(const std::vector<LogicTerm>& clauseVars) {
    return NaiveAtLeastOne(clauseVars) && NaiveAtMostOne(clauseVars);
}

LogicTerm NaiveAtLeastOne(const std::vector<LogicTerm>& clauseVars) {
    LogicTerm naiveAtLeastOne = LogicTerm(false);
    for (const LogicTerm& x: clauseVars) {
        naiveAtLeastOne = naiveAtLeastOne || x;
    }
    return naiveAtLeastOne;
}

LogicTerm NaiveAtMostOne(const std::vector<LogicTerm>& clauseVars) {
    LogicTerm naiveAtMostOne = LogicTerm(true);
    for (std::size_t i = 0U; i < clauseVars.size() - 1; i++) {
        for (std::size_t j = i + 1U; j < clauseVars.size(); j++) {
            naiveAtMostOne = naiveAtMostOne && (!clauseVars[i] || !clauseVars[j]);
        }
    }
    return naiveAtMostOne;
}

LogicTerm AtMostOneBiMander(const std::vector<LogicTerm>& vars) {
    auto                   subords = groupVarsBimander(vars, vars.size() / 2);
    LogicTerm              ret     = LogicTerm(true);
    std::vector<LogicTerm> binary_vars{};
    auto                   m = subords.size();
    for (int j = 0; j < std::ceil(std::log2(m)); j++) {
        binary_vars.emplace_back();
    }
    for (std::size_t i = 0U; i < m; i++) {
        LogicTerm binary = LogicTerm(true);
        for (std::size_t h = 0U; h < subords[i].size(); h++) {
            LogicTerm b2 = LogicTerm(true);
            for (std::size_t j = 0U; j < static_cast<std::size_t>(std::ceil(std::log2(m))); j++) {
                if (i & 1U << j) {
                    b2 = b2 && (!subords[i][h] || binary_vars[static_cast<int>(j)]);
                } else {
                    b2 = b2 && (!subords[i][h] || !binary_vars[static_cast<int>(j)]);
                }
            }
            binary = binary && b2;
        }
        ret = ret && binary && NaiveAtMostOne(vars);
    }
    return ret;
}

LogicTerm ExactlyOneCMDR(const std::vector<NestedVar>& subords, const LogicTerm& cmdrVar, LogicBlock* logic) {
    auto                   ret = LogicTerm(true);
    std::vector<LogicTerm> clauseVars{};
    clauseVars.reserve(subords.size());
    for (const auto& it: subords) {
        if (it.var.getOpType() != OpType::None) {
            clauseVars.emplace_back(it.var);
        } else {
            LogicTerm localCdr = logic->makeVariable("cdr_var");
            clauseVars.emplace_back(localCdr);
            ret = ret && ExactlyOneCMDR(it.list, localCdr, logic);
        }
    }
    if (cmdrVar.getOpType() == OpType::Variable) {
        clauseVars.emplace_back(!cmdrVar);
    }
    return ret && NaiveExactlyOne(clauseVars);
}

LogicTerm AtMostOneCMDR(const std::vector<NestedVar>& subords, const LogicTerm& cmdrVar, LogicBlock* logic) {
    auto                   ret = LogicTerm(true);
    std::vector<LogicTerm> clauseVars;
    clauseVars.reserve(subords.size());
    for (auto& it: subords) {
        if (it.var.getOpType() != OpType::None) {
            clauseVars.emplace_back(it.var);
        } else {
            LogicTerm localCdr = logic->makeVariable("cdr_var");
            clauseVars.emplace_back(localCdr);
            ret = ret && AtMostOneCMDR(it.list, localCdr, logic);
        }
    }
    if (cmdrVar.getOpType() == OpType::Variable) {
        clauseVars.emplace_back(!cmdrVar);
    }
    return ret && NaiveAtMostOne(clauseVars);
}

std::vector<NestedVar> groupVars(const std::vector<LogicTerm>& vars, std::size_t maxSize) {
    std::vector<NestedVar> vVars;
    vVars.reserve(vars.size());
    for (const auto& var: vars) {
        vVars.emplace_back(var);
    }
    if (vVars.size() <= 6U) {
        return vVars;
    }
    return groupVarsAux(vVars, maxSize);
}

std::vector<NestedVar> groupVarsAux(const std::vector<NestedVar>& vars, std::size_t maxSize) {
    auto numVars = vars.size();
    if (numVars <= maxSize) {
        return vars;
    }
    std::vector<NestedVar> ret{};
    size_t                 numGr = numVars / maxSize;
    ret.reserve(numGr);
    for (unsigned int i = 0U; i < numGr; i++) {
        auto from = vars.begin();
        std::advance(from, static_cast<long>(i * numVars / numGr));
        auto to = from;
        if (i == numGr - 1U) {
            to = vars.end();
        } else {
            std::advance(to, static_cast<long>(numVars / numGr));
        }
        ret.emplace_back(LogicTerm::noneTerm(), std::vector<NestedVar>(from, to));
    }
    return groupVarsAux(ret, maxSize);
}

std::vector<std::vector<LogicTerm>>
groupVarsBimander(const std::vector<LogicTerm>& vars, std::size_t groupCount) {
    std::vector<std::vector<LogicTerm>> result{};
    auto                                chunkSize = vars.size() / groupCount;

    for (size_t i = 0U; i < vars.size(); i += chunkSize) {
        auto from = vars.begin();
        std::advance(from, static_cast<long>(i));
        auto to  = vars.begin();
        auto end = std::min(vars.size(), i + chunkSize);
        std::advance(to, static_cast<long>(end));
        result.emplace_back(from, to);
    }

    return result;
}

LogicTerm BuildBDD(const std::set<WeightedVar>&  inputLiterals,
                   const std::vector<LogicTerm>& vars, int leq) {
    std::vector<WeightedVar> literals(inputLiterals.begin(), inputLiterals.end());
    history.clear();
    long k      = leq;
    long maxSum = 0;
    for (const auto& l: literals) {
        maxSum += l.weight;
    }
    LogicTerm true_lit{};
    LogicTerm formula{};
    LogicTerm result = BuildBDD(0U, 0, maxSum, k, literals, vars, formula, true_lit);
    return result && formula;
}

LogicTerm BuildBDD(unsigned long index, long curSum, long maxSum, long k,
                   const std::vector<WeightedVar>& inputLiterals,
                   const std::vector<LogicTerm>& vars, LogicTerm& formula,
                   LogicTerm& true_lit) {
    if (curSum + maxSum < k) {
        return true_lit;
    }
    if (curSum >= k) {
        return !(true_lit);
    }
    if (history.count({inputLiterals[index].var.getID(), curSum}) > 0) {
        const SavedLit& l = history[{inputLiterals[index].var.getID(), curSum}];
        if (l.type == Type::ProgramVar) {
            return !(l.var);
        } else {
            return l.var;
        }
    }

    LogicTerm high = BuildBDD(index + 1U, curSum + inputLiterals[index].weight,
                              maxSum - inputLiterals[index].weight, k,
                              inputLiterals, vars, formula, true_lit);
    LogicTerm low  = BuildBDD(index + 1U, curSum, maxSum - inputLiterals[index].weight, k,
                              inputLiterals, vars, formula, true_lit);

    if (high.deepEquals(low)) {
        return high;
    }

    LogicTerm node = LogicTerm(true);

    if (high.deepEquals(!true_lit) && low.deepEquals(true_lit)) {
        node                                                              = !(inputLiterals[index].var);
        history[std::make_pair(inputLiterals[index].var.getID(), curSum)] = SavedLit(Type::ProgramVar, inputLiterals[index].var);
    } else {
        node = LogicTerm();
        if (!low.deepEquals(true_lit)) {
            formula = formula && (low || !(node));
        }
        if (high.deepEquals(!true_lit)) {
            formula = formula && (!(inputLiterals[index].var) || !(node));
        } else {
            formula = formula && (high || !(inputLiterals[index].var) || !(node));
        }
        history[std::make_pair(inputLiterals[index].var.getID(), curSum)] = SavedLit(Type::AuxVar, node);
    }
    return node;
}
