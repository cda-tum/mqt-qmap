#ifndef Encodings_hpp
#define Encodings_hpp

#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#include <z3++.h>

using namespace z3;

struct NestedVar {
    explicit NestedVar(unsigned long varID):
        varID(varID), list(){};
    NestedVar(unsigned long varID, std::vector<NestedVar> list):
        varID(varID), list(std::move(list)) {
    }
    unsigned long          varID = std::numeric_limits<unsigned long>::max();
    std::vector<NestedVar> list;
};

struct WeightedVar {
    WeightedVar(unsigned long varID, int weight):
        varID(varID), weight(weight) {}
    unsigned long varID  = std::numeric_limits<unsigned long>::max();
    int           weight = 0;
};
inline bool operator<(const WeightedVar& rhs, const WeightedVar& lhs) {
    return rhs.weight < lhs.weight;
}
inline bool operator==(const WeightedVar& rhs, const WeightedVar& lhs) {
    return rhs.weight == lhs.weight && rhs.varID == lhs.varID;
}

enum class Type { Uninitialized,
                  AuxVar,
                  ProgramVar };
struct SavedLit {
    SavedLit():
        type(Type::Uninitialized), varID(0) {}
    SavedLit(Type type, unsigned long varID):
        type(type), varID(varID) {}
    Type          type  = Type::Uninitialized;
    unsigned long varID = 0;
};

expr varAlloc(expr_vector& auxvars, z3::context& c);

expr AtMostOneCMDR(const std::vector<z3::expr>& vars, const std::vector<NestedVar>& subords, int cmdrVar, expr_vector& auxvars, z3::context& c);

expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, const std::vector<NestedVar>& subords, int cmdrVar, expr_vector& auxvars, z3::context& c);

expr NaiveExactlyOne(const std::vector<z3::expr>& clauseVars, z3::context& c);

expr NaiveAtMostOne(const std::vector<z3::expr>& clauseVars, z3::context& c);
expr NaiveAtMostOne(const std::vector<z3::expr>& vars, const std::vector<unsigned long>& varIDs, z3::context& c);

expr NaiveAtLeastOne(const std::vector<z3::expr>& clauseVars, z3::context& c);

expr AtMostOneBiMander(const std::vector<z3::expr>& vars, const std::vector<unsigned long>& varIDs, expr_vector& auxvars, z3::context& c);

expr BuildBDD(const std::set<WeightedVar>& inputLiterals, const std::vector<z3::expr>& vars, expr_vector& auxVars, int leq, z3::context& c);
expr BuildBDD(unsigned long index, long curSum, long maxSum, long k, const std::vector<WeightedVar>& inputLiterals, const std::vector<z3::expr>& vars, expr_vector& auxVars, expr& formula, expr& true_lit, z3::context& c);

std::vector<NestedVar> groupVars(const std::vector<z3::expr>& vars, std::size_t maxSize);
std::vector<NestedVar> groupVarsAux(const std::vector<NestedVar>& vars, std::size_t maxSize);

std::vector<std::vector<unsigned long>> groupVarsBimander(const std::vector<unsigned long>& vars, std::size_t groupCount);

std::string printBimanderVars(const std::vector<std::vector<unsigned long>>& vars);
std::string printNestedVars(const std::vector<NestedVar>& vars, int level = 0);
std::string printWeightedVars(const std::vector<WeightedVar>& wVars, const expr_vector& vars);

int  findLongestPath(const CouplingMap cm, int nQubits);
int  findLongestPath(const CouplingMap cm, int nQubits, const std::set<unsigned short>& qubitChoice);
void findLongestPath(unsigned short node, int curSum, const std::vector<std::vector<unsigned short>>& connections, std::vector<int>& d, std::vector<bool>& visited);
#endif
