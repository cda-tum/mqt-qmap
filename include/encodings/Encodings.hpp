#ifndef Encodings_hpp
#define Encodings_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <z3++.h>
#include <cmath>
#include <map>
#include <limits>

#include "utils.hpp"

using namespace z3;

struct NestedVar
{
	NestedVar(unsigned long varID) : varID(varID), list() {};
	NestedVar(unsigned long varID, const std::vector<NestedVar>& list) : varID(varID), list(list)
	{
	}
	unsigned long varID = std::numeric_limits<unsigned long>::max();
	std::vector<NestedVar> list;
};

struct WeightedVar
{
	WeightedVar(unsigned long varID, int weight) : varID(varID), weight(weight) {}
	unsigned long varID = std::numeric_limits<unsigned long>::max();
	int weight = 0;
};
inline
bool operator<(const WeightedVar& rhs, const WeightedVar& lhs) {
	return rhs.weight < lhs.weight;
}
inline
bool operator==(const WeightedVar& rhs, const WeightedVar& lhs) {
	return rhs.weight == lhs.weight && rhs.varID == lhs.varID;
}

struct SavedLit
{
	SavedLit() :type(-1), varID(0) {}
	SavedLit(int type, unsigned long varID) : type(type), varID(varID) {}
	int type = -1;
	unsigned long varID = 0;
};

expr varAlloc(expr_vector& auxvars, z3::context& c);

expr AtMostOneCMDR(const std::vector<z3::expr>& vars, expr_vector& auxvars, z3::context& c);
expr AtMostOneCMDR(const std::vector<z3::expr>& vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, z3::context& c);

expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, expr_vector& auxvars, z3::context& c);
expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, z3::context& c);

expr NaiveExactlyOne(const std::vector<z3::expr>& clauseVars, z3::context& c);
expr NaiveExactlyOne(const std::vector<z3::expr>& vars, std::vector<int> varIDs, z3::context& c);
expr NaiveExactlyOne(const std::vector<z3::expr>& vars, std::vector<NestedVar> clauseVars, z3::context& c);

expr NaiveAtMostOne(const std::vector<z3::expr>& clauseVars, z3::context& c);
expr NaiveAtMostOne(const std::vector<z3::expr>& vars, std::vector<int> varIDs, z3::context& c);
expr NaiveAtMostOne(const std::vector<z3::expr>& vars, std::vector<NestedVar> varIDs, z3::context& c);

expr NaiveAtLeastOne(const std::vector<z3::expr>& clauseVars, z3::context& c);

expr AtMostOneBiMander(const std::vector<z3::expr>& vars, std::vector<int> varIDs, expr_vector& auxvars, z3::context& c);
expr ExactlyOneBiMander(const std::vector<z3::expr>& vars, std::vector<int> varIDs, expr_vector& auxvars, z3::context& c);

expr BuildBDD(const std::set<WeightedVar> &inputLiterals, const std::vector<z3::expr>& vars, expr_vector& auxVars, int leq, z3::context& c);
expr BuildBDD(unsigned long index, long curSum, long maxSum, long k, const std::vector<WeightedVar>& inputLiterals, const std::vector<z3::expr>& vars, expr_vector auxVars, expr& formula, expr& true_lit, z3::context& c);

std::vector<NestedVar> groupVars(const std::vector<z3::expr>& vars, int maxSize);
std::vector<NestedVar> groupVars(const std::vector<NestedVar>& vars, int maxSize);
std::vector<NestedVar> groupVarsAux(const std::vector<NestedVar>& vars, int maxSize);

std::vector<std::vector<int>> groupVarsBimander(expr_vector vars, int groupCount);
std::vector<std::vector<int>> groupVarsBimander(std::vector<int> vars, int groupCount);

std::string printBimanderVars(const std::vector<std::vector<int>>& vars);
std::string printNestedVars(const std::vector<NestedVar>& vars, int level = 0);
std::string printWeightedVars(const std::vector<WeightedVar>& wVars, expr_vector vars);
#endif