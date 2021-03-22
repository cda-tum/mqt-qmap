#ifndef Encodings_hpp
#define Encodings_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <chrono>
#include <z3++.h>
#include <ciso646>
#include <math.h>
#include <map>

using namespace z3;

struct NestedVar
{
	NestedVar(unsigned long id) {
		varID = id;
	};
	NestedVar(unsigned long id, std::vector<NestedVar> list) {
		this->varID = id;
		this->list = list;
	}
	unsigned long varID = ULONG_MAX;
	std::vector<NestedVar> list;
} ;

struct WeightedVar
{
	WeightedVar(unsigned long varID, int weight) {
		this->varID = varID;
		this->weight = weight;
	}
	unsigned long varID;
	int weight = 0;
};

struct SavedLit {
	SavedLit() {
		this->type = -1;
		this->id = 0;
	}
	SavedLit(int type, unsigned long id) {
		this->type = type;
		this->id = id;
	}
	int type = -1;
	unsigned long id = 0;
};

expr varAlloc(expr_vector& auxvars, context& c);

expr atMostOneCMDR(std::vector<z3::expr> vars, expr_vector& auxvars, context& c);
expr atMostOneCMDR(std::vector<z3::expr> vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, context& c);

expr exactlyOneCMDR(std::vector<z3::expr> vars, expr_vector& auxvars, context& c);
expr exactlyOneCMDR(std::vector<z3::expr> vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, context& c);

std::vector<NestedVar> groupVars(std::vector<expr> vars, int maxSize);
std::vector<NestedVar> groupVars(std::vector<NestedVar> vars, int maxSize);
std::vector<NestedVar> groupVarsAux(std::vector<NestedVar> vars, int maxSize);

expr NaiveExactlyOne(std::vector<z3::expr> clauseVars, context& c);
expr NaiveExactlyOne(std::vector<z3::expr> vars, std::vector<int> varIDs, context& c);
expr NaiveExactlyOne(std::vector<z3::expr> vars, std::vector<NestedVar> clauseVars, context& c);

expr NaiveAtMostOne(std::vector<z3::expr> clauseVars, context& c);
expr NaiveAtMostOne(std::vector<z3::expr> vars, std::vector<int> varIDs, context& c);
expr NaiveAtMostOne(std::vector<z3::expr> vars, std::vector<NestedVar> varIDs, context& c);

expr NaiveAtLeastOne(std::vector<expr> clauseVars, context& c);

expr AtMostOneBiMander(std::vector<z3::expr> vars, std::vector<int> varIDs, expr_vector& auxvars, context& c);
expr ExactlyOneBiMander(std::vector<z3::expr> vars, std::vector<int> varIDs, expr_vector& auxvars, context& c);

expr buildBDD(std::vector<WeightedVar> inputLiterals, expr_vector vars, expr_vector& auxVars, int leq, context& c);
expr buildBDD(unsigned long index, long curSum, long maxSum, std::vector<WeightedVar> inputLiterals, expr_vector vars, expr_vector auxVars, expr& formula, expr& true_lit, context& c);

std::vector<std::vector<int>> groupVarsBimander(expr_vector vars, int groupCount);
std::vector<std::vector<int>> groupVarsBimander(std::vector<int> vars, int groupCount);

std::string printBimanderVars(std::vector<std::vector<int>> vars);
std::string printNestedVars(std::vector<NestedVar> vars, int level = 0);
std::string printWeightedVars(std::vector<WeightedVar> wVars, expr_vector vars);
#endif