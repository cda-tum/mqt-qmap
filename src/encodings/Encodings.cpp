
#include "encodings/Encodings.hpp"

std::map<std::pair<unsigned long, long>, SavedLit> history;
std::vector<std::vector<unsigned short>> connections;
std::vector<int> d;
std::vector<bool> visited;
int maxSum;

bool sortWeightedVar(WeightedVar v1, WeightedVar v2) {
	return v1.weight < v2.weight;
}

bool compWeightedVar(WeightedVar v1, WeightedVar v2) {
	return v1.varID == v2.varID && v1.weight == v2.weight;
}


expr NaiveExactlyOne(std::vector<z3::expr> vars, std::vector<NestedVar> varIDs, context& c) {
	std::vector<expr> clauseVars;
	for (auto var : varIDs) {
		clauseVars.push_back(vars[var.varID]);
	}
	return NaiveExactlyOne(clauseVars, c);
}

expr NaiveExactlyOne(std::vector<z3::expr> vars, std::vector<int> varIDs, context& c) {
	std::vector<expr> clauseVars;
	for (auto var : varIDs) {
		clauseVars.push_back(vars[var]);
	}
	return NaiveExactlyOne(clauseVars, c);
}

expr NaiveExactlyOne(std::vector<expr> clauseVars, context& c) {
	return NaiveAtLeastOne(clauseVars, c) and NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtLeastOne(std::vector<expr> clauseVars, context& c) {
	expr naiveAtLeastOne = c.bool_val(false);
	for (expr x : clauseVars)
	{
		naiveAtLeastOne = naiveAtLeastOne or x;
	}
	return naiveAtLeastOne;
}

expr NaiveAtMostOne(std::vector<z3::expr> vars, std::vector<NestedVar> varIDs, context& c) {
	std::vector<expr> clauseVars;
	for (auto var : varIDs) {
		clauseVars.push_back(vars[var.varID]);
	}
	return NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtMostOne(std::vector<z3::expr> vars, std::vector<int> varIDs, context& c) {
	std::vector<expr> clauseVars;
	for (auto var : varIDs) {
		clauseVars.push_back(vars[var]);
	}
	return NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtMostOne(std::vector<expr> clauseVars, context& c) {
	expr naiveAtMostOne = c.bool_val(true);
	for (int i = 0; i < clauseVars.size() - 1; i++) {
		for (int j = i + 1; j < clauseVars.size(); j++)
		{
			naiveAtMostOne = naiveAtMostOne and (not clauseVars[i] or not clauseVars[j]);
		}
	}
	return naiveAtMostOne;
}

expr AtMostOneBiMander(std::vector<z3::expr> vars, std::vector<int> varIDs, expr_vector& auxvars, context& c) {
	std::vector < std::vector <int>> subords = groupVarsBimander(varIDs, varIDs.size() / 2);
	expr ret = c.bool_val(true);
	expr_vector binary_vars(c);
	size_t m = subords.size();
	for (int j = 0; j < ceil(log2(m)); j++)
	{
		binary_vars.push_back(varAlloc(auxvars, c));
	}
	for (int i = 0; i < m; i++) {
		expr binary = c.bool_val(true);
		for (int h = 0; h < subords[i].size(); h++) {
			expr b2 = c.bool_val(true);
			for (int j = 0; j < ceil(log2(m)); j++)
			{
				if (i & 1 << j) {
					b2 = b2 and (not vars[subords[i][h]] or binary_vars[j]);
				}
				else {
					b2 = b2 and (not vars[subords[i][h]] or not binary_vars[j]);
				}
			}
			binary = binary and b2;
		}
		ret = ret and binary and NaiveAtMostOne(vars, subords[i], c);
	}
	//std::printf("%s\n",ret.to_string().c_str());
	return ret;
}


expr ExactlyOneBiMander(std::vector<z3::expr> vars, std::vector<int> varIDs, expr_vector& auxvars, context& c) {
	std::vector<NestedVar> nVars;
	for (int id : varIDs) {
		NestedVar x = {((unsigned long)id)};
		nVars.push_back(x);
	}
	return exactlyOneCMDR(vars, groupVars(nVars, 3), -1, auxvars, c);
}

expr exactlyOneCMDR(std::vector<z3::expr> vars, expr_vector& auxvars, context& c) {
	return exactlyOneCMDR(vars, groupVars(vars, 3), -1, auxvars, c);
}

expr exactlyOneCMDR(std::vector<z3::expr> vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, context& c) {
	expr ret = c.bool_val(true);
	std::vector<expr> clauseVars;
	for (int i = 0; i < subords.size(); i++)
	{
		if (subords[i].varID != ULONG_MAX) {
			clauseVars.push_back(vars[subords[i].varID]);
		}
		else 
		{
			clauseVars.push_back(varAlloc(auxvars, c));
			ret = ret and exactlyOneCMDR(vars, subords[i].list, auxvars.size()-1, auxvars, c);
		}
	}
	if (cmdrVar>=0) {
		clauseVars.push_back(not auxvars[cmdrVar]);
	}
	return ret and NaiveExactlyOne(clauseVars, c);
}

expr atMostOneCMDR(std::vector<z3::expr> vars, expr_vector& auxvars, context& c) {
	return atMostOneCMDR(vars, groupVars(vars, 3), -1, auxvars, c);
}

expr atMostOneCMDR(std::vector<z3::expr> vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, context& c) {
	expr ret = c.bool_val(true);
	std::vector<expr> clauseVars;
	for (int i = 0; i < subords.size(); i++)
	{
		if (subords[i].varID != ULONG_MAX) {
			clauseVars.push_back(vars[subords[i].varID]);
		}
		else
		{
			clauseVars.push_back(varAlloc(auxvars, c));
			ret = ret and atMostOneCMDR(vars, subords[i].list, auxvars.size() - 1, auxvars, c);
		}
	}
	if (cmdrVar >= 0) {
		clauseVars.push_back(not auxvars[cmdrVar]);
	}
	return ret and NaiveAtMostOne(clauseVars, c);
}

std::vector<NestedVar> groupVars(std::vector<NestedVar> vars, int maxSize) {
	if (vars.size() <= 6)
		return vars;
	return groupVarsAux(vars, maxSize);
}

std::vector<NestedVar> groupVars(std::vector<z3::expr> vars, int maxSize) {
	std::vector<NestedVar> vVars;
	for (unsigned long i = 0; i < vars.size(); i++)
	{
		vVars.push_back(NestedVar{ i });
	}
	if (vVars.size() <= 6)
		return vVars;
	return groupVarsAux(vVars, maxSize);
}

std::vector<NestedVar> groupVarsAux(std::vector<NestedVar> vars, int maxSize) {
	size_t numVars = vars.size();
	if (numVars <= maxSize)
		return vars;
	std::vector<NestedVar> ret;
	size_t numGr = numVars / maxSize;
	for (unsigned int i = 0; i < numGr; i++)
	{
		ret.push_back(NestedVar{ ULONG_MAX, std::vector<NestedVar>(vars.begin() + i * numVars / numGr, vars.begin() + (i + 1) * numVars / numGr) });
	}
	return groupVarsAux(ret, maxSize);
}

std::vector<std::vector<int>> groupVarsBimander(std::vector<int> vars, int groupCount) {
	std::vector<std::vector<int>> result;
	auto chunkSize = vars.size() / groupCount;

	for (size_t i = 0; i < vars.size(); i += chunkSize) {
		auto end = std::min(vars.size(), i + chunkSize);
		result.emplace_back(vars.begin() + i, vars.begin() + end);

	}

	return result;
}

std::vector<std::vector<int>> groupVarsBimander(expr_vector vars, int groupCount) {
	std::vector<std::vector<int>> result;
	std::vector<int> v;
	int maxSize = ceil(vars.size() / float(groupCount));
	size_t i = 0;
	while (i < vars.size()) {
		if (v.size() == maxSize)
		{
			result.push_back(std::vector<int>(v.begin(), v.end()));
			v.clear();
		}
		v.push_back(i++);
	}
	if (v.size() != 0)
		result.push_back(v);
	return result;
}
int k = 0;
expr buildBDD(std::vector<WeightedVar> inputLiterals, std::vector<z3::expr> vars, expr_vector& auxVars, int leq, context& c) {
	std::sort(inputLiterals.begin(), inputLiterals.end(), sortWeightedVar);
	inputLiterals.erase(std::unique(inputLiterals.begin(), inputLiterals.end(), compWeightedVar), inputLiterals.end());
	//std::printf(printWeightedVars(inputLiterals, vars).c_str());
	history.clear();
	k = leq;
	long maxSum = 0;
	for (auto l : inputLiterals) {
		maxSum += l.weight;
	}
	expr true_lit = varAlloc(auxVars, c);
	expr formula = varAlloc(auxVars, c);
	expr result = buildBDD(0, 0, maxSum, inputLiterals, vars, auxVars, formula, true_lit, c);
	return result and formula;
}

expr buildBDD(unsigned long index, long curSum, long maxSum, std::vector<WeightedVar> inputLiterals, std::vector<z3::expr> vars, expr_vector auxVars, expr& formula, expr& true_lit, context& c) {
	if (curSum + maxSum < k)
		return true_lit;
	if (curSum >= k)
		return not(true_lit);
	if (history.count(std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)) > 0)
	{
		SavedLit l = history[std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)];
		if (l.type == 0) {
			return not(vars[l.id]);
		}
		else {
			return auxVars[l.id];
		}
	}

	expr high = buildBDD(index + 1, curSum + inputLiterals[index].weight, maxSum - inputLiterals[index].weight, inputLiterals, vars, auxVars, formula, true_lit, c);
	expr low = buildBDD(index + 1, curSum, maxSum - inputLiterals[index].weight, inputLiterals, vars, auxVars, formula, true_lit, c);
	 
	if (eq(high, low))
		return high;

	expr node = c.bool_val(true);

	if (eq(high, not(true_lit)) && eq(low, true_lit)) {
		node = not(vars[inputLiterals[index].varID]);
		history[std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)] = SavedLit(0, inputLiterals[index].varID);
	}
	else {
		node = varAlloc(auxVars, c);
		if (!eq(low, true_lit)) {
			formula = formula and (low or not(node));
		}
		if (eq(high, not(true_lit))) {
			formula = formula and (not(vars[inputLiterals[index].varID]) or not(node));
		}
		else {
			formula = formula and (high or not(vars[inputLiterals[index].varID]) or not(node));
		}
		history[std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)] = SavedLit(1, auxVars.size()-1);
	}
	return node;
}


expr varAlloc(expr_vector& auxvars, context& c) {
	static int nextVar = 0;
	std::stringstream out;
	out << "c_" << nextVar++;
	auxvars.push_back(c.bool_const(out.str().c_str()));
	//std::cout << "c_" << nextVar << std::endl;
	return auxvars[auxvars.size() - 1];
}

std::string printBimanderVars(std::vector<std::vector<int>> vars) {
	std::stringstream out;
	for (auto& vec : vars) {
		out << "[" << std::endl;
		out << "\t";
		for (auto& var : vec) {
			out << var << " ";
		}
		out << std::endl;
		out << "]" << std::endl;
	}

	return out.str();
}


std::string printNestedVars(std::vector<NestedVar> vars, int level) {
	std::stringstream out;

	int num = 1;
	for (auto& var : vars)
	{
		if (var.varID != ULONG_MAX) 
		{
			out << var.varID << "-";
		}
		else
		{
			for (int i = 0; i < level && num>1; i++)
			{
				out << "\t";
			}
			out << " [ " << level << ":" << num++ << std::endl;
			for (int i = 0; i < level+1; i++)
			{
				out << "\t";
			}
			out << printNestedVars(var.list, level + 1);
			for (int i = 0; i < level; i++)
			{
				out << "\t";
			}
			out << " ] " << std::endl;
		}

	}
	out << std::endl;
	return out.str();
}

std::string printWeightedVars(std::vector<WeightedVar> wVars, expr_vector vars) {
	std::stringstream out;
	for (auto var : wVars) {
		out << vars[var.varID].to_string() << " - " << var.weight << std::endl;
	}
	return out.str();
}


int findLongestPath(const CouplingMap cm, int nQubits){
	connections.clear();
	connections.resize(nQubits);
	maxSum = -1;
	for (auto edge: cm){
		connections.at(edge.first).emplace_back(edge.second);
	}
	for (int q=0; q < nQubits; ++q){
		d.clear();
		d.resize(nQubits);
		std::fill(d.begin(), d.end(), 0);
		visited.clear();
		visited.resize(nQubits);
		std::fill(visited.begin(), visited.end(), false);
		findLongestPath(q, 0);
		auto it = std::max_element(d.begin(), d.end());
		if ((*it)>maxSum)
			maxSum = (*it);
	}
	return maxSum;
}

void findLongestPath(unsigned short node, int curSum){
	if (visited.at(node))
		return;
	visited[node] = true;

	if (d.at(node)<curSum)
		d[node] = curSum;
	if (connections.at(node).empty()){
		visited[node]=false;
		return;
	}

	for (auto child: connections.at(node)){
		findLongestPath(child, curSum+1);
	}

	visited[node]=false;
}