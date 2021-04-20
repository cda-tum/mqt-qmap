
#include "Encodings.hpp"

std::map<std::pair<unsigned long, long>, SavedLit> history;

expr NaiveExactlyOne(const std::vector<expr>& clauseVars, z3::context& c)
{
	return NaiveAtLeastOne(clauseVars, c) and NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtLeastOne(const std::vector<expr>& clauseVars, z3::context& c)
{
	expr naiveAtLeastOne = c.bool_val(false);
	for (const expr& x : clauseVars)
	{
		naiveAtLeastOne = naiveAtLeastOne or x;
	}
	return naiveAtLeastOne;
}

expr NaiveAtMostOne(const std::vector<z3::expr>& vars, const std::vector<unsigned long>& varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (const auto& var : varIDs)
	{
		clauseVars.emplace_back(vars[var]);
	}
	return NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtMostOne(const std::vector<expr>& clauseVars, z3::context& c)
{
	expr naiveAtMostOne = c.bool_val(true);
	for (std::size_t i = 0; i < clauseVars.size() - 1; i++)
	{
		for (std::size_t j = i + 1; j < clauseVars.size(); j++)
		{
			naiveAtMostOne = naiveAtMostOne and (not clauseVars[i] or not clauseVars[j]);
		}
	}
	return naiveAtMostOne;
}

expr AtMostOneBiMander(const std::vector<z3::expr>& vars, const std::vector<unsigned long>& varIDs, expr_vector& auxvars, z3::context& c)
{
	auto subords = groupVarsBimander(varIDs, varIDs.size() / 2);
	expr ret = c.bool_val(true);
	expr_vector binary_vars(c);
	auto m = subords.size();
	for (int j = 0; j < std::ceil(std::log2(m)); j++)
	{
		binary_vars.push_back(varAlloc(auxvars, c));
	}
	for (std::size_t i = 0; i < m; i++)
	{
		expr binary = c.bool_val(true);
		for (std::size_t h = 0; h < subords[i].size(); h++)
		{
			expr b2 = c.bool_val(true);
			for (std::size_t j = 0; j < static_cast<std::size_t>(std::ceil(std::log2(m))); j++)
			{
				if (i & 1 << j)
				{
					b2 = b2 and (not vars[subords[i][h]] or binary_vars[static_cast<int>(j)]);
				}
				else
				{
					b2 = b2 and (not vars[subords[i][h]] or not binary_vars[static_cast<int>(j)]);
				}
			}
			binary = binary and b2;
		}
		ret = ret and binary and NaiveAtMostOne(vars, subords[i], c);
	}
	return ret;
}

expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, const std::vector<NestedVar>& subords, int cmdrVar, expr_vector& auxvars, z3::context& c)
{
	expr ret = c.bool_val(true);
	std::vector<expr> clauseVars;
	clauseVars.reserve(subords.size());
	for (auto& it : subords)
	{
		if (it.varID != std::numeric_limits<unsigned long>::max())
		{
			clauseVars.emplace_back(vars[it.varID]);
		}
		else
		{
			clauseVars.emplace_back(varAlloc(auxvars, c));
			ret = ret and ExactlyOneCMDR(vars, it.list, static_cast<int>(auxvars.size() - 1), auxvars, c);
		}
	}
	if (cmdrVar >= 0)
	{
		clauseVars.emplace_back(not auxvars[cmdrVar]);
	}
	return ret and NaiveExactlyOne(clauseVars, c);
}

expr AtMostOneCMDR(const std::vector<z3::expr>& vars, const std::vector<NestedVar>& subords, int cmdrVar, expr_vector& auxvars, z3::context& c)
{
	expr ret = c.bool_val(true);
	std::vector<expr> clauseVars;
	clauseVars.reserve(subords.size());
	for (auto& it : subords)
	{
		if (it.varID != std::numeric_limits<unsigned long>::max())
		{
			clauseVars.emplace_back(vars[it.varID]);
		}
		else
		{
			clauseVars.emplace_back(varAlloc(auxvars, c));
			ret = ret and AtMostOneCMDR(vars, it.list, static_cast<int>(auxvars.size() - 1), auxvars, c);
		}
	}
	if (cmdrVar >= 0)
	{
		clauseVars.emplace_back(not auxvars[cmdrVar]);
	}
	return ret and NaiveAtMostOne(clauseVars, c);
}

std::vector<NestedVar> groupVars(const std::vector<z3::expr>& vars, std::size_t maxSize)
{
	std::vector<NestedVar> vVars;
	vVars.reserve(vars.size());
	for (unsigned long i = 0; i < vars.size(); i++)
	{
		vVars.emplace_back(NestedVar{ i });
	}
	if (vVars.size() <= 6)
		return vVars;
	return groupVarsAux(vVars, maxSize);
}

std::vector<NestedVar> groupVarsAux(const std::vector<NestedVar>& vars, std::size_t maxSize)
{
	auto numVars = vars.size();
	if (numVars <= maxSize)
		return vars;
	std::vector<NestedVar> ret;
	size_t numGr = numVars / maxSize;
	ret.reserve(numGr);
	for (unsigned int i = 0; i < numGr; i++)
	{
		auto from = vars.begin();
		std::advance(from, static_cast<long>(i * numVars / numGr));
		auto to = from;
		std::advance(to, static_cast<long>(numVars / numGr));
		ret.emplace_back(std::numeric_limits<unsigned long>::max(), std::vector<NestedVar>(from, to));
	}
	return groupVarsAux(ret, maxSize);
}

std::vector<std::vector<unsigned long>> groupVarsBimander(const std::vector<unsigned long>& vars, std::size_t groupCount)
{
	std::vector<std::vector<unsigned long>> result;
	auto chunkSize = vars.size() / groupCount;

	for (size_t i = 0; i < vars.size(); i += chunkSize)
	{
		auto from = vars.begin();
		std::advance(from, static_cast<long>(i));
		auto to = vars.begin();
		auto end = std::min(vars.size(), i + chunkSize);
		std::advance(to, static_cast<long>(end));
		result.emplace_back(from , to);
	}

	return result;
}

expr BuildBDD(const std::set<WeightedVar> &inputLiterals, const std::vector<z3::expr>& vars, expr_vector& auxVars, int leq, z3::context& c)
{
	std::vector<WeightedVar> literals (inputLiterals.begin(), inputLiterals.end());
	history.clear();
	long k = leq;
	long maxSum = 0;
	for (auto &l : literals)
	{
		maxSum += l.weight;
	}
	expr true_lit = varAlloc(auxVars, c);
	expr formula = varAlloc(auxVars, c);
	expr result = BuildBDD(0, 0, maxSum, k, literals, vars, auxVars, formula, true_lit, c);
	return result and formula;
}

expr BuildBDD(unsigned long index, long curSum, long maxSum, long k, const std::vector<WeightedVar>& inputLiterals, const std::vector<z3::expr>& vars, expr_vector& auxVars, expr& formula, expr& true_lit, z3::context& c)
{
	if (curSum + maxSum < k)
		return true_lit;
	if (curSum >= k)
		return not(true_lit);
	if (history.count({inputLiterals[index].varID, curSum}) > 0)
	{
		const SavedLit& l = history[{inputLiterals[index].varID, curSum}];
		if (l.type == Type::ProgramVar)
		{
			return not(vars[l.varID]);
		}
		else
		{
			return auxVars[static_cast<int>(l.varID)];
		}
	}

	expr high = BuildBDD(index + 1, curSum + inputLiterals[index].weight, maxSum - inputLiterals[index].weight, k, inputLiterals, vars, auxVars, formula, true_lit, c);
	expr low = BuildBDD(index + 1, curSum, maxSum - inputLiterals[index].weight, k, inputLiterals, vars, auxVars, formula, true_lit, c);

	if (eq(high, low))
		return high;

	expr node = c.bool_val(true);

	if (eq(high, not(true_lit)) && eq(low, true_lit))
	{
		node = not(vars[inputLiterals[index].varID]);
		history[std::make_pair(inputLiterals[index].varID, curSum)] = SavedLit(Type::ProgramVar, inputLiterals[index].varID);
	}
	else
	{
		node = varAlloc(auxVars, c);
		if (!eq(low, true_lit))
		{
			formula = formula and (low or not(node));
		}
		if (eq(high, not(true_lit)))
		{
			formula = formula and (not(vars[inputLiterals[index].varID]) or not(node));
		}
		else
		{
			formula = formula and (high or not(vars[inputLiterals[index].varID]) or not(node));
		}
		history[std::make_pair(inputLiterals[index].varID, curSum)] = SavedLit(Type::AuxVar, auxVars.size() - 1);
	}
	return node;
}

expr varAlloc(expr_vector& auxvars, z3::context& c)
{
	static int nextVar = 0;
	std::stringstream out;
	out << "c_" << nextVar++;
	auxvars.push_back(c.bool_const(out.str().c_str()));
	return auxvars[static_cast<int>(auxvars.size() - 1)];
}

std::string printBimanderVars(const std::vector<std::vector<unsigned long>>& vars)
{
	std::stringstream out;
	for (const auto& vec : vars)
	{
		out << "[" << std::endl;
		out << "\t";
		for (const auto& var : vec)
		{
			out << var << " ";
		}
		out << std::endl;
		out << "]" << std::endl;
	}

	return out.str();
}

std::string printNestedVars(const std::vector<NestedVar>& vars, int level)
{
	std::stringstream out;

	int num = 1;
	for (const auto& var : vars)
	{
		if (var.varID != std::numeric_limits<unsigned long>::max())
		{
			out << var.varID << "-";
		}
		else
		{
			for (int i = 0; i < level && num > 1; i++)
			{
				out << "\t";
			}
			out << " [ " << level << ":" << num++ << std::endl;
			for (int i = 0; i < level + 1; i++)
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

std::string printWeightedVars(const std::vector<WeightedVar>& wVars, const expr_vector& vars)
{
	std::stringstream out;
	for (const auto& var : wVars)
	{
		out << vars[static_cast<int>(var.varID)].to_string() << " - " << var.weight << std::endl;
	}
	return out.str();
}

int findLongestPath(const CouplingMap cm, int nQubits){
	std::vector<std::vector<unsigned short>> connections;
	std::vector<int> d;
	std::vector<bool> visited;
	connections.resize(nQubits);
	int maxSum = -1;
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
		findLongestPath(q, 0, connections, d, visited);
		auto it = std::max_element(d.begin(), d.end());
		if ((*it)>maxSum)
			maxSum = (*it);
	}
	return maxSum;
}


int findLongestPath(const CouplingMap cm, int nQubits, const std::set<unsigned short> &qubitChoice){
	std::vector<std::vector<unsigned short>> connections;
	std::vector<int> d;
	std::vector<bool> visited;
	connections.resize(nQubits);
	int maxSum = -1;
	for (auto edge: cm){
		if (qubitChoice.count(edge.first) && qubitChoice.count(edge.second))
			connections.at(edge.first).emplace_back(edge.second);
	}
	for (int q=0; q < nQubits; ++q){
		if (connections.at(q).empty())
			continue;
		d.clear();
		d.resize(nQubits);
		std::fill(d.begin(), d.end(), 0);
		visited.clear();
		visited.resize(nQubits);
		std::fill(visited.begin(), visited.end(), false);
		findLongestPath(q, 0, connections, d, visited);
		auto it = std::max_element(d.begin(), d.end());
		if ((*it)>maxSum)
			maxSum = (*it);
	}
	return maxSum;
}

void findLongestPath(unsigned short node, int curSum, const std::vector<std::vector<unsigned short>> &connections, std::vector<int> &d, std::vector<bool> &visited){
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
		findLongestPath(child, curSum+1, connections, d, visited);
	}

	visited[node]=false;
} 