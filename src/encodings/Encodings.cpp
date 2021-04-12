
#include "encodings/Encodings.hpp"

std::map<std::pair<unsigned long, long>, SavedLit> history;

expr NaiveExactlyOne(const std::vector<z3::expr>& vars, const std::vector<NestedVar>& varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (const auto& var : varIDs)
	{
		clauseVars.emplace_back(vars[var.varID]);
	}
	return NaiveExactlyOne(clauseVars, c);
}

expr NaiveExactlyOne(const std::vector<z3::expr>& vars, const std::vector<unsigned long>& varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (const auto& var : varIDs)
	{
		clauseVars.emplace_back(vars[var]);
	}
	return NaiveExactlyOne(clauseVars, c);
}

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

expr NaiveAtMostOne(const std::vector<z3::expr>& vars, const std::vector<NestedVar>& varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (const auto& var : varIDs)
	{
		clauseVars.emplace_back(vars[var.varID]);
	}
	return NaiveAtMostOne(clauseVars, c);
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

expr ExactlyOneBiMander(const std::vector<z3::expr>& vars, const std::vector<unsigned long>& varIDs, expr_vector& auxvars, z3::context& c)
{
	std::vector<NestedVar> nVars;
	nVars.reserve(varIDs.size());
	for (const auto& id : varIDs)
	{
		nVars.emplace_back(id);
	}
	return ExactlyOneCMDR(vars, groupVars(nVars, 3), -1, auxvars, c);
}

expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, expr_vector& auxvars, z3::context& c)
{
	return ExactlyOneCMDR(vars, groupVars(vars, 3), -1, auxvars, c);
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

expr AtMostOneCMDR(const std::vector<z3::expr>& vars, expr_vector& auxvars, z3::context& c)
{
	return AtMostOneCMDR(vars, groupVars(vars, 3), -1, auxvars, c);
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

std::vector<NestedVar> groupVars(const std::vector<NestedVar>& vars, std::size_t maxSize)
{
	if (vars.size() <= 6) //Since for n<=5 commander is no faster
		return vars;
	return groupVarsAux(vars, maxSize);
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

std::vector<std::vector<unsigned long>> groupVarsBimander(const expr_vector& vars, std::size_t groupCount)
{
	std::vector<std::vector<unsigned long>> result;
	std::vector<unsigned long> v;
	std::size_t maxSize = std::ceil(vars.size() / double(groupCount));
	size_t i = 0;
	while (i < vars.size())
	{
		if (v.size() == maxSize)
		{
			result.emplace_back(v.begin(), v.end());
			v.clear();
		}
		v.emplace_back(i++);
	}
	if (!v.empty())
		result.emplace_back(v);
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
		if (l.type == 0)
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
		history[std::make_pair(inputLiterals[index].varID, curSum)] = SavedLit(0, inputLiterals[index].varID);
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
		history[std::make_pair(inputLiterals[index].varID, curSum)] = SavedLit(1, auxVars.size() - 1);
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
