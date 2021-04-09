
#include "encodings/Encodings.hpp"

std::map<std::pair<unsigned long, long>, SavedLit> history;

expr NaiveExactlyOne(const std::vector<z3::expr>& vars, std::vector<NestedVar> varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (auto var : varIDs)
	{
		clauseVars.emplace_back(vars[var.varID]);
	}
	return NaiveExactlyOne(clauseVars, c);
}

expr NaiveExactlyOne(const std::vector<z3::expr>& vars, std::vector<int> varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (auto var : varIDs)
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
	for (expr x : clauseVars)
	{
		naiveAtLeastOne = naiveAtLeastOne or x;
	}
	return naiveAtLeastOne;
}

expr NaiveAtMostOne(const std::vector<z3::expr>& vars, std::vector<NestedVar> varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (auto var : varIDs)
	{
		clauseVars.emplace_back(vars[var.varID]);
	}
	return NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtMostOne(const std::vector<z3::expr>& vars, std::vector<int> varIDs, z3::context& c)
{
	std::vector<expr> clauseVars;
	clauseVars.reserve(varIDs.size());
	for (auto var : varIDs)
	{
		clauseVars.emplace_back(vars[var]);
	}
	return NaiveAtMostOne(clauseVars, c);
}

expr NaiveAtMostOne(const std::vector<expr>& clauseVars, z3::context& c)
{
	expr naiveAtMostOne = c.bool_val(true);
	for (int i = 0; i < clauseVars.size() - 1; i++)
	{
		for (int j = i + 1; j < clauseVars.size(); j++)
		{
			naiveAtMostOne = naiveAtMostOne and (not clauseVars[i] or not clauseVars[j]);
		}
	}
	return naiveAtMostOne;
}

expr AtMostOneBiMander(const std::vector<z3::expr>& vars, std::vector<int> varIDs, expr_vector& auxvars, z3::context& c)
{
	auto subords = groupVarsBimander(varIDs, varIDs.size() / 2);
	expr ret = c.bool_val(true);
	expr_vector binary_vars(c);
	size_t m = subords.size();
	for (int j = 0; j < ceil(log2(m)); j++)
	{
		binary_vars.push_back(varAlloc(auxvars, c));
	}
	for (int i = 0; i < m; i++)
	{
		expr binary = c.bool_val(true);
		for (int h = 0; h < subords[i].size(); h++)
		{
			expr b2 = c.bool_val(true);
			for (int j = 0; j < ceil(log2(m)); j++)
			{
				if (i & 1 << j)
				{
					b2 = b2 and (not vars[subords[i][h]] or binary_vars[j]);
				}
				else
				{
					b2 = b2 and (not vars[subords[i][h]] or not binary_vars[j]);
				}
			}
			binary = binary and b2;
		}
		ret = ret and binary and NaiveAtMostOne(vars, subords[i], c);
	}
	return ret;
}

expr ExactlyOneBiMander(const std::vector<z3::expr>& vars, std::vector<int> varIDs, expr_vector& auxvars, z3::context& c)
{
	std::vector<NestedVar> nVars;
	nVars.reserve(varIDs.size());
	for (int id : varIDs)
	{
		NestedVar x = { ((unsigned long)id) };
		nVars.emplace_back(x);
	}
	return ExactlyOneCMDR(vars, groupVars(nVars, 3), -1, auxvars, c);
}

expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, expr_vector& auxvars, z3::context& c)
{
	return ExactlyOneCMDR(vars, groupVars(vars, 3), -1, auxvars, c);
}

expr ExactlyOneCMDR(const std::vector<z3::expr>& vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, z3::context& c)
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
			ret = ret and ExactlyOneCMDR(vars, it.list, auxvars.size() - 1, auxvars, c);
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

expr AtMostOneCMDR(const std::vector<z3::expr>& vars, std::vector<NestedVar> subords, int cmdrVar, expr_vector& auxvars, z3::context& c)
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
			ret = ret and AtMostOneCMDR(vars, it.list, auxvars.size() - 1, auxvars, c);
		}
	}
	if (cmdrVar >= 0)
	{
		clauseVars.emplace_back(not auxvars[cmdrVar]);
	}
	return ret and NaiveAtMostOne(clauseVars, c);
}

std::vector<NestedVar> groupVars(const std::vector<NestedVar>& vars, int maxSize)
{
	if (vars.size() <= 6)
		return vars;
	return groupVarsAux(vars, maxSize);
}

std::vector<NestedVar> groupVars(const std::vector<z3::expr>& vars, int maxSize)
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

std::vector<NestedVar> groupVarsAux(const std::vector<NestedVar>& vars, int maxSize)
{
	size_t numVars = vars.size();
	if (numVars <= maxSize)
		return vars;
	std::vector<NestedVar> ret;
	size_t numGr = numVars / maxSize;
	ret.reserve(numGr);
	for (unsigned int i = 0; i < numGr; i++)
	{
		ret.emplace_back(NestedVar{ std::numeric_limits<unsigned long>::max(), std::vector<NestedVar>(vars.begin() + i * numVars / numGr, vars.begin() + (i + 1) * numVars / numGr) });
	}
	return groupVarsAux(ret, maxSize);
}

std::vector<std::vector<int>> groupVarsBimander(std::vector<int> vars, int groupCount)
{
	std::vector<std::vector<int>> result;
	auto chunkSize = vars.size() / groupCount;

	for (size_t i = 0; i < vars.size(); i += chunkSize)
	{
		auto end = std::min(vars.size(), i + chunkSize);
		result.emplace_back(vars.begin() + i, vars.begin() + end);
	}

	return result;
}

std::vector<std::vector<int>> groupVarsBimander(expr_vector vars, int groupCount)
{
	std::vector<std::vector<int>> result;
	std::vector<int> v;
	int maxSize = ceil(vars.size() / float(groupCount));
	size_t i = 0;
	while (i < vars.size())
	{
		if (v.size() == maxSize)
		{
			result.emplace_back(std::vector<int>(v.begin(), v.end()));
			v.clear();
		}
		v.emplace_back(i++);
	}
	if (v.size() != 0)
		result.emplace_back(v);
	return result;
}

expr BuildBDD(std::vector<WeightedVar> inputLiterals, const std::vector<z3::expr>& vars, expr_vector& auxVars, int leq, z3::context& c)
{
	inputLiterals.erase(std::unique(inputLiterals.begin(), inputLiterals.end()), inputLiterals.end());
	history.clear();
	long k = leq;
	long maxSum = 0;
	for (auto l : inputLiterals)
	{
		maxSum += l.weight;
	}
	expr true_lit = varAlloc(auxVars, c);
	expr formula = varAlloc(auxVars, c);
	expr result = BuildBDD(0, 0, maxSum, k, inputLiterals, vars, auxVars, formula, true_lit, c);
	return result and formula;
}

expr BuildBDD(unsigned long index, long curSum, long maxSum, long k, std::vector<WeightedVar> inputLiterals, const std::vector<z3::expr>& vars, expr_vector auxVars, expr& formula, expr& true_lit, z3::context& c)
{
	if (curSum + maxSum < k)
		return true_lit;
	if (curSum >= k)
		return not(true_lit);
	if (history.count(std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)) > 0)
	{
		SavedLit l = history[std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)];
		if (l.type == 0)
		{
			return not(vars[l.varID]);
		}
		else
		{
			return auxVars[l.varID];
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
		history[std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)] = SavedLit(0, inputLiterals[index].varID);
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
		history[std::pair<unsigned long, long>(inputLiterals[index].varID, curSum)] = SavedLit(1, auxVars.size() - 1);
	}
	return node;
}

expr varAlloc(expr_vector& auxvars, z3::context& c)
{
	static int nextVar = 0;
	std::stringstream out;
	out << "c_" << nextVar++;
	auxvars.push_back(c.bool_const(out.str().c_str()));
	return auxvars[auxvars.size() - 1];
}

std::string printBimanderVars(std::vector<std::vector<int>> vars)
{
	std::stringstream out;
	for (auto& vec : vars)
	{
		out << "[" << std::endl;
		out << "\t";
		for (auto& var : vec)
		{
			out << var << " ";
		}
		out << std::endl;
		out << "]" << std::endl;
	}

	return out.str();
}

std::string printNestedVars(std::vector<NestedVar> vars, int level)
{
	std::stringstream out;

	int num = 1;
	for (auto& var : vars)
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

std::string printWeightedVars(std::vector<WeightedVar> wVars, expr_vector vars)
{
	std::stringstream out;
	for (auto var : wVars)
	{
		out << vars[var.varID].to_string() << " - " << var.weight << std::endl;
	}
	return out.str();
}
