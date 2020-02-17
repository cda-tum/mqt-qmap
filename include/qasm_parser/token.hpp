/*
DD-based simulator by JKU Linz, Austria

Developer: Alwin Zulehner, Robert Wille

With code from the QMDD implementation provided by Michael Miller (University of Victoria, Canada)
and Philipp Niemann (University of Bremen, Germany).

For more information, please visit http://iic.jku.at/eda/research/quantum_simulation

If you have any questions feel free to contact us using
alwin.zulehner@jku.at or robert.wille@jku.at

If you use the quantum simulator for your research, we would be thankful if you referred to it
by citing the following publication:

@article{zulehner2018simulation,
    title={Advanced Simulation of Quantum Computations},
    author={Zulehner, Alwin and Wille, Robert},
    journal={IEEE Transactions on Computer Aided Design of Integrated Circuits and Systems (TCAD)},
    year={2018},
    eprint = {arXiv:1707.00865}
}
*/

#ifndef TOKEN_H_
#define TOKEN_H_

#include <map>
#include <string>

class Token {
 public:

	enum class Kind {include, none, identifier, number, plus, semicolon, eof, lpar, rpar, lbrack, rbrack, lbrace, rbrace, comma, minus, times, nninteger, real, qreg, creg, ugate, cxgate, gate, pi, measure, openqasm, probabilities, sin, cos, tan, exp, ln, sqrt, div, power, string, gt, barrier, opaque, _if, eq, reset, snapshot};

	Token(Kind kind, int line, int col) {
		this->kind = kind;
		this->line = line;
		this->col = col;
		this->val = 0;
		this->valReal = 0.0;
	}

	Token() : Token(Kind::none, 0, 0) {
	}

	static std::map<Kind, std::string> KindNames;
	Kind kind;
	int line;
	int col;
	int val;
	double valReal;
	std::string str;

 };

#endif /* TOKEN_H_ */
