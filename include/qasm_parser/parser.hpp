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

#ifndef QASM_SIMULATOR_H_
#define QASM_SIMULATOR_H_

#include <vector>
#include <set>

#include "scanner.hpp"
#include "token.hpp"

using fp = double;

static constexpr fp PI = 3.141592653589793238462643383279502884197169399375105820974L;

class QASMparser {
public:
	explicit QASMparser(const std::string& fname);
	virtual ~QASMparser();

	void Parse();

    struct gate {
        int target;
        int control;
        char type[128];
    };

    std::vector<std::vector<gate> > getLayers() {
        return layers;
    }

    unsigned int getNqubits() {
        return nqubits;
    }

    unsigned int getNgates() {
        return ngates;
    }

	std::vector<gate> getGates() {
		return gates;
	}

	void clear();
private:
	class Expr {
	public:
		enum class Kind {number, plus, minus, sign, times, sin, cos, tan, exp, ln, sqrt, div, power, id};
		fp num;
		Kind kind;
		Expr* op1 = nullptr;
		Expr* op2 = nullptr;
		std::string id;

		Expr(Kind kind,  Expr* op1, Expr* op2, const fp& num, const std::string& id) {
			this->kind = kind;
			this->op1 = op1;
			this->op2 = op2;
			this->num = num;
			this->id = id;
		}

		Expr(const Expr& expr) {
			this->kind = expr.kind;
			this->num = expr.num;
			this->id = expr.id;
			if(expr.op1 != nullptr) {
				this->op1 = new Expr(*expr.op1);
			}
			if(expr.op2 != nullptr) {
				this->op2 = new Expr(*expr.op2);
			}
		}

		~Expr() {
			delete op1;
			delete op2;
		}
	};

	class BasisGate {
	public:
		virtual ~BasisGate() = default;
	};

	class Ugate : public BasisGate {
	public:
		Expr* theta;
		Expr* phi;
		Expr* lambda;
		std::string target;

		Ugate(Expr* theta, Expr* phi, Expr* lambda, const std::string& target) {
			this->theta = theta;
			this->phi = phi;
			this->lambda = lambda;
			this->target = target;
		}

		~Ugate() override {
			delete theta;
			delete phi;
			delete lambda;

		}
	};

	class CXgate : public BasisGate {
	public:
		std::string control;
		std::string target;

		CXgate(const std::string& control, const std::string& target) {
			this->control = control;
			this->target = target;
		}
	};

	class CompoundGate {
	public:
		std::vector<std::string> parameterNames;
		std::vector<std::string> argumentNames;
		std::vector<BasisGate*> gates;
		bool opaque;
	};

	class Snapshot {
	public:
		~Snapshot() {
			delete[] probabilities;
			delete[] statevector;
		}

		unsigned long long len;
		double* probabilities;
		std::string* statevector;
		std::map<std::string, double> probabilities_ket;
	};

	void scan();
	void check(Token::Kind expected);

	Token la,t;
	Token::Kind sym = Token::Kind::none;

	std::string fname;
  	std::istream* in;
	QASMscanner* scanner;
	std::map<std::string, std::pair<int ,int> > qregs;
	std::map<std::string, std::pair<int, int*> > cregs;
	std::pair<int, int> QASMargumentQreg();
	std::pair<std::string, int> QASMargumentCreg();
	Expr* QASMexponentiation();
	Expr* QASMfactor();
	Expr* QASMterm();
	Expr* QASMexp();
	void QASMgateDecl();
	void QASMopaqueGateDecl();
	void QASMgate(bool execute = true);
	void QASMqop(bool execute = true);
	void QASMexpList(std::vector<Expr*>& expressions);
	void QASMidList(std::vector<std::string>& identifiers);
	void QASMargsList(std::vector<std::pair<int, int> >& arguments);
	std::set<Token::Kind> unaryops {Token::Kind::sin,Token::Kind::cos,Token::Kind::tan,Token::Kind::exp,Token::Kind::ln,Token::Kind::sqrt};

	std::map<std::string, CompoundGate> compoundGates;
	static Expr* RewriteExpr(Expr* expr, std::map<std::string, Expr*>& exprMap);
	static void printExpr(Expr* expr);

	std::vector<std::vector<gate> > layers;
        std::vector<gate> gates;

	unsigned int nqubits = 0;
    int* last_layer;
    unsigned int ngates = 0;

    void addUgate(int target, const fp& theta, const fp& phi, const fp& lambda);
    void addCXgate(int target, int control);
};

#endif /* QASM_SIMULATOR_H_ */
