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

#include "token.hpp"

std::map<Token::Kind, std::string> Token::KindNames = {
			{Kind::none,"none"},
			{Kind::include,"include"},
			{Kind::identifier,"<identifier>"},
			{Kind::number,"<number>"},
			{Kind::plus,"+"},
			{Kind::semicolon,";"},
			{Kind::eof,"EOF"},
			{Kind::lpar,"("},
			{Kind::rpar,")"},
			{Kind::lbrack,"["},
			{Kind::rbrack,"]"},
			{Kind::lbrace,"{"},
			{Kind::rbrace,"}"},
			{Kind::comma,","},
			{Kind::minus,"-"},
			{Kind::times,"*"},
			{Kind::nninteger,"<nninteger>"},
			{Kind::real,"<real>"},
			{Kind::qreg,"qreg"},
			{Kind::creg,"creg"},
			{Kind::ugate,"U"},
			{Kind::cxgate,"CX"},
			{Kind::gate,"gate"},
			{Kind::pi,"pi"},
			{Kind::measure,"measure"},
			{Kind::openqasm,"openqasm"},
			{Kind::probabilities,"probabilities"},
			{Kind::opaque,"opaque"},
			{Kind::sin,"sin"},
			{Kind::cos,"cos"},
			{Kind::tan,"tan"},
			{Kind::exp,"exp"},
			{Kind::ln,"ln"},
			{Kind::sqrt,"sqrt"},
			{Kind::div,"/"},
			{Kind::power,"^"},
			{Kind::string,"string"},
			{Kind::gt,">"},
			{Kind::barrier,"barrier"},
			{Kind::_if,"if"},
			{Kind::eq,"=="},
			{Kind::reset,"reset"}
	};
