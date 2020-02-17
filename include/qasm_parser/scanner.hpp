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

#ifndef QASMSCANNER_HPP_
#define QASMSCANNER_HPP_

#include <iostream>
#include <fstream>
#include <map>
#include <cctype>
#include <cwctype>
#include <sstream>
#include <stack>
#include <string>


#include "token.hpp"


class QASMscanner {


public:
    explicit QASMscanner(std::istream& in_stream);
    Token next();
    void addFileInput(const std::string& fname);

private:
  	std::istream& in;
  	std::stack<std::istream*> streams;
  	char ch;
  	std::map<std::string, Token::Kind> keywords;
  	int line;
    int col;
    void nextCh();
    void readName(Token& t);
    void readNumber(Token& t);
    void readString(Token& t);
    void skipComment();

    class LineInfo {
    public:
    	char ch;
    	int line, col;

    	LineInfo(char ch, int line, int col) {
    		this->ch = ch;
    		this->line = line;
    		this->col = col;
    	}
    };

  	std::stack<LineInfo> lines;

};

#endif /* QASMSCANNER_HPP_ */
