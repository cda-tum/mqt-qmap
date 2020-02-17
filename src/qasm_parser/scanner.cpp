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

#include <iostream>
#include <fstream>
#include <istream>
#include <map>
#include <wctype.h>
#include <ctype.h>
#include <sstream>

#include "scanner.hpp"

QASMscanner::QASMscanner(std::istream& in_stream) : in(in_stream) {
        // initialize error handling support
        keywords["qreg"] = Token::Kind::qreg;
        keywords["creg"] = Token::Kind::creg;
        keywords["gate"] = Token::Kind::gate;
        keywords["measure"] = Token::Kind::measure;
        keywords["U"] = Token::Kind::ugate;
        keywords["CX"] = Token::Kind::cxgate;
        keywords["pi"] = Token::Kind::pi;
        keywords["OPENQASM"] = Token::Kind::openqasm;
        keywords["show_probabilities"] = Token::Kind::probabilities;
        keywords["sin"] = Token::Kind::sin;
        keywords["cos"] = Token::Kind::cos;
        keywords["tan"] = Token::Kind::tan;
        keywords["exp"] = Token::Kind::exp;
        keywords["ln"] = Token::Kind::ln;
        keywords["sqrt"] = Token::Kind::sqrt;
        keywords["include"] = Token::Kind::include;
        keywords["barrier"] = Token::Kind::barrier;
        keywords["opaque"] = Token::Kind::opaque;
        keywords["if"] = Token::Kind::_if;
        keywords["reset"] = Token::Kind::reset;
        keywords["snapshot"] = Token::Kind::snapshot;
        line = 1;
        col = 0;
        ch = 0;
        nextCh();
}

void QASMscanner::addFileInput(const std::string& fname) {
	auto ifs = new std::ifstream (fname, std::ifstream::in);
	if(ifs->fail()) {
		std::cerr << "Failed to open file '" << fname << "'!" << std::endl;
	} else {
		streams.push(ifs);
		lines.push(LineInfo(ch, line, col));
		line = 0;
		col = 0;
	}
	nextCh();
}

void QASMscanner::nextCh() {
	if(!streams.empty() && streams.top()->eof()) {
		delete streams.top();
		streams.pop();
		ch = lines.top().ch;
		col = lines.top().col;
		line = lines.top().line;
		lines.pop();
		return;
	}
	if(!streams.empty()) {
		col++;
		streams.top()->get(ch);
	} else {
		if(!in.eof()) {
			col++;
			in.get(ch);
		} else {
			ch = (char) -1;
		}
	}
	if(ch == '\n') {
		col = 0;
		line++;
	}
}

Token QASMscanner::next() {
	while(std::isspace(ch)) {
		nextCh();
    }

	Token t = Token(Token::Kind::none, line, col);

	switch(ch) {
        case 'a':
        case 'b':
        case 'c':
        case 'd':
        case 'e':
        case 'f':
        case 'g':
        case 'h':
        case 'i':
        case 'j':
        case 'k':
        case 'l':
        case 'm':
        case 'n':
        case 'o':
        case 'p':
        case 'q':
        case 'r':
        case 's':
        case 't':
        case 'u':
        case 'v':
        case 'w':
        case 'x':
        case 'y':
        case 'z':
        case 'A':
        case 'B':
        case 'C':
        case 'D':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'I':
        case 'J':
        case 'K':
        case 'L':
        case 'M':
        case 'N':
        case 'O':
        case 'P':
        case 'Q':
        case 'R':
        case 'S':
        case 'T':
        case 'U':
        case 'V':
        case 'W':
        case 'X':
        case 'Y':
        case 'Z':
            readName(t);break;
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case '.':
            readNumber(t);break;
        case ';': t.kind = Token::Kind::semicolon; nextCh(); break;
        case (char) -1: t.kind = Token::Kind::eof; break;
        case '(': t.kind = Token::Kind::lpar; nextCh(); break;
        case ')': t.kind = Token::Kind::rpar; nextCh(); break;
        case '[': t.kind = Token::Kind::lbrack; nextCh(); break;
        case ']': t.kind = Token::Kind::rbrack; nextCh(); break;
        case '{': t.kind = Token::Kind::lbrace; nextCh(); break;
        case '}': t.kind = Token::Kind::rbrace; nextCh(); break;
        case ',': t.kind = Token::Kind::comma; nextCh(); break;
        case '+': nextCh(); t.kind = Token::Kind::plus; break;
        case '-': nextCh(); t.kind = Token::Kind::minus; break;
        case '*': nextCh(); t.kind = Token::Kind::times; break;
        case '/': nextCh(); if(ch == '/') {
        						skipComment();
        						t = next();
        					} else {
        						t.kind = Token::Kind::div;
        					}
        					break;
        case '^': nextCh(); t.kind = Token::Kind::power; break;
        case '"': nextCh(); readString(t); nextCh(); break;
        case '>': nextCh(); t.kind = Token::Kind::gt; break;
        case '=': nextCh(); if(ch == '=') {
        						nextCh();
        						t.kind = Token::Kind::eq;
        					} else {
        			            std::cerr << "ERROR: UNEXPECTED CHARACTER: '" << ch << "'! " << std::endl;
        					}
        					break;
        default:
            std::cerr << "ERROR: UNEXPECTED CHARACTER: '" << ch << "'! at line: " << line << ", col: " << col << std::endl;
            nextCh();
        }

        return t;
    }

void QASMscanner::readString(Token& t) {
    	std::stringstream ss;
        while(ch != '"') {
            ss << ch;
            nextCh();
        }
        t.str = ss.str();
        t.kind = Token::Kind::string;
    }

void QASMscanner::skipComment() {
	while(ch != '\n' && ch != (char) -1) {
		nextCh();
	}
}

void QASMscanner::readName(Token& t) {
    	std::stringstream ss;
        while(isdigit(ch) || (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || ch == '_') {
            ss << ch;
            nextCh();
        }
        t.str = ss.str();
        auto it = keywords.find(t.str);
        if(it != keywords.end()) {
            t.kind = it->second;
        } else {
            t.kind = Token::Kind::identifier;
        }
    }

void QASMscanner::readNumber(Token& t) {
        std::stringstream ss;
        while(isdigit(ch)) {
            ss << ch;
            nextCh();
        }
        t.kind = Token::Kind::nninteger;
        t.str = ss.str();
        if(ch != '.') {
        	ss >> t.val;
        	return;
        }
        t.kind = Token::Kind::real;
        ss << ch;
        nextCh();
        while(isdigit(ch)) {
        	ss << ch;
        	nextCh();
        }
        if(ch != 'e' && ch != 'E') {
        	ss >> t.valReal;
        	return;
        }
        ss << ch;
        nextCh();
        if(ch == '-' || ch == '+') {
        	ss << ch;
        	nextCh();
        }
        while(isdigit(ch)) {
        	ss << ch;
        	nextCh();
        }
        ss >> t.valReal;
    }



