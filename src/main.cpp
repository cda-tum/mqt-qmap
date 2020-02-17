#include "mapper.hpp"


#if MINIMAL_MAPPING_ENABLED
#include "circuit.hpp"
#endif 

#include <boost/program_options.hpp>
#include <regex>
#include <math.h>

namespace po = boost::program_options;

/**
 * Global variables
 */
architecture   arch;
unsigned long  ngates                  = 0;
unsigned long  current_depth           = 0;
unsigned int   nqubits                 = 0;

std::set<edge>                                                               graph;
std::vector<std::vector<QASMparser::gate>>                                   layers;
unique_priority_queue<node, cleanup_node, node_cost_greater, node_func_less> nodes;

double get_pi_div(double val) {
	if(val == 0) {
		return 0;
	}
	const int precision = 10000;
	return round(M_PI / val * precision) / precision;
}

int main(int argc, char** argv) {
#if MINIMAL_MAPPING_ENABLED
	bool   exact       = false;
#endif
	bool   verbose     = false;
	bool   real_format = false;
	std::string input,  input_coupling;
	std::string output, output_statistics;

	// argument handling
	try {
		po::options_description desc{"Options"};
    	desc.add_options()
			("help,h",                                                        "help screen")
#if MINIMAL_MAPPING_ENABLED
			("exact,i",         po::bool_switch(&exact),                      "exact mapping is used rather than heuristic")
#endif
			("input,i",         po::value<std::string>(&input)->required(),   "input file")
			("output,o",        po::value<std::string>(&output),              "output file                           (only for heuristic)")
			("statistic,s",     po::value<std::string>(&output_statistics),   "output statistics file                (only for heuristic)")
			("coupling_file,c", po::value<std::string>(&input_coupling),      "coupling graph - file                 (only for heuristic)")
			("verbose,v",       po::bool_switch(&verbose),                    "verbose                               (only for heuristic)")
			("real,r",          po::bool_switch(&real_format),                "output the circuit in the real format (only for heuristic)");

		po::positional_options_description p;
		p.add("input",     1);
		p.add("statistic", 1);
		p.add("output",    1);

		po::variables_map vm; 
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		
		if (vm.count("help")) {
			std::cout << desc << std::endl;
			return 0;
		}
		po::notify(vm);
	} catch (const po::error &ex) {
		std::cerr << ex.what() << std::endl;
		exit(ERROR);
	}

#if MINIMAL_MAPPING_ENABLED
dfdfd
	std::cout << "ERREREER" << std::endl;
	if(exact) {
		return exact_mapping(input);
	}
#endif 

	if(verbose) {
		std::cout << "Input:        " << input             << std::endl;
		std::cout << "Output:       " << output            << std::endl;
		std::cout << "Statistic:    " << output_statistics << std::endl;
		std::cout << "CouplingFile: " << input_coupling    << std::endl;
		std::cout << "Verbose:      " << verbose           << std::endl;
	}

	// parsing	
	QASMparser* parser = new QASMparser(input.c_str());
	parser->Parse();

	std::vector<QASMparser::gate> gates = parser->getGates();
	nqubits = parser->getNqubits();
	ngates  = parser->getNgates();

	parser->clear();
	delete parser;

	
	// graph handling
	if(!create_architecture_properties(input_coupling)) {
		std::cout << "Error while generating the graph" << std::endl;
		exit(ERROR);
	}
	if(arch.positions > 0 ? nqubits > (unsigned int)arch.positions : (int)nqubits > arch.positions) {
        std::cerr << "ERROR before mapping: more logical qubits than physical ones!" << std::endl;
        exit(ERROR);
    }

		
	// print infos
	const char* bName = basename(input.c_str());
	if(verbose) {
		std::cout << "Circuit name: " << bName << " (requires " << nqubits << " qubits)" << std::endl;
		std::cout << std::endl;
		std::cout << "Before mapping: "                      << std::endl;
		std::cout << "  elementary gates: " << ngates        << std::endl;
		std::cout << "  depth:            " << layers.size() << std::endl;
	} else {
    	std::cout << bName << ',' << nqubits << ',' << ngates << ',' << layers.size() << ',' << std::flush;
	}
	
	// start mapping algorithm
	clock_t begin_time = clock();

	int                                        total_swaps = 0;	
	circuit_properties                         properties  = create_circuit_properties();
    std::vector<QASMparser::gate>              all_gates;
	std::vector<std::vector<QASMparser::gate>> mapped_circuit;

	mapping(gates, mapped_circuit, all_gates, total_swaps, properties);

	double    time     = double(clock() - begin_time) / CLOCKS_PER_SEC;
	int       depth    = mapped_circuit.size();
	int       cost     = all_gates.size()-total_swaps;
#if SPECIAL_OPT	
	long long workload = workload_cost(properties.workload);
	double    fidelity = fidelity_cost(properties.fidelities);
#else
	long long workload = 0;
	double    fidelity = 0;
#endif
	// print statistics
	if(verbose) {
		std::cout << std::endl << "After mapping (no post mapping optimizations are conducted): " << std::endl;
		std::cout << "  elementary gates: " << cost  << std::endl;
		std::cout << "  depth:            " << depth << std::endl;

		std::cout << "\nThe mapping required " << time << " seconds" << std::endl;

		std::cout << "\nInitial mapping of the logical qubits (q) to the physical qubits (Q) of the IBM QX5 architecture: " << std::endl;

		for(unsigned int i = 0; i < nqubits; i++) {
			std::cout << "  q" << i << " is initially mapped to Q" << properties.locations[i] << std::endl;
		} 
	} else {
    	std::cout << time << ',' << cost << ',' << depth << "," << fidelity << std::endl;
	}

	// dump resulting circuit
	if(!output.empty()) {
		std::ofstream of(output);
		if(real_format) {
			of << ".numvars "   << nqubits << std::endl;
			of << ".variables";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << " q" << i;
			}
			of << std::endl;
			of << ".constants ";
			for(unsigned int i = 0; i < nqubits; i++) {
				of << "0";
			}
			of << std::endl;
			of << ".begin" << std::endl;
			for (std::vector<std::vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
					it != mapped_circuit.end(); it++) {
				std::vector<QASMparser::gate> v = *it;
				for (std::vector<QASMparser::gate>::iterator it2 = v.begin(); it2 != v.end(); it2++) {
					std::string hadamard = "U(pi/2,0,pi)";
					if(it2->control != -1) {
						of << "t2 " << "q" << it2->control << " q" << it2->target << std::endl;
					} else if(hadamard.compare(it2->type) == 0) {
						of << "h1 q" << it2->target << std::endl;
					} else {
						std::string s(it2->type);
						std::regex rgx("U\\(([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+), ([+-]?([0-9]*[.])?[0-9]+)\\).*");
						std::smatch match;
						
						if(std::regex_search(s, match, rgx)) {							
							double theta = stof(match[1]);
							double phi   = stof(match[3]);
							double delta = stof(match[5]);

							double theta_div = get_pi_div(theta); 
							double phi_div   = get_pi_div(phi); 
							double delta_div = get_pi_div(delta);

							// conversion to rotation gates
							if(phi_div == 0) {
								of << "rz1:" << 1                     << " q" << it2->target << std::endl; //1.0 / 3
							} else {
								of << "rz1:" << (int)(phi_div / (1 + 3 * phi_div)) << " q" << it2->target << std::endl;
							}
							of << "rx1:" << 2                               << " q" << it2->target << std::endl;
							if(theta_div == 0) {
								of << "rz1:" << 1                           << " q" << it2->target << std::endl;
							} else {
								of << "rz1:" << (int)(theta_div / (1 + theta_div)) << " q" << it2->target << std::endl;
							}
							of << "rx1:" << 2                               << " q" << it2->target << std::endl;
							if(delta_div != 0) {
								of << "rz1:" << delta_div                   << " q" << it2->target << std::endl;
							}
						}
					}
				}
			}
			of << ".end" << std::endl;
		} else {
			of << "OPENQASM 2.0;"              << std::endl;
			of << "include \"qelib1.inc\";"    << std::endl;
			of << "qreg q[" << nqubits << "];" << std::endl;
			of << "creg c[" << nqubits << "];" << std::endl;

			for (std::vector<std::vector<QASMparser::gate> >::iterator it = mapped_circuit.begin();
				it != mapped_circuit.end(); it++) {
				for (std::vector<QASMparser::gate>::iterator it2 = it->begin(); it2 != it->end(); it2++) {
					of << it2->type << " ";
					if (it2->control != -1) {
						of << "q[" << it2->control << "],";
					}
					of << "q[" << it2->target << "];" << std::endl;
				}
			}
		}
	}

	// store timing
	if(!output_statistics.empty()) {
		std::ofstream ofstat (output_statistics, std::ofstream::app);
		//ofstat << bName << " : " << time << " " << depth << " " << cost << " " << workload << " " << alloc_tries << " " << total_swaps << std::endl;
		ofstat << bName << " : " << time << " " << depth << " " << cost << " " << workload << " " << total_swaps << " " << fidelity << std::endl;
	}

	delete_circuit_properties(properties);
	delete_architecture_properties();


	return 0;
}
