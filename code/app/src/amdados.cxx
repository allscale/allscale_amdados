//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include <iostream>

#include "amdados/app/debugging.h"

namespace amdados {

void PrintHelp()
{
	std::cout << std::endl << "(c) IBM Research Ireland, 2017-2018";
	std::cout << std::endl << std::endl;
	std::cout << "Help:";
	std::cout << std::endl;
	std::cout << "This program demonstrates solution to advection-diffusion problem";
	std::cout << std::endl;
	std::cout << "with data assimilation. Two separate scenarios are available:";
	std::cout << std::endl;
	std::cout << "1) auxiliary scenario 'sensors' (option: --scenario sensors)";
	std::cout << std::endl;
	std::cout << "   is used to generate random sensors locations.";
	std::cout << std::endl;
	std::cout << "2) default scenario 'simulation' expects the file of sensor";
	std::cout << std::endl;
	std::cout << "   locations and the file of observations generated by Python";
	std::cout << std::endl;
	std::cout << "   solver in the output directory; this is the main scenario";
	std::cout << std::endl;
	std::cout << "   where the advection-diffusion problem is being solved.";
	std::cout << std::endl;
	std::cout << "The same configuration file (default is 'amdados.conf') must be";
	std::cout << std::endl;
	std::cout << "used by either scenario as well as Python scripts for consistency.";
	std::cout << std::endl;
	std::cout << "The option '--config path/to/config_file' allows different";
	std::cout << std::endl;
	std::cout << "configuration files. See README.md for further details.";
	std::cout << std::endl;
	std::cout << "The option '--help' or '-h' prints this help.";
	std::cout << std::endl << std::endl;
}

void ScenarioSimulation(const std::string &);
void ScenarioSensors(const std::string &);
void ScenarioBenchmark(const std::string&, int size);

} // namespace amdados

int main(int argc, char ** argv)
{
    MY_TRY
    {
        std::string scenario = "simulation";
        std::string config_file = "amdados.conf";

        // Parse command-line options.
        for (int a = 0; a < argc; ++a) {
            std::string token = argv[a];
            if (token == "--scenario") {
                if (++a < argc) {
                    scenario = argv[a];
                }
            } else if (token == "--config") {
                if (++a < argc) {
                    config_file = argv[a];
                }
            } else if ((token == "--help") || (token == "-h")) {
                amdados::PrintHelp();
                return EXIT_SUCCESS;
            }
        }

        if (scenario == "sensors") {
            MY_INFO("%s", "SCENARIO: 'sensors'")
            amdados::ScenarioSensors(config_file);
        } else if (scenario.substr(0,9) == "benchmark") {
            MY_INFO("%s", "SCENARIO: 'benchmark'")
            int N = 10;
            if (scenario.size() > 9 && scenario[9] == ':') {
                N = atoi(scenario.c_str() + 10);
            }
            amdados::ScenarioBenchmark(config_file, N);
        } else {
            MY_INFO("%s", "SCENARIO: 'simulation'")
            amdados::ScenarioSimulation(config_file);
        }
        return EXIT_SUCCESS;
    }
    MY_CATCH
    return EXIT_FAILURE;
}

