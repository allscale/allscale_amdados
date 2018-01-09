//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <iostream>
#include "../include/debugging.h"

// Components to be included
// 1) Grid structures specific to each subdomain
//    Each grid structure contains information
//    a) structures for three different layers of resolution:
//       (100m (1) ; 20m(2); 4m(3))
//    b) Solution on each layer
//
// 2) Mechanism to switch from layers 1, 2 & 3
// 3) Advection diffusion solver translated from api-prototype
//  ) Boundary synchronizations translated from api
//    Iterative check on error convergence across each subdomain
//    until error norm of all 4 boundaries < threshold
// 4) Data assimilation structures translated from api-prototype
// 5) Matrix operation structures for DA
// 6) Data assimilation solution
// 7) File reads for initial conditions (simple ascii format)
// 8) File reads for flow fields and observation data
//    (larger files in structured format)

// Model Initialization & file read
// Create structures here for read from file following variables:
// Ndom, nelems;

// Data structures initialization for grid, advection diffusion and DA

// Advection diffusion solution
// Data assimilation check on observation
// Switch to appropriate grid resolution
// Data assimilation solver
// File output at periodic intervals

void PrintHelp()
{
    std::cout << "TODO: help" << std::endl;
}

namespace amdados {

void ScenarioSimulation(const std::string &);
void ScenarioSensors(const std::string &);

} // namespace amdados

int main(int argc, char ** argv)
{
    MY_TRY
    {
        std::string scenario;
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
                PrintHelp();
                return EXIT_SUCCESS;
            }
        }

        if (scenario == "sensors") {
            MY_INFO("%s", "SCENARIO: 'sensors'")
            amdados::ScenarioSensors(config_file);
        } else {
            MY_INFO("%s", "SCENARIO: 'simulation'")
            amdados::ScenarioSimulation(config_file);
        }
        return EXIT_SUCCESS;
    }
    MY_CATCH
    return EXIT_FAILURE;
}

