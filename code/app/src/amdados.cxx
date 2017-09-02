//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------
#include <iostream>
#include <sstream>

// Components to be included
// 1) Grid structures specific to each subdomain
//    Each grid structure contains information
//    a) structures for three different layers of resolution: (100m (1) ; 20m(2); 4m(3))
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
// 7) File reads for initial conditions (simple assci format)
// 8) File reads for flowfields and observation data (larger files in structured format)

// Model Initialization & file read
// Create structures here for read from file following variables:
// Ndom, nelems;

// Data structures initialization for grid, advection diffusion and DA

// Advection diffusion solution
// Data assimilation check on observation
// Switch to appropriate grid resolution
// Data assimilation solver
// File output at periodic intervals

int main(int argc, char ** argv)
{
    std::cout << std::endl << std::endl << std::endl;
    int scenario = 0;
    if (argc > 1) {
        if (argc > 2) {
            std::cout << "ERROR: at most 1 input argument is expected" << std::endl;
            return 1;
        }
        if (!(std::istringstream(argv[1]) >> scenario)) {
            std::cout << "Failed to read scenario ID" << std::endl;
            return 1;
        }
    }
    std::cout << "Scenario: " << scenario << std::endl << std::endl << std::flush;

    switch (scenario) {
        case 0: {
            int Amdados2DMain(void);
            return Amdados2DMain();
        }
        break;
        default: std::cout << "ERROR: unknown scenario: " << scenario << std::endl;
    }
    return 1;
}

