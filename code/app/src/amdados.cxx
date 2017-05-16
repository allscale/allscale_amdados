#include <iostream>

#include "amdados/app/answer.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/parameters.h"
#include "amdados/app/solver.h"
#include "amdados/app/utils/filter.h"

#include "allscale/api/user/data/grid.h"


#include "amdados/app/utils/matrix.h"
//#include "amdados/app/cholesky.h"
//#include "amdados/app/kalman_filter.h"


using namespace amdados::app;

int main()
{
    Compute(zero, size_global);

	std::cout << "The answer is " << num_domains_x << std::endl;
	std::cout << " observation data = " << obsv_glob[{1,1,1}] << std::endl;
	std::cout << "size=  " << zero << std::endl;
	std::cout << "active layer=  " << A[{1,1}].getActiveLayer() << std::endl;


    std::vector<int> nums{3, 4, 2, 8, 15, 267};

    auto print = [](const int& n) { std::cout << " " << n; };

    std::cout << "before:";
    for_each(nums.begin(), nums.end(), print);
    {

    }
    std::cout << '\n';

	return EXIT_SUCCESS;

	// Components to be included
	// 1) Grid structures specific to each subdomain
    //    Each grid structure contains information
    //    a) structures for three different layers of resolution: (100m (1) ; 20m(2); 4m(3))
    //    b) Solution on each layer
    //
    // 2) Mechanism to switch from layers 1, 2 & 3
    // 3) Advection diffusion solver translated from api-prototype
    //  ) Boundary synchronizations translated from api
    //    Iterative check on error convergence across each subdomain until error norm of all 4 boundaries < threshold
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



}
