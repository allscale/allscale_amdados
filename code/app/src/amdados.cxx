//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include "amdados/app/amdados_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/sparse_matrix.h"
#include "amdados/app/utils/configuration.h"
#include "amdados/app/model/i_model.h"
#include "amdados/app/model/euler_finite_diff.h"

using namespace amdados::app::utils;

namespace amdados {
namespace app {

void AmdadosSolver(const Configuration & conf);

} // namespace app
} // namespace amdados

int main()
{
    Configuration conf;
    conf.ReadConfigFile("../../amdados.conf");
    conf.PrintParameters();
    MakeDirectory(conf.asCString("output_dir"));
    amdados::app::AmdadosSolver(conf);
    return EXIT_SUCCESS;

//Compute(zero, size_global);
//std::cout << "active layer: " << A[{1,1}].getActiveLayer() << std::endl;  // XXX ???

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
