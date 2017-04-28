#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/operator/pfor.h"

#include "amdados/app/amdados_grid.h"
#include "amdados/app/parameters.h"

namespace amdados {
namespace app {

// steps for the solver;
// 1) Need to initialize all domains to zero (& all layers I would think L100, L20 and L4)
// 2) Need to distribute work across processors for solution of advection diffusion
//      In this manner it would be relatively straightforward task to implement the solver within each subdomain
// 3) Need to integrate data assimilation schemes from kalman_filter.h and filter.h
//      Data structures are in place and simply needs to be integrated with appropriate data from advection-diffusion process
// 4) Need to integrate adaptive meshing capabilities
//      Right now this links with (3): If data is available then resolve solution to higher fidelity grid (down to 20m) and
//      assimilate at this scale. Then solution returns to 100m grid and solution continues. For time being; advection-diffusion
//      always happens on the 100m grid and we resolve at higher resolution only for data assimilation

// --- Initialize ---
  // pfor(zero, size, [&](const  allscale::api::user::data::GridPoint<2>& pos) {

       // initialize all cells on the 100m resolution
    //   A[pos].setActiveLayer(L_100m);
    //   A[pos].DiscretizeElements();
       // here we compute quadrature for each grid

  // });
//    allscale::api::user::data::GridPoint<3> size_grd = {10, 10,2};
//    allscale::api::user::data::Grid<double,3> obsv_glob(size_grd);

} // end namespace app
} // end namespace amdados
