#pragma once

#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/operator/pfor.h"

#include "amdados/app/amdados_grid.h"
#include "amdados/app/parameters.h"

namespace amdados {
namespace app {

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
