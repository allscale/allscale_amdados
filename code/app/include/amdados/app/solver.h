#pragma once

namespace amdados {
namespace app {

    allscale::api::user::data::GridPoint<3> size_grd = {10, 10,2};
    allscale::api::user::data::Grid<double,3> obsv_glob(size_grd);

} // end namespace app
} // end namespace amdados
