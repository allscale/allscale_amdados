#pragma once
#include "amdados/app/utils/matrix.h"


using namespace allscale::api::user;

namespace amdados {
namespace app {
namespace utils {

template<size_t NELEMS_X, size_t NELEMS_Y>
bool CheckNoNan(const MATRIX<NELEMS_X,NELEMS_Y> & grid)
{
    // TODO: check grid sizes: must be NELEMS_X by NELEMS_Y
    for (int i = 0; i < static_cast<int>(NELEMS_X); i++) {
        for (int j = 0; j < static_cast<int>(NELEMS_Y); j++) {
            if (std::isnan(grid[{i,j}]))
                return false;
        }
    }
    return true;
}

double mu1(double timestep){
    return  -0.6 * sin(timestep/10 - M_PI) * 0.2;
}

double mu2(double timestep){
    return  -1.2 * sin(timestep/5 - M_PI) * 0.2;
}

} // end namespace utils
} // end namespace app
} // end namespace allscale
