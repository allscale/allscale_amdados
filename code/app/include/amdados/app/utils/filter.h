#pragma once
#include "allscale/api/user/data/grid.h"
#include "amdados/app/static_grid.h"

using namespace allscale::api::user;

namespace amdados {
namespace app {
namespace utils {

		// compute model covariance structure as function of exponential distance
		// Follow same structure as that implemented for Runga Kutta solver
		template<size_t SizeX, size_t SizeY>
		void getModelCovar(allscale::utils::grid<double,SizeX,SizeY>& mod_covar) {
			int m = 0;
			int n = 0;
			int a = 1; int R1 = 1; int R2 = 1; // data assimilation parametrization coefficients
			for(int i=0; i<(int)SizeX; i++) {      //SizeX and SizeY are number of elements in X and Y
				n = n + 1;
				for(int j=0; j<(int)SizeX; j++) {
				mod_covar[{i,j}] = a*exp(-  abs(  ((( i - (i + j -n) +2  )/R1)^2) + (((j - i + j - n +2)/R2)^2) )  );
				}
				m = m + 1;
			}
		}

} // end namespace utils
} // end namespace app
} // end namespace amdados
