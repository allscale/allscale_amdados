#pragma once
#include "allscale/api/user/data/grid.h"
#include "amdados/app/static_grid.h"

using namespace allscale::api::user;

namespace amdados {
namespace app {
namespace utils {

		// compute model covariance structure as function of exponential distance
		// Follow same structure as that implemented for Runga Kutta solver
		template<int SizeX, int SizeY>
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


		template<int SizeX, int SizeY>
		void ComputeH(allscale::utils::grid<double,SizeX,SizeY>& HMatrix)
		{
		    for(int i=0; i<(int)SizeX; ++i) {      //SizeX and SizeY are number of elements in X and Y
		        HMatrix[{i,i}] = 1;   // identity matrix to map from equivalent grids
		    }
		}

		template<int SizeX, int SizeY>
		void ComputeR(allscale::utils::grid<double,SizeX,SizeY>& RMatrix) {
		    for(int i=0; i<(int)SizeX; ++i) {      //SizeX and SizeY are number of elements in X and Y
		        //RMatrix[{i,i}] = ((double) rand() / (RAND_MAX)) ;   // Random noise between 0 and 1 for now
		        // Albert: Random noise between 0.5 and 1 for now
		        RMatrix[{i,i}] = 1.0 + static_cast<double>(std::rand()) /
		                               static_cast<double>(RAND_MAX);
		    }
		}

		template<int SizeX>
		void getObservation(allscale::utils::grid<double,SizeX>& obsvloc,data::Grid<double,3>& obsver,
		                    int nx, int ny, int timestep, const data::GridPoint<2>& domain)
		{
		    for(int i = 0; i < nx; i++)
		    {
		        int igl = i + (domain[0])*nx;  // map from local to global indices
		        for(int j = 0; j < ny; j++)
		        {
		            int jgl = j + (domain[1])*ny;
		            obsvloc[{i * nx + j}] = obsver[{igl,jgl,timestep}];
		        }
		    }
		}

} // end namespace utils
} // end namespace app
} // end namespace amdados
