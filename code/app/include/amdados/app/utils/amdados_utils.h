#pragma once
#include "amdados/app/utils/matrix.h"


using namespace allscale::api::user;

namespace amdados {
namespace app {
namespace utils {

			//-------------------------------------------------------------------------------------------------
			// Function checks there is no NAN values on the grid.
			//-------------------------------------------------------------------------------------------------
			template<size_t LENGTH>
			bool CheckNoNan(const VECTOR<LENGTH> & grid)
			{
			    // TODO: check size: must be LENGTH
			    for (int i = 0; i < static_cast<int>(LENGTH); i++) {
			        if (std::isnan(grid[{i}]))
			            return false;
			    }
			    return true;
			}

			double mu1(double timestep){
				return  -0.6 * sin(timestep/10 - M_PI) * 0.2;
			}

			double mu2(double timestep){
				return  -1.2 * sin(timestep/5 - M_PI) * 0.2;
			}

			//-------------------------------------------------------------------------------------------------
			// Function checks there is no NAN values on the grid.
			//-------------------------------------------------------------------------------------------------
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

			//-------------------------------------------------------------------------------------------------
			// Function reshapes a vector into 2D grid structure,
			// in Matlab notation: grid = reshape(vec, [NELEMS_X, NELEMS_Y]).
			//-------------------------------------------------------------------------------------------------
			template<size_t NELEMS_X, size_t NELEMS_Y, typename GRID>
			void Reshape1Dto2D(GRID & grid, const VECTOR<NELEMS_X * NELEMS_Y> & vec)
			{
			    // TODO: check grid sizes: must be NELEMS_X by NELEMS_Y
			    for (int i = 0; i < static_cast<int>(NELEMS_X); i++) {
			        for (int j = 0; j < static_cast<int>(NELEMS_Y); j++) {
			            grid[{i,j}] = vec[{i * static_cast<int>(NELEMS_X) + j}];
			        }
			    }
			}

			//-------------------------------------------------------------------------------------------------
			// Function unrolls 2D grid structure into a vector, in Matlab notation: vec = grid(:).
			//-------------------------------------------------------------------------------------------------
			template<size_t NELEMS_X, size_t NELEMS_Y, typename GRID>
			void Reshape2Dto1D(VECTOR<NELEMS_X * NELEMS_Y> & vec, const GRID & grid)
			{
				// TODO: check grid sizes: must be NELEMS_X by NELEMS_Y
				for (int i = 0; i < static_cast<int>(NELEMS_X); i++) {
					for (int j = 0; j < static_cast<int>(NELEMS_Y); j++) {
						vec[{i * static_cast<int>(NELEMS_X) + j}] = grid[{i,j}];
					}
				}
			}

} // end namespace utils
} // end namespace app
} // end namespace allscale
