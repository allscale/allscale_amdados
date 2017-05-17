#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <limits>
#include <mutex>
#include <random>
#include <memory>
#include <sstream>
#include <iomanip>


#include "allscale/api/user/data/grid.h"
#include "allscale/api/user/operator/pfor.h"

#include "allscale/utils/assert.h"

#include "amdados/app/parameters.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/ibm_apply_runge_kutta.h"
#include "amdados/app/utils/filter.h"
#include "amdados/app/utils/amdados_utils.h"


using namespace allscale::api::user;

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

	// A utility function to save the current state.
	auto SaveGrid2D = [](const std::string & path, const std::string & title, int t,
						  const data::Grid<sub_domain,2> & grid) {
		std::stringstream ss;
        ss << path << "/" << title << std::setfill('0') << std::setw(5) << t << ".txt";
        std::fstream stateFile(ss.str(), std::ios::out | std::ios::trunc);
		assert_true(stateFile.good()) << "failed to open the file for writing: " + ss.str();

        stateFile << "# Layout: [1] dimensionality, [2..dim+1] sizes per dimension, [dim+2...] values"
                  << std::endl;
        size_t dim = 2;
        stateFile << dim << std::endl << (nelems_x * num_domains_x) << std::endl
                                      << (nelems_y * num_domains_y) << std::endl;
        for (int k = 0; k < nelems_x * num_domains_x; ++k) {
        for (int j = 0; j < nelems_y * num_domains_y; ++j) {
            double value = grid[{k/nelems_x, j/nelems_y}].getLayer<L_100m>()
                               [{k%nelems_x, j%nelems_y}];
            stateFile << value << std::endl;
        }}
        stateFile.flush();
	};

void Compute(data::GridPoint<2>& zero, data::GridPoint<2> size_global)
{
    ReadObservations(obsv_glob,filename,nelems_glob_x,nelems_glob_y);
    A[{spot_x,spot_y}].getLayer<L_100m>()[{8,8}] = spot_density;

	pfor(zero, size_global, [&](const data::GridPoint<2>& pos) {
		// initialize all cells on the 100m resolution
		  A[pos].setActiveLayer(L_100m);
		  // initialize the concentration
		  
		  A[pos].forAllActiveNodes([](double& value) {
			 value = 0.0;        // initialize rho with 0
		  });
	  });

    pfor(zero, size_global, [&](const data::GridPoint<2>& idx) {
        auto& tempvar = P[idx];
        utils::getModelCovar(tempvar);
    });

    // bring in a high concentration at some spot
      A[{spot_x,spot_y}].getLayer<L_100m>()[{8,8}] = spot_density;

      // --- Run Simulation ---
      data::Grid<sub_domain,2> tmp_before_RK(size_global);  // debugging
      data::Grid<sub_domain,2> tmp_after_RK(size_global);   // debugging
      data::Grid<sub_domain,2> tmp_before_KF(size_global);  // debugging
      data::Grid<sub_domain,2> tmp_after_KF(size_global);   // debugging
      data::Grid<sub_domain,2> tmp_obvs(size_global);       // debugging
	  // TODO: assert return value
      (void)std::system("mkdir -p output");     // make the output directory
      for (int t = 0; t <= T; t++) {
          std::cout << "Time = " << t << std::endl;
          SaveGrid2D("output", "state", t, A);

          // print state
          if (output_every_nth_time_step != 0 && t % output_every_nth_time_step == 0) printState(A,t);

          // go from time t to t+1

          // go from time t to t+1

             // compute next time step => store it in B
             pfor(zero, size_global, [&](const data::GridPoint<2>& idx) {

                 // compute next step
                 const auto& cur = A[idx];  //cur is defined as grid A[i,j] and distributed across shared memory
                 auto& res = B[idx];         // hence cur is domain of size cur[xElCount,yElCount]

                 // init result with current state
                 res = cur;

                 assert(utils::CheckNoNan(A[idx].getLayer<L_100m>()));
                 assert(utils::CheckNoNan(B[idx].getLayer<L_100m>()));
                 assert(utils::CheckNoNan(res.getLayer<L_100m>()));

                 double timept = t*delta;
                 int nx = 1;
                 int ny = 1;
                 double flowu = utils::mu1(timept);
                 double flowv = utils::mu2(timept);


                 for (Direction dir : { Up, Down, Left, Right }) {

                        // skip global boarder
                        if (dir == Up    && idx[0] == 0)         continue; // if direction == up and no neighbour to south
                        if (dir == Down  && idx[0] == size_global[0]-1) continue;
                        if (dir == Left  && idx[1] == 0)         continue;
                        if (dir == Right && idx[1] == size_global[1]-1) continue;

                        // obtain the local boundary
                        auto local_boundary = cur.getBoundary(dir);

                        // obtain the neighboring boundary
                        auto remote_boundary =
                            (dir == Up)   ? A[idx + data::GridPoint<2>{-1,0}].getBoundary(Down)  : // remote boundary is bottom strip of neighbour
                            (dir == Down) ? A[idx + data::GridPoint<2>{ 1,0}].getBoundary(Up)    : // remote boundary is top of neighbour
                            (dir == Left) ? A[idx + data::GridPoint<2>{0,-1}].getBoundary(Right) : // remote boundary is left of domain
                                            A[idx + data::GridPoint<2>{0, 1}].getBoundary(Left);

                        // compute local flow in domain to decide if flow in or out of domain
                        if (dir == Down) ny = -1; // flow into domain from below
                        if (dir == Left) nx = -1; // flow into domain from left
                        double flow_boundary =  nx*flowu + ny*flowv;
                        // TODO: scale the boundary vectors to the same resolution

                        // compute updated boundary
                        assert(local_boundary.size() == remote_boundary.size());
                        if (flow_boundary < 0) {  // then flow into domain need to update boundary with neighbour value
                            for(size_t i = 0; i<local_boundary.size(); i++) {
                                // for now, we just take the average
                                // need to update this to account for flow direction (Fearghal)
                                local_boundary[i] = remote_boundary[i];
                            }
                        }
				//	std::cout << "DIRECTION: "
				//	<< (dir == Up ? "Up" : (dir == Down ? "Down" : (dir == Left ? "Left" : "Right"))) << std::endl;
					assert(utils::CheckNoNan(res.getLayer<L_100m>()));
					// update boundary in result
					res.setBoundary(dir,local_boundary);

					assert(utils::CheckNoNan(res.getLayer<L_100m>()));
				}
                assert(utils::CheckNoNan(res.getLayer<L_100m>()));

                // 2) run iterations of runge-kutta
				double error = 1000;
				while (error > 0.1) {           // runge kutta iteration on sub-domain
					error = applyRungeKutta(res,flowu,flowv,delta,stepsize);
				}
                 assert(utils::CheckNoNan(res.getLayer<L_100m>()));


                 if (t % 1 == 0) // For now assimilation at periodic timesteps
                 {
                     // Compute Mapping matrix H to project from observation to model grid [nelems_tot,nelems_tot]
                     utils::ComputeH(H[idx]);
                     // Estimate noise/uncertainty metrics of observations [nelems_tot,nelems_tot]
                     utils::ComputeR(R[idx]);
                     // Extract appropriate observation data for subdomain and time [nelems_tot * nelems_tot]
                     utils::getObservation(Obvs[idx],obsv_glob,nelems_x,nelems_y,t,idx);

                  //   assert(utils::CheckNoNan(Obvs[idx]));

                     tmp_before_KF[idx] = res;
                     A[{spot_x,spot_y}].getLayer<L_100m>()[{8,8}] = spot_density;
                     utils::Reshape1Dto2D<nelems_x, nelems_y>(tmp_obvs[idx].getLayer<L_100m>(), Obvs[idx]);
					 utils::Reshape2Dto1D<nelems_x, nelems_y>(forecast[idx], res.getLayer<L_100m>());
            //         std::cout << " this layer 3=  " <<   Obvs[{spot_x,spot_y}][{8}] <<  std::endl;


                     utils::Reshape2Dto1D<nelems_x, nelems_y>(forecast[idx], res.getLayer<L_100m>());


                     assert(utils::CheckNoNan(H[idx]));
                     assert(utils::CheckNoNan(P[idx]));
                     assert(utils::CheckNoNan(R[idx]));
                     assert(utils::CheckNoNan(Obvs[idx]));
                     assert(utils::CheckNoNan(forecast[idx]));
                     assert(utils::CheckNoNan(BLUE[idx]));
                 }

             });

      }

}

//allscale::api::user::pfor(zero, size, [&](const  allscale::api::user::data::GridPoint<2>& pos)

//  allscale::api::user::pfor(zero, size, [&](const  allscale::api::user::data::GridPoint<2>& pos) {

       // initialize all cells on the 100m resolution
    //   A[pos].setActiveLayer(L_100m);
    //   A[pos].DiscretizeElements();
       // here we compute quadrature for each grid

 //  });
     allscale::api::user::data::GridPoint<3> tempvar = {10, 10,2};
//    allscale::api::user::data::Grid<double,3> obsv_glob(size_grd);

} // end namespace app
} // end namespace amdados
