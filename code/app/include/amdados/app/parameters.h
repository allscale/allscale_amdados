#pragma once
#include <cstdlib>
#include <iostream>

#include "allscale/api/user/data/grid.h"
#include "amdados/app/amdados_grid.h"


//using namespace allscale::api::user::data;

namespace amdados {
namespace app {


// need to update this section to read these parameters from input file
    const int num_domains_x = 10;
    const int num_domains_y = 10;

    const int nelems_x = 10;
    const int nelems_y = 10;
    const int nelems_glob_x = nelems_x * num_domains_x;
    const int nelems_glob_y = nelems_y * num_domains_y;

    int T = 2;
    const int timestep = T;
    int output_every_nth_time_step = 5;
    double delta = 0.1; // not used yet

    // initial spot
    int spot_x = 0;   // start at bottom left point and then obtain analytical expression of evolution
    int spot_y = 0;
    double spot_density = 10000;

    // set up the configuration of a grid cell (static)
    using sub_domain_config = CellConfig<
            layers<                 //  1000m x 1000m  each sub-domain covers
                layer<nelems_x,nelems_y>,       //  10 x 10  100m nodes each consisting of
                layer<5,5>,         //   5 x  5   20m nodes each consisting of
                layer<5,5>          //   5 x  5    4m nodes
            >
    >;

    //assign each layer a level corresponding to coarsest to finest
    enum {
        L_100m = 2,
        L_20m = 1,
        L_4m = 0,
    };

    allscale::api::user::data::GridPoint<2> zero = 0;
    allscale::api::user::data::GridPoint<2> size_global = {num_domains_x, num_domains_y };

    // create the type of a grid cell
    using sub_domain = Cell<double,sub_domain_config>;
    // create the overall grid
    allscale::api::user::data::Grid<sub_domain,2> A(size_global);  // A is of form A[{ndox,ndomy}].layer[{xElCount,yElcount}]
    allscale::api::user::data::Grid<sub_domain,2> B(size_global);


    const std::string filename = "..//..//Observation.txt";
    int observint = 1; // number of timesteps between observation availability
// all initialization parameters - move to input file


    // data structures for observation data
    allscale::api::user::data::GridPoint<3> size_grd = {nelems_glob_x, nelems_glob_y,timestep + 1};
    allscale::api::user::data::Grid<double,3> obsv_glob(size_grd);


           // initialize all cell on the 100m resolution

// Call to read observation data
void ReadObservations(allscale::api::user::data::Grid<double,3>& obsver, const std::string filename, int nobspts_x ,int nobspts_y)
{
    /// Get points of observations and real values
    // Need to be cognizant of timestep availability
    std::ifstream in;
    in.open(filename);
    if(in.is_open())
    {
        int t;
        while(in >> t)   // header time stamp
        {
            for(int i = 0; i < nobspts_x; i++)
            {
                for(int j = 0; j < nobspts_y; j++)
                {
                    int i1, j1;
                    double val;
                    in >> i1 >> j1 >> val;
                    obsver[{i1,j1,t}] = val;   //store in 3D array; extract 1D slice from appropriate location fo DA in 1D vector with mapping
                }
            }
            if(t  == 0)
            {// right now assume data available every timestep; need to update
            }
        }
        in.close();
    }
    else
    {
        std::cout << "ERROR: observation file: " << filename << "  was not found;" << std::endl;
        throw std::exception();
    }
}

// introduce a function to read flowfield data
// Read of form I,J,U,V and map to grid of amdados application
void ReadFlows(allscale::api::user::data::Grid<double,3>& flowfield, const std::string filename_flow, int nflopts_x ,int nflopts_y)
{


}


   // ReadObservations(obsv_glob,filename,nelems_glob_x,nelems_glob_y,observint);


} // end namespace app
} // end namespace amdados
