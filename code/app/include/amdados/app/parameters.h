#pragma once
#include "allscale/api/user/data/grid.h"
#include "amdados/app/solver.h"

namespace amdados {
namespace app {

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

    const std::string filename = "..//Observation.txt";
    int observint = 1; // number of timesteps between observation availability


    double ReadObservations(allscale::api::user::data::Grid<double,3>& obsver, const std::string filename, int nobspts_x ,int nobspts_y, int observint)
{
    /// Get points of observations and real values
    // Need to be cognizant of timestep availability
    std::ifstream in;
    in.open(filename);
    if(in.is_open())
    {
        size_t t;
        while(in >> t)   // header time stamp
        {
            for(size_t i = 0; i < nobspts_x; i++)
            {
                for(size_t j = 0; j < nobspts_y; j++)
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


   // ReadObservations(obsv_glob,filename,nelems_glob_x,nelems_glob_y,observint);


} // end namespace app
} // end namespace amdados
