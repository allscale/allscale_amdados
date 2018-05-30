//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include <chrono>
#include <iostream>

#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/core/io.h"
#include "allscale/utils/assert.h"

#include "../include/geometry.h"
#include "../include/amdados_utils.h"
#include "../include/configuration.h"
#include "../include/matrix.h"
#include "../include/debugging.h"

namespace amdados {

using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;


// forward declaration - implemented in simulator

point2d_t Glo2CellIndex(const point2d_t & p, const size2d_t & cell_size);

point2d_t Glo2Sub(const point2d_t & p, const size2d_t & cell_size);

void RunDataAssimilation(const Configuration         & conf,
                         const Grid<point_array_t,2> & sensors,
                         const Grid<Matrix,2>        & observations);

void InitialGuess(const Configuration & conf,
                  double_array_t      & x,
                  double_array_t      & y,
                  const point2d_t     & idx);

void OptimizePointLocations(double_array_t & x, double_array_t & y);

namespace {

    void GenerateSensorData(const Configuration& conf, Grid<point_array_t,2>& sensors, Grid<Matrix,2>& observations) {

        // Clear the data structure.
        sensors.forEach([](point_array_t & arr) { arr.clear(); });
        observations.forEach([](Matrix & m) { m.Clear(); });

        // Define useful constants.
        const point2d_t GridSize = GetGridSize(conf);
        const double fraction =
                Bound(conf.asDouble("sensor_fraction"), 0.001, 0.75);

        const size2d_t finest_layer_size(conf.asUInt("subdomain_x"),
                                         conf.asUInt("subdomain_y"));

        // --- initialize sensor positions ---

        // Function scales a coordinate from [0..1] range to specified size.
        auto ScaleCoord = [](double v, int size) -> int {
            int i = static_cast<int>(std::floor(v * size));
            i = std::min(std::max(i, 0), size - 1);
            return i;
        };

        // Global (whole domain) sizes.
        const int Nx = conf.asInt("subdomain_x") * GridSize.x;
        const int Ny = conf.asInt("subdomain_y") * GridSize.y;
        const int problem_size = Nx * Ny;

        // Generate pseudo-random sensor locations.
        const int Nobs = std::max(Round(fraction * problem_size), 1);
        std::cout << "\tfraction of sensor points = " << fraction
                  << ", #observations = " << Nobs << std::endl;
        double_array_t x(Nobs), y(Nobs);
        InitialGuess(conf, x, y, point2d_t(0,0));
        OptimizePointLocations(x, y);

        // Save (scaled) sensor locations.
        for (int k = 0; k < Nobs; ++k) {
            int xk = ScaleCoord(x[k], Nx);
            int yk = ScaleCoord(y[k], Ny);

            // insert sensor position
            point2d_t pt{xk,yk};
            point2d_t idx = allscale::utils::elementwiseDivision(pt,finest_layer_size);
            assert_true((0 <= idx.x) && (idx.x < GridSize.x));
            assert_true((0 <= idx.y) && (idx.y < GridSize.y));
            sensors[idx].push_back(pt % finest_layer_size);
        }


        // --- initialize observations ---
        int Nt = conf.asUInt("Nt");
        for(int x = 0; x < Nx/GridSize.x; x++) {
            for(int y = 0; y < Ny/GridSize.y; y++) {

                // get index of current sensor / observation
                point2d_t idx{x,y};

                // the the positions of the observations
                auto& positions = sensors[idx];

                // insert observation
                auto num_observations = positions.size();
                Matrix& m = observations[idx];
                if (m.Empty()) {
                    m.Resize(Nt,num_observations);
                }

                // if there are no observations, continue processing
                if (num_observations == 0) continue;

                // fill in observation values
                int t_step = Nt / num_observations;
                for(std::size_t cnt=0; cnt < num_observations; cnt++) {
                    m(t_step * cnt,cnt) = 1.0f;
                }
            }
        }

    }

}


/**
 * Function implements a special scenario of Amdados application,
 * where it creates and saves pseudo-randomly distributed sensor locations.
 * N O T E, the sensor points are stored unordered and their coordinates are
 * global, i.e. defined with respected to the whole domain's coordinate system.
 */
void ScenarioBenchmark(const std::string& config_file, int problem_size)
{
    MY_TIME_IT("Running scenario 'benchmark' ...")

    // --- load configuration ---

    // Read the configuration file.
    Configuration conf;
    conf.ReadConfigFile(config_file.c_str());

    // over-ride problem size
    conf.SetInt("num_subdomains_x",problem_size);
    conf.SetInt("num_subdomains_y",problem_size);

    // initialize dependent parameters
    InitDependentParams(conf);
    conf.PrintParameters();

    // print some status info for the user
    int steps = conf.asUInt("Nt");
    std::cout << "Running benchmark based on configuration file \"" << config_file << "\" with domain size " << problem_size << "x" << problem_size << " for " << steps << " time steps ...\n";

    // --- generate sensor data ---
    std::cout << "Generating artificial sensory input data ...\n";

    Grid<point_array_t,2> sensors(GetGridSize(conf));
    Grid<Matrix,2> observations(GetGridSize(conf));
    GenerateSensorData(conf,sensors,observations);

    // --- run simulation ---
    std::cout << "Running benchmark simulation ...\n";
    auto start = std::chrono::high_resolution_clock::now();

	// Run the simulation with data assimilation. Important: by this time
	// some parameters had been initialized in InitDependentParams(..), so
	// we can safely proceed to the main part of the simulation algorithm.
	MY_TIME_IT("Running the simulation with data assimilation ...")
	RunDataAssimilation(conf, sensors, observations);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;

    // --- summarize performance data ---
    double time = (std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() / 1000.0f);
    std::cout << "Simulation took " << time << "s\n";

    double throughput = (problem_size * problem_size * steps) / time;
    std::cout << "Throughput: " << throughput << " sub-domains/s\n";

}


} // namespace amdados

