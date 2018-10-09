//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <chrono>

#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/core/io.h"
#include "allscale/utils/assert.h"
#include "allscale/utils/vector.h"

#include "amdados/app/amdados_utils.h"
#include "amdados/app/configuration.h"
#include "amdados/app/geometry.h"
#include "amdados/app/matrix.h"
#include "amdados/app/debugging.h"
#include "amdados/app/sensors_generator.h"

namespace amdados {

using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;

// Defined in "scenario_simulation.cpp":
void InitDependentParams(Configuration & conf);
void RunDataAssimilation(const Configuration         & conf,
                         const Grid<point_array_t,2> & sensors,
                         const Grid<Matrix,2>        & observations);

// Defined in "scenario_sensors.cpp":
void OptimizePointLocations(double_array_t & x, double_array_t & y);
void InitialGuess(const Configuration & conf,
                  double_array_t      & x,
                  double_array_t      & y,
                  const point2d_t     & idx);

namespace {

//-----------------------------------------------------------------------------
// Generator or synthetic data for benchmarking.
//-----------------------------------------------------------------------------
void GenerateSensorData(const Configuration   & conf,
		                Grid<point_array_t,2> & sensors,
						Grid<Matrix,2>        & observations)
{
	// Define useful constants.
	const point2d_t GridSize = GetGridSize(conf);

//	const double fraction =
//			Bound(conf.asDouble("sensor_fraction"), 0.001, 0.75);

	const size2d_t finest_layer_size(conf.asUInt("subdomain_x"),
									 conf.asUInt("subdomain_y"));

	// --- initialize sensor positions ---

//	// Function scales a coordinate from [0..1] range to specified size.
//	auto ScaleCoord = [](double v, index_t size) -> index_t {
//		index_t i = static_cast<index_t>(std::floor(v * size));
//		i = std::min(std::max(i, index_t(0)), size - 1);
//		return i;
//	};
//
//	// Global (whole domain) sizes.
//	const auto Nx = conf.asInt("subdomain_x") * GridSize.x;
//	const auto Ny = conf.asInt("subdomain_y") * GridSize.y;
//	const auto problem_size = Nx * Ny;
//
//	// Generate pseudo-random sensor locations.
//	const int Nobs = std::max(Round(fraction * problem_size), 1);
//	std::cout << "\tfraction of sensor points = " << fraction;
//	std::cout << ", #observations = " << Nobs << std::endl;
//	double_array_t x(Nobs), y(Nobs);
//	InitialGuess(conf, x, y, point2d_t(0,0));
//	OptimizePointLocations(x, y);

	// temporary sensor location storage
	std::map< point2d_t, std::vector<point2d_t> > locations;

//	// Save (scaled) sensor locations to temporary storage.
//	for (int k = 0; k < Nobs; ++k) {
//		auto xk = ScaleCoord(x[k], Nx);
//		auto yk = ScaleCoord(y[k], Ny);
//		// insert sensor position
//		point2d_t pt{xk,yk};
//		point2d_t idx = allscale::utils::elementwiseDivision(
//													pt, finest_layer_size);
//		assert_true((0 <= idx.x) && (idx.x < GridSize.x));
//		assert_true((0 <= idx.y) && (idx.y < GridSize.y));
//
//		locations[idx].push_back(pt % finest_layer_size);
//	}

	point_array_t sensor_positions;
	SensorsGenerator().MakeSensors(conf, sensor_positions);

    // Save sensor locations to temporary storage.
	for (size_t k = 0; k < sensor_positions.size(); ++k) {
		auto x = sensor_positions[k].x;
		auto y = sensor_positions[k].y;
		// insert sensor position
		point2d_t pt{x,y};
		point2d_t idx =
		        allscale::utils::elementwiseDivision(pt, finest_layer_size);
		assert_true((0 <= idx.x) && (idx.x < GridSize.x));
		assert_true((0 <= idx.y) && (idx.y < GridSize.y));
		locations[idx].push_back(pt % finest_layer_size);
	}

	// --- initialize observations ---

	const index_t Nt = static_cast<index_t>(conf.asUInt("Nt"));

	assert_eq(sensors.size(), observations.size());
	assert_eq(sensors.size(), GridSize);

	::allscale::api::user::algorithm::pfor({0,0}, GridSize,
			[&sensors, &observations, locations, Nt](const auto & idx) {

		allscale::api::core::sema::needs_write_access_on(sensors[idx]);
		allscale::api::core::sema::needs_write_access_on(observations[idx]);

		// Clear the data structures.
		sensors[idx].clear();
		observations[idx].Clear();

		// Copy sensor positions from temporary to simulation storage.
		if (locations.find(idx) != locations.end()) {
			sensors[idx] = locations.at(idx);
		}

		// Insert observations.
		const auto num_observations = sensors[idx].size();
		Matrix & m = observations[idx];
		m.Resize(Nt, num_observations);

		// Fill in observation values.
		if (num_observations > 0) {
			index_t t_step = Nt / num_observations;
			for (std::size_t cnt = 0; cnt < num_observations; cnt++) {
				m(t_step * cnt, cnt) = 1.0f;
			}
		}
	});
}

}	// anonymous namespace

//-----------------------------------------------------------------------------
// Function implements performance testing scenario for Amdados application.
// @param config_file  name of configuration file.
// @param problem_size number of subdomains in either dimension.
//-----------------------------------------------------------------------------
void ScenarioBenchmark(const std::string & config_file, int problem_size)
{
    MY_TIME_IT("Running scenario 'benchmark' ...")

    // --- load configuration ---

    // Read the configuration file.
    Configuration conf;
    conf.ReadConfigFile(config_file.c_str());

    // over-ride problem size
    conf.SetInt("num_subdomains_x", problem_size);
    conf.SetInt("num_subdomains_y", problem_size);

    // initialize dependent parameters
    InitDependentParams(conf);
    conf.PrintParameters();

    // print some status info for the user
    int steps = conf.asInt("Nt");
	std::cout << "Running benchmark based on configuration file \""
			  << config_file;
	std::cout << "\" with domain size " << problem_size << "x" << problem_size;
	std::cout << " for " << steps << " time steps ...\n";

    // --- generate sensor data ---

    std::cout << "Generating artificial sensory input data ...\n";

    Grid<point_array_t,2> sensors(GetGridSize(conf));
    Grid<Matrix,2> observations(GetGridSize(conf));
    GenerateSensorData(conf, sensors, observations);

    // --- run simulation ---

    std::cout << "Running benchmark simulation ...\n";
    auto start = std::chrono::high_resolution_clock::now();

	// Run the simulation with data assimilation. Important: by this time
	// some parameters had been initialized in InitDependentParams(..), so
	// we can safely proceed to the main part of the simulation algorithm.
	RunDataAssimilation(conf, sensors, observations);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;

    // --- summarize performance data ---

    double time = std::chrono::duration_cast<std::chrono::milliseconds>(
    					duration).count() / 1000.0;
    std::cout << "Simulation took " << time << "s\n";

    double throughput = (problem_size * problem_size * steps) / time;
    std::cout << "Throughput: " << throughput << " sub-domains/s\n";
}

} // namespace amdados

