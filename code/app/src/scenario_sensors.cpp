//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>

#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/async.h"
#include "allscale/api/core/io.h"
#include "allscale/utils/assert.h"

#include "amdados/app/amdados_utils.h"
#include "amdados/app/configuration.h"
#include "amdados/app/geometry.h"
#include "amdados/app/matrix.h"
#include "amdados/app/debugging.h"
#include "amdados/app/sensors_generator.h"

namespace amdados {

namespace {

using ::allscale::api::user::data::Grid;

/**
 * Function converts global coordinates to the subdomain index on the grid.
 */
inline point2d_t Glo2CellIndex(const point2d_t & p, const size2d_t & cell_size)
{
    return point2d_t(p.x / cell_size.x,
                     p.y / cell_size.y);
}

/**
 * Function converts global coordinates to the local ones inside a subdomain.
 */
inline point2d_t Glo2Sub(const point2d_t & p, const size2d_t & cell_size)
{
    return point2d_t(p.x % cell_size.x,
                     p.y % cell_size.y);
}

///**
// * Function evaluates objective function and its gradient.
// */
//void EvaluateObjective(double & J, double_array_t & gradJ,
//                             const double_array_t & x,
//                             const double_array_t & y)
//{
//    const int N = static_cast<int>(x.size());   // short-hand alias
//    const int NN = N * N;
//    assert_true(x.size() == y.size());
//
//    const double EPS = std::sqrt(std::numeric_limits<double>::epsilon());
//
//    J = 0.0;
//    gradJ.resize(2*N);
//    std::fill(gradJ.begin(), gradJ.end(), 0.0);
//
//    for (int i = 0; i < N; ++i) {
//        // Reciprocal distances to subdomain borders.
//        const double r_x1 = 1.0 / (std::pow(      x[i],2) + EPS);
//        const double r_x2 = 1.0 / (std::pow(1.0 - x[i],2) + EPS);
//        const double r_y1 = 1.0 / (std::pow(      y[i],2) + EPS);
//        const double r_y2 = 1.0 / (std::pow(1.0 - y[i],2) + EPS);
//
//        J += (r_x1 + r_x2 +
//              r_y1 + r_y2);
//
//        double gx = 0.0, gy = 0.0;
//        for (int j = 0; j < N; ++j) {
//            double dx = x[i] - x[j];
//            double dy = y[i] - y[j];
//            double sqdist = dx*dx + dy*dy + EPS;
//            J  += 1.0 / sqdist;
//            gx -= dx / std::pow(sqdist,2);
//            gy -= dy / std::pow(sqdist,2);
//        }
//        gradJ[i  ] = 2.0 * (gx - x[i]  * std::pow(r_x1,2) +
//                          (1.0 - x[i]) * std::pow(r_x2,2));
//        gradJ[i+N] = 2.0 * (gy - y[i]  * std::pow(r_y1,2) +
//                          (1.0 - y[i]) * std::pow(r_y2,2));
//    }
//
//    J /= NN;
//    std::transform(gradJ.begin(), gradJ.end(), gradJ.begin(),
//                                    [=](double x){ return x/NN; });
//}
//
///**
// * Function scales a coordinate from [0..1] range to specified size.
// */
//inline index_t ScaleCoord(double v, index_t size) {
//    index_t i = static_cast<index_t>(std::floor(v * size));
//    i = std::min(std::max(i, index_t(0)), size - 1);
//    return i;
//}

}   // anonymous namespace

/**
 * Generates initial space distribution of sensor points.
 */
void InitialGuess(const Configuration & conf,
                  double_array_t      & x,
                  double_array_t      & y,
                  const point2d_t     & idx)
{
    std::mt19937_64 gen(RandomSeed() +
                        	idx.x * conf.asInt("num_subdomains_y") + idx.y);
    std::uniform_real_distribution<double> distrib(0.0, 1.0);
    std::generate(x.begin(), x.end(), [&](){ return distrib(gen); });
    std::generate(y.begin(), y.end(), [&](){ return distrib(gen); });
}

///**
// * Minimizes the objective function by gradient descent.
// */
//void OptimizePointLocations(double_array_t & x, double_array_t & y)
//{
//    const int N = static_cast<int>(x.size());   // short-hand alias
//    assert_true(x.size() == y.size());
//
//    const double DOWNSCALE = 0.1;
//    const double INITIAL_STEP = 0.1;
//    const double TOL = std::numeric_limits<double>::epsilon() * std::log(N);
//
//    double J = 0.0, J_new = 0.0;    // values of objective function
//    double step = INITIAL_STEP;     // step in gradient descent
//
//    double_array_t x_new(N), y_new(N);           // sensors' coordinates
//    double_array_t gradJ(2*N), gradJ_new(2*N);   // gradients of J
//
//    EvaluateObjective(J, gradJ, x, y);
//    for (bool proceed = true; proceed && (step > TINY);) {
//        bool is_inside = true;
//        for (int k = 0; (k < N) && is_inside; ++k) {
//            x_new[k] = x[k] - step * gradJ[k  ];
//            y_new[k] = y[k] - step * gradJ[k+N];
//            is_inside = ((0.0 <= x_new[k]) && (x_new[k] <= 1.0) &&
//                         (0.0 <= y_new[k]) && (y_new[k] <= 1.0));
//        }
//        if (!is_inside) {
//            step *= DOWNSCALE;
//            continue;
//        }
//        EvaluateObjective(J_new, gradJ_new, x_new, y_new);
//        if (J < J_new) {
//            step *= DOWNSCALE;
//            continue;
//        }
//        proceed = (J - J_new > J * TOL);
//        x = x_new;
//        y = y_new;
//        J = J_new;
//        gradJ = gradJ_new;
//        step *= 2.0;
//    }
//}

/**
 * Function implements a special scenario of Amdados application,
 * where it creates and saves pseudo-randomly distributed sensor locations.
 * N O T E, the sensor points are stored unordered and their coordinates are
 * global, i.e. defined with respected to the whole domain's coordinate system.
 */
void ScenarioSensors(const std::string & config_file)
{
    MY_TIME_IT("Running scenario 'sensors' ...")

    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;
    using ::allscale::api::user::algorithm::async;

    // Read configuration file.
    Configuration conf;
    conf.ReadConfigFile(config_file.c_str());

#if 1

    point_array_t sensors;
    SensorsGenerator().MakeSensors(conf, sensors);

    // Open file manager and the output file for writing, save sensor locations.
    std::string filename = MakeFileName(conf, "sensors");
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto out = manager.openOutputStream(e);
    for (size_t k = 0; k < sensors.size(); ++k) {
        out << sensors[k].x << " " << sensors[k].y << "\n";
    }
    manager.close(out);

#else
    // Define useful constants.
    const point2d_t GridSize = GetGridSize(conf);
    const double fraction =
            Bound(conf.asDouble("sensor_fraction"), 0.001, 0.75);

    if (conf.IsExist("sensor_per_subdomain") &&
            (conf.asInt("sensor_per_subdomain") > 0)) {
        // Local (subdomain) sizes at the finest resolution.
        const int      Sx = conf.asInt("subdomain_x");
        const int      Sy = conf.asInt("subdomain_y");
        const int      sub_problem_size = Sx * Sy;
        const size2d_t subdomain_size(Sx, Sy);

        // Save sensor locations.
        std::string filename = MakeFileName(conf, "sensors");
        allscale::api::user::algorithm::async(
            [conf,fraction,GridSize,subdomain_size,
             sub_problem_size,Sx,Sy,filename]() {
            // Open file manager and the output file for writing.
            FileIOManager & manager = FileIOManager::getInstance();
            Entry e = manager.createEntry(filename, Mode::Text);
            auto out = manager.openOutputStream(e);

            // Generate pseudo-random sensor locations.
            const int Nobs = std::max(Round(fraction * sub_problem_size), 1);
            double_array_t x(Nobs), y(Nobs);    // sensors' coordinates

            for(int i = 0; i < GridSize.x; ++i) {
            for(int j = 0; j < GridSize.y; ++j) {
                point2d_t idx(i,j);
                InitialGuess(conf, x, y, idx);
                OptimizePointLocations(x, y);

                // Save (scaled) pseudo-randomly distributed sensor locations.
                for(int k = 0; k < Nobs; ++k) {
                    point2d_t loc(ScaleCoord(x[k], Sx), ScaleCoord(y[k], Sy));
                    point2d_t glo = Sub2Glo(loc, idx, subdomain_size);
                    out << glo.x << " " << glo.y << "\n";
                }
            }}
            manager.close(out);
        });
    } else {
        // Global (whole domain) sizes.
        const auto Nx = conf.asInt("subdomain_x") * GridSize.x;
        const auto Ny = conf.asInt("subdomain_y") * GridSize.y;
        const auto problem_size = Nx * Ny;

        // Generate pseudo-random sensor locations.
        const int Nobs = std::max(Round(fraction * problem_size), 1);
        MY_LOG(INFO) << "fraction of sensor points = " << fraction
                     << ", #observations = " << Nobs;
        double_array_t x(Nobs), y(Nobs);
        InitialGuess(conf, x, y, point2d_t(0,0));
        OptimizePointLocations(x, y);

        // Open file manager and the output file for writing.
        std::string filename = MakeFileName(conf, "sensors");
        FileIOManager & manager = FileIOManager::getInstance();
        Entry e = manager.createEntry(filename, Mode::Text);
        auto out = manager.openOutputStream(e);

        // Save (scaled) sensor locations.
        for (int k = 0; k < Nobs; ++k) {
            auto xk = ScaleCoord(x[k], Nx);
            auto yk = ScaleCoord(y[k], Ny);
            out << xk << " " << yk << "\n";
        }
        manager.close(out);
    }
#endif
}

/**
 * Function sequentially (!) reads the file of sensor locations.
 */
void LoadSensorLocations(const Configuration   & conf,
                         Grid<point_array_t,2> & sensors)
{
    MY_TIME_IT("Loading sensor locations ...")

    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;
    using ::allscale::api::user::algorithm::pfor;

    const size2d_t finest_layer_size(conf.asUInt("subdomain_x"),
                                     conf.asUInt("subdomain_y"));

    // Define useful constants.
    const point2d_t GridSize = GetGridSize(conf);
    assert_true(sensors.size() == GridSize);

    // Clear the data structure.
    pfor(point2d_t(0, 0), GridSize, [&sensors](const point2d_t & idx) {
        sensors[idx].clear();
    });

    // Read the sensor file sequentially.
    std::string filename = MakeFileName(conf, "sensors");
    CheckFileExists(conf, filename);
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto in = manager.openInputStream(e);
    while (1) {
        point2d_t pt;
        in.atomic([&](auto & file) { file >> pt.x >> pt.y; });
        if (!in) {
            break;
        }
        point2d_t idx = Glo2CellIndex(pt, finest_layer_size);
        assert_true((0 <= idx.x) && (idx.x < GridSize.x));
        assert_true((0 <= idx.y) && (idx.y < GridSize.y));
        sensors[idx].push_back(Glo2Sub(pt, finest_layer_size));
    }
    manager.close(in);

#ifdef AMDADOS_DEBUGGING
    size_t num_sensors = 0;
    sensors.forEach([&](point_array_t & arr) { num_sensors += arr.size(); });
    MY_LOG(INFO) << "Average number of sensors per subdomain: " <<
            (double(num_sensors) / double(GridSize[0] * GridSize[1]));
#endif
}

/**
 * Function sequentially (!) reads the file of sensor measurements.
 */
void LoadSensorMeasurements(const Configuration         & conf,
                            const Grid<point_array_t,2> & sensors,
                            Grid<Matrix,2>              & observations)
{
    MY_TIME_IT("Loading sensor measurements ...")
    MY_LOG(INFO) << "-------------------------------------------------------\n"
                 << "B E W A R E: if you had run: amdados --scenario sensors\n"
                 << "then you have to rerun:\n"
                 << "        python3 python/ObservationsGenerator.py\n"
                 << "otherwise there will be a mismatch error between\n"
                 << "the new sensors and the old observation locations.\n"
                 << "------------------------------------------------------\n";

    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;
    using ::allscale::api::user::algorithm::pfor;

    const size2d_t finest_layer_size(conf.asInt("subdomain_x"),
                                     conf.asInt("subdomain_y"));

    // Define useful constants.
    const point2d_t GridSize = GetGridSize(conf);
    const int       Nt = conf.asInt("Nt");

    Grid<int,2> counters(GridSize);
    assert_decl(int last_timestamp = -1);

    assert_true(sensors.size() == GridSize);
    assert_true(observations.size() == GridSize);

    // Clear the data structure.
    pfor(point2d_t(0, 0), GridSize, [&observations](const point2d_t & idx) {
        observations[idx].Clear();
    });

    // Read the sensor file sequentially.
    std::string filename = MakeFileName(conf, "analytic");
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto in = manager.openInputStream(e);
    while (1) {
        int t = 0, num = 0;
        point2d_t pt;
        float val = 0.0f;

        // Read the header of a new time-slice (timestamp and num. of records).
        in.atomic([&](auto & file) { file >> t >> num; });
        if (!in) {
            break;
        }
        assert_true(t < Nt);
        assert_true(last_timestamp + 1 == t);

        // Reset all the counters upon arrival of a new time-slice.
        pfor(point2d_t(0,0), GridSize, [&counters](const point2d_t & idx) {
            counters[idx] = 0;
        });

        // Read all the records of the time-slice.
        for (int i = 0; i < num; ++i) {
            // Read the global coordinates of a sensor and a measurement.
            in.atomic([&](auto & file) { file >> pt.x >> pt.y >> val; });

            // Get subdomain position on the grid.
            point2d_t idx = Glo2CellIndex(pt, finest_layer_size);
            assert_true((0 <= idx.x) && (idx.x < GridSize.x));
            assert_true((0 <= idx.y) && (idx.y < GridSize.y));

            // Get the data matrix of the subdomain, and the counter.
            Matrix & m = observations[idx];
            int    & cnt = counters[idx];
            int      nsensors = static_cast<int>(sensors[idx].size());

            // Allocate the data matrix, if necessary.
            if (m.Empty()) {
                m.Resize(Nt, nsensors);
            }

            // Check that sensor coordinates appear in exactly the same order.
            assert_true(cnt < nsensors);
            assert_true(Glo2Sub(pt, finest_layer_size) == sensors[idx][cnt]);

            // Save the measurement in the data matrix.
            m(t,cnt) = val;
            ++cnt;
        }
        assert_true(in);
        assert_decl(last_timestamp = t);

        // Check that all the entries of data matrix were set given a timestamp.
        assert_decl(pfor(point2d_t(0,0), GridSize, [&](const point2d_t & idx) {
            assert_true(counters[idx] == static_cast<int>(sensors[idx].size()));
        }));
    }
    manager.close(in);
    assert_true(last_timestamp + 1 == Nt);
}

} // namespace amdados

