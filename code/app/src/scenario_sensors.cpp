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

// amdados_utils.cpp:
point2d_t GetGridSize(const Configuration & conf);

namespace {

using ::allscale::api::user::data::Grid;

/**
 * Function evaluates objective function and its gradient.
 */
void EvaluateObjective(double & J, double_array_t & gradJ,
                             const double_array_t & x,
                             const double_array_t & y)
{
    const int N = static_cast<int>(x.size());   // short-hand alias
    const int NN = N * N;
    assert_true(x.size() == y.size());

    const double EPS = std::sqrt(std::numeric_limits<double>::epsilon());

    J = 0.0;
    gradJ.resize(2*N);
    std::fill(gradJ.begin(), gradJ.end(), 0.0);

    for (int i = 0; i < N; ++i) {
        // Reciprocal distances to subdomain borders.
        const double r_x1 = 1.0 / (std::pow(      x[i],2) + EPS);
        const double r_x2 = 1.0 / (std::pow(1.0 - x[i],2) + EPS);
        const double r_y1 = 1.0 / (std::pow(      y[i],2) + EPS);
        const double r_y2 = 1.0 / (std::pow(1.0 - y[i],2) + EPS);

        J += (r_x1 + r_x2 +
              r_y1 + r_y2);

        double gx = 0.0, gy = 0.0;
        for (int j = 0; j < N; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double sqdist = dx*dx + dy*dy + EPS;
            J  += 1.0 / sqdist;
            gx -= dx / std::pow(sqdist,2);
            gy -= dy / std::pow(sqdist,2);
        }
        gradJ[i  ] = 2.0 * (gx - x[i]  * std::pow(r_x1,2) +
                          (1.0 - x[i]) * std::pow(r_x2,2));
        gradJ[i+N] = 2.0 * (gy - y[i]  * std::pow(r_y1,2) +
                          (1.0 - y[i]) * std::pow(r_y2,2));
    }

    J /= NN;
    std::transform(gradJ.begin(), gradJ.end(), gradJ.begin(),
                                    [=](double x){ return x/NN; });
}

/**
 * Generates initial space distribution of sensor points.
 */
void InitialGuess(const Configuration & conf,
                  double_array_t      & x,
                  double_array_t      & y,
                  const point2d_t     & idx)
{
    std::mt19937 gen(RandomSeed() +
                        idx[_X_] * conf.asInt("num_subdomains_y") + idx[_Y_]);
    std::uniform_real_distribution<double> distrib(0.0, 1.0);
    std::generate(x.begin(), x.end(), [&](){ return distrib(gen); });
    std::generate(y.begin(), y.end(), [&](){ return distrib(gen); });
}

/**
 * Minimizes the objective function by gradient descent.
 */
void OptimizePointLocations(double_array_t & x, double_array_t & y)
{
    const int N = static_cast<int>(x.size());   // short-hand alias
    assert_true(x.size() == y.size());

    const double DOWNSCALE = 0.1;
    const double INITIAL_STEP = 0.1;
    const double TOL = std::numeric_limits<double>::epsilon() * std::log(N);

    double J = 0.0, J_new = 0.0;    // values of objective function
    double step = INITIAL_STEP;     // step in gradient descent

    double_array_t x_new(N), y_new(N);           // sensors' coordinates
    double_array_t gradJ(2*N), gradJ_new(2*N);   // gradients of J

    EvaluateObjective(J, gradJ, x, y);
    for (bool proceed = true; proceed && (step > TINY);) {
        bool is_inside = true;
        for (int k = 0; (k < N) && is_inside; ++k) {
            x_new[k] = x[k] - step * gradJ[k  ];
            y_new[k] = y[k] - step * gradJ[k+N];
            is_inside = ((0.0 <= x_new[k]) && (x_new[k] <= 1.0) &&
                         (0.0 <= y_new[k]) && (y_new[k] <= 1.0));
        }
        if (!is_inside) {
            step *= DOWNSCALE;
            continue;
        }
        EvaluateObjective(J_new, gradJ_new, x_new, y_new);
        if (J < J_new) {
            step *= DOWNSCALE;
            continue;
        }
        proceed = (J - J_new > J * TOL);
        x = x_new;
        y = y_new;
        J = J_new;
        gradJ = gradJ_new;
        step *= 2.0;
    }
}

} // anonymous namespace

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
    using ::allscale::api::user::algorithm::pfor;

    // Read configuration file.
    Configuration conf;
    conf.ReadConfigFile(config_file.c_str());

    // Define useful constants.
    const point2d_t GridSize = GetGridSize(conf);
    const double SensorFraction =
                        Bound(conf.asDouble("sensor_fraction"), 0.01, 0.75);

    // Open file manager and the output file for writing.
    std::string filename = MakeFileName(conf, "sensors");
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto out = manager.openOutputStream(e);

    // Save sensor locations. Note, order is not guaranteed.
    pfor(point2d_t(0,0), GridSize, [&](const point2d_t & idx) {
        // Generate pseudo-random sensor locations.
        const int Nobs = std::max(Round(SensorFraction * SUB_PROBLEM_SIZE), 1);
        double_array_t x(Nobs), y(Nobs);    // sensors' coordinates
        InitialGuess(conf, x, y, idx);
        OptimizePointLocations(x, y);

        // Save pseudo-randomly distributed sensor locations.
        for (int k = 0; k < Nobs; ++k) {
            // Scale coordinates from [0..1] range to subdomain sizes.
            int xk = static_cast<int>(std::floor(x[k] * SUBDOMAIN_X));
            int yk = static_cast<int>(std::floor(y[k] * SUBDOMAIN_Y));

            xk = std::min(std::max(xk, 0), SUBDOMAIN_X - 1);
            yk = std::min(std::max(yk, 0), SUBDOMAIN_Y - 1);

            // Make global coordinates and save.
            xk = Sub2GloX(idx, xk);
            yk = Sub2GloY(idx, yk);
            out.atomic([=](auto & file) { file << xk << " " << yk << "\n"; });
        }
    });
    manager.close(out);
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

    // Define useful constants.
    const point2d_t GridSize = GetGridSize(conf);
    assert_true(sensors.size() == GridSize);

    // Clear the data structure.
    sensors.forEach([](point_array_t & arr) { arr.clear(); });

    // Read the sensor file sequentially.
    std::string filename = MakeFileName(conf, "sensors");
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto in = manager.openInputStream(e);
    while (1) {
        point2d_t pt;
        in.atomic([&](auto & file) { file >> pt.x >> pt.y; });
        if (!in) {
            break;
        }
        point2d_t idx = Glo2CellIndex(pt);
        assert_true((0 <= idx.x) && (idx.x < GridSize.x));
        assert_true((0 <= idx.y) && (idx.y < GridSize.y));
        sensors[idx].push_back(Glo2Sub(pt));
    }
    manager.close(in);

#ifdef AMDADOS_DEBUGGING
    size_t num_sensors = 0;
    sensors.forEach([&](point_array_t & arr) { num_sensors += arr.size(); });
    MY_INFO("Average number of sensors per subdomain: %f",
            (double(num_sensors) / double(GridSize[0] * GridSize[1])));
#endif
}

/**
 * Function sequentially (!) reads the file of sensor locations.
 */
void LoadSensorMeasurements(const Configuration         & conf,
                            const Grid<point_array_t,2> & sensors,
                            Grid<Matrix,2>              & observations)
{
    MY_TIME_IT("Loading sensor measurements ...")
    MY_INFO("%s", "---------------------------------------------------------")
    MY_INFO("%s", "B E W A R E: if you had run: amdados --scenario sensors")
    MY_INFO("%s", "then you have to rerun:")
    MY_INFO("%s", "        python3 python/ObservationsGenerator.py")
    MY_INFO("%s", "otherwise there will be a mismatch error between")
    MY_INFO("%s", "the new sensors and the old observation locations.")
    MY_INFO("%s", "---------------------------------------------------------")

    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;
    using ::allscale::api::user::algorithm::pfor;

    // Define useful constants.
    const point2d_t GridSize = GetGridSize(conf);
    const int       Nt = conf.asInt("Nt");

    Grid<int,2> counters(GridSize);
    int         last_timestamp = -1;

    assert_true(sensors.size() == GridSize);
    assert_true(observations.size() == GridSize);

    // Clear the data structure.
    observations.forEach([](Matrix & m) { m.Clear(); });

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
        counters.forEach([](int & c) { c = 0; });

        // Read all the records of the time-slice.
        for (int i = 0; i < num; ++i) {
            // Read the global coordinates of a sensor and a measurement.
            in.atomic([&](auto & file) { file >> pt.x >> pt.y >> val; });

            // Get subdomain position on the grid.
            point2d_t idx = Glo2CellIndex(pt);
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
            assert_true(Glo2Sub(pt) == sensors[idx][cnt]);

            // Save the measurement in the data matrix.
            m(t,cnt) = val;
            ++cnt;
        }
        assert_true(in);
        last_timestamp = t;

        // Check that all the entries of data matrix were set given a timestamp.
        pfor(point2d_t(0,0), GridSize, [&](const point2d_t & idx) {
            assert_true(counters[idx] == static_cast<int>(sensors[idx].size()));
        });
    }
    manager.close(in);
    assert_true(last_timestamp + 1 == Nt);
}

} // namespace amdados

