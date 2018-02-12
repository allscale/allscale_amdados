#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include <allscale/api/core/io.h>
#include <allscale/utils/assert.h>

#include "../include/debugging.h"
#include "../include/geometry.h"
#include "../include/amdados_utils.h"
#include "../include/configuration.h"
#include "../include/matrix.h"

namespace amdados {

// amdados_utils.cpp:
point2d_t GetOrigin(const Configuration & conf);
point2d_t GetGridSize(const Configuration & conf);

namespace {

using ::allscale::api::user::data::Grid;

typedef std::vector<double> dbl_arr_t;

const double TINY = numeric_limits<double>::min() /
           std::pow(numeric_limits<double>::epsilon(),3);

 /**
  * Function evaluates objective function and its gradient.
  */
void EvaluateObjective(double & J, dbl_arr_t & gradJ,
                             const dbl_arr_t & x,
                             const dbl_arr_t & y)
{
    const int N = NUM_SUBDOMAIN_OBSERVATIONS;       // short-hand alias
    const int NN = N * N;

    const double EPS = std::sqrt(numeric_limits<double>::epsilon());

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
void InitialGuess(dbl_arr_t & x, dbl_arr_t & y, const point2d_t & idx)
{
    std::mt19937 gen(RandomSeed() + idx[_X_] * NUM_DOMAINS_Y + idx[_Y_]);
    std::uniform_real_distribution<double> distrib(0.0, 1.0);
    std::generate(x.begin(), x.end(), [&](){ return distrib(gen); });
    std::generate(y.begin(), y.end(), [&](){ return distrib(gen); });
}

/**
 * Minimizes the objective function by gradient descent.
 */
void OptimizePointLocations(dbl_arr_t & x, dbl_arr_t & y)
{
    const int N = NUM_SUBDOMAIN_OBSERVATIONS;       // short-hand alias

    const double DOWNSCALE = 0.1;
    const double INITIAL_STEP = 0.1;
    const double TOL = numeric_limits<double>::epsilon() * std::log(N);

    double J = 0.0, J_new = 0.0;    // values of objective function
    double step = INITIAL_STEP;     // step in gradient descent

    dbl_arr_t x_new(N), y_new(N);           // sensors' coordinates
    dbl_arr_t gradJ(2*N), gradJ_new(2*N);   // gradients of J

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
int ScenarioSensors(const std::string & config_file)
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
    const point2d_t Origin = GetOrigin(conf);
    const point2d_t GridSize = GetGridSize(conf);

    // Open file manager and the output file for writing.
    std::string filename = MakeFileName(conf, "sensors");
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto out = manager.openOutputStream(e);

    // Save sensor locations. Note, order is not guaranteed.
    pfor(Origin, GridSize, [&](const point2d_t & idx) {
        // Generate pseudo-random sensor locations.
        const int N = NUM_SUBDOMAIN_OBSERVATIONS;   // short-hand alias
        dbl_arr_t x(N), y(N);                       // sensors' coordinates
        InitialGuess(x, y, idx);
        OptimizePointLocations(x, y);

        // Save pseudo-randomly distributed sensor locations.
        for (int k = 0; k < N; ++k) {
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
    return EXIT_SUCCESS;
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
}

/**
 * Function sequentially (!) reads the file of sensor locations.
 */
void LoadSensorMeasurements(const Configuration    & conf,
                            const Grid<point_array_t,2> & sensors,
                            Grid<Matrix,2>         & observations)
{
    MY_TIME_IT("Loading sensor measurements ...")
    MY_INFO("-----------------------------------------------------------------")
    MY_INFO("B E W A R E: if you had run: 'amdados --scenario sensors'")
    MY_INFO("then you have to rerun: 'python3 python/ObservationsGenerator.py'")
    MY_INFO("otherwise there will be a mismatch error between the new sensors")
    MY_INFO("and the old observation locations.")
    MY_INFO("-----------------------------------------------------------------")

    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;
    using ::allscale::api::user::algorithm::pfor;

    // Define useful constants.
    const point2d_t Origin = GetOrigin(conf);
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
        pfor(Origin, GridSize, [&](const point2d_t & idx) {
            assert_true(counters[idx] == static_cast<int>(sensors[idx].size()));
        });
    }
    manager.close(in);
    assert_true(last_timestamp + 1 == Nt);
}

} // namespace amdados










///**
// * Functor compares 2D points using "less" semantic.
// */
//struct CmpPoint2d : public std::binary_function<point2d_t,point2d_t,bool>
//{
//    result_type operator()(const first_argument_type  & a,
//                           const second_argument_type & b)
//    {
//        return ((a.x < b.x) || ((a.x == b.x) && (a.y < b.y)));
//    }
//};

///**
// * Function makes the data structure more compact and sorts sensor locations.
// */
//void CompactAndSort(Grid<point_array_t,2> & data)
//{
//    data.pforEach([](point_array_t & arr) {
//        if (arr.size() < arr.capacity()) {
//            point_array_t tmp = arr;
//            arr.swap(tmp);
//        }
//        std::sort(arr.begin(), arr.end(), CmpPoint2d());
//    });
//}

///**
// * Testing against race condition.
// */
//void Test(const Configuration & conf)
//{
//    using ::allscale::api::core::FileIOManager;
//    using ::allscale::api::core::Entry;
//    using ::allscale::api::core::Mode;
//    using ::allscale::api::user::algorithm::pfor;
//
//    // Define useful constants.
//    const int Nx = conf.asInt("num_subdomains_x") * SUBDOMAIN_X;  // full domain size
//    const int Ny = conf.asInt("num_subdomains_y") * SUBDOMAIN_Y;  // full domain size
//    const point2d_t Origin = {0, 0};
//    const point2d_t GridSize = { conf.asInt("num_subdomains_x"),
//                                 conf.asInt("num_subdomains_y") };
//
//    // Make sensor's filename.
//    std::stringstream filename;
//    filename << conf.asString("output_dir") << PathSep << "sensors"
//             << "_Nx" << Nx << "_Ny" << Ny << ".txt";
//
//    // Data after the first pass. It must be the same in any subsequent pass.
//    Grid<point_array_t,2> data_0(GridSize);
//
//    // Testing passes.
//    for (int test = 0; test < 1000; ++test) {
//        Grid<point_array_t,2> data(GridSize);
//
//#if 1
//        // Read the sensor file sequentially.
//        {
//            FileIOManager & manager = FileIOManager::getInstance();
//            Entry e = manager.createEntry(filename.str(), Mode::Text);
//            auto in = manager.openInputStream(e);
//            for (int y = 0; y < GridSize.y; ++y) {
//            for (int x = 0; x < GridSize.x; ++x) {
//                point2d_t pt;
//                in.atomic([&](auto & file) { file >> pt.x >> pt.y; });
//                if (in) {
//                    int cell_x = pt.x / SUBDOMAIN_X;
//                    int cell_y = pt.y / SUBDOMAIN_Y;
//                    pt.x %= SUBDOMAIN_X;
//                    pt.y %= SUBDOMAIN_Y;
//                    assert_true((0 <= cell_x) && (cell_x < GridSize.x));
//                    assert_true((0 <= cell_y) && (cell_y < GridSize.y));
//                    data[{cell_x, cell_y}].push_back(pt);
//                }
//            }}
//            manager.close(in);
//
//            data.pforEach([](point_array_t & arr) {
//                if (arr.size() < arr.capacity()) {
//                    point_array_t tmp = arr;
//                    arr.swap(tmp);
//                }
//
//                std::sort(arr.begin(), arr.end(),
//                    [](const point2d_t & a, const point2d_t & b) {
//                        return ((a.x < b.x) || ((a.x == b.x) && (a.y < b.y)));
//                    }
//                );
//            });
//        }
//#else
//        // Read the sensor file in parallel.
//        {
//            FileIOManager & manager = FileIOManager::getInstance();
//            Entry e = manager.createEntry(filename.str(), Mode::Text);
//            auto input_stream = manager.openInputStream(e);
//            pfor(Origin, GridSize, [&](const point2d_t & idx) {
//                auto in = manager.getInputStream(e);
//                while (in) {
//                    point2d_t pt;
//                    in.atomic([&](auto & file) { file >> pt.x >> pt.y; });
//                    if (in) {
//                        int cell_x = pt.x / SUBDOMAIN_X;
//                        int cell_y = pt.y / SUBDOMAIN_Y;
//                        if (idx == point2d_t(cell_x, cell_y)) {
//                            pt.x %= SUBDOMAIN_X;
//                            pt.y %= SUBDOMAIN_Y;
//                            data[idx].push_back(pt);
//                        }
//                    }
//                }
//            });
//            manager.close(input_stream);
//        }
//#endif
//
//        // First time we just copy, then compare.
//        pfor(Origin, GridSize, [&](const point2d_t & idx) {
//            auto & a = data_0[idx];
//            auto & b = data  [idx];
//            if (test == 0) {
//                a = b;
//            } else {
//                assert_true(a.size() == b.size());
//                assert_true(std::equal(a.begin(), a.end(), b.begin()));
//            }
//        });
//        std::cout << "+";
//    }
//    std::cout << std::endl;
//}



//    // Make the data structure more compact.
//    sensors.forEach([](point_array_t & arr) {
//        if (arr.size() < arr.capacity()) {
//            point_array_t tmp = arr;
//            arr.swap(tmp);
//        }
//    });
