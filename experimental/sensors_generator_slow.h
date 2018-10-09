//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Class generates pseudo-random sensor positions inside the domain.
//=============================================================================
class SensorsGenerator
{
public:
    typedef std::vector<double> double_array_t;

//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
SensorsGenerator()
{
}

//-----------------------------------------------------------------------------
// Destructor.
//-----------------------------------------------------------------------------
virtual ~SensorsGenerator()
{
}

//-----------------------------------------------------------------------------
// Function generates pseudo-randomly distributed sensor locations.
// N O T E, the sensor points are stored unordered and their coordinates are
// global, i.e. defined with respected to the whole domain's coordinate system.
//-----------------------------------------------------------------------------
virtual void Generate(const Configuration    & conf,
                      ::std::vector<Point2D> & sensors) const
{
    MY_TIME_IT("Running sensors generator ...")

    // Function scales a coordinate from [0..1] range to specified size.
    auto ScaleCoord = [](double v, int size) -> int {
        int i = static_cast<int>(std::floor(v * size));
        i = std::min(std::max(i, 0), size - 1);
        return i;
    };

    // Global (whole domain) sizes.
    const int Nx = conf.asInt("subdomain_x") * conf.asInt("num_subdomains_x");
    const int Ny = conf.asInt("subdomain_y") * conf.asInt("num_subdomains_y");

    const double fraction =
            std::min(std::max(conf.asDouble("sensor_fraction"), 0.001), 0.75);

    // Generate pseudo-random sensor locations.
    const size_t Nsensors = std::max(
        static_cast<size_t>(std::floor(fraction * Nx * Ny + 0.5)), size_t(1));
    MY_LOG(INFO) << "fraction of sensor points = " << fraction
                 << ", #observations = " << Nsensors;

    double_array_t x(Nsensors), y(Nsensors);
    InitialGuess(conf, x, y, point2d_t(0,0));
    OptimizePointLocations(x, y);

    // Save (scaled) sensor locations: x[] and y[] lie inside the range [0..1].
    sensors.resize(Nsensors);
    for (size_t k = 0; k < Nsensors; ++k) {
        sensors[k] = point2d_t(ScaleCoord(x[k], Nx),
                               ScaleCoord(y[k], Ny));
    }
}

private:
//-----------------------------------------------------------------------------
// Function evaluates objective function and its gradient.
//-----------------------------------------------------------------------------
virtual void EvaluateObjective(double & J, double_array_t & gradJ,
                                     const double_array_t & x,
                                     const double_array_t & y) const
{
    const size_t N = x.size();      // short-hand alias
    const size_t NN = N * N;
    assert_true(x.size() == y.size());

    const double EPS = std::sqrt(std::numeric_limits<double>::epsilon());

    J = 0.0;
    gradJ.resize(2*N);
    std::fill(gradJ.begin(), gradJ.end(), 0.0);

    for (size_t i = 0; i < N; ++i) {
        // Reciprocal distances to subdomain borders.
        const double r_x1 = 1.0 / (std::pow(      x[i],2) + EPS);
        const double r_x2 = 1.0 / (std::pow(1.0 - x[i],2) + EPS);
        const double r_y1 = 1.0 / (std::pow(      y[i],2) + EPS);
        const double r_y2 = 1.0 / (std::pow(1.0 - y[i],2) + EPS);

        J += (r_x1 + r_x2 +
              r_y1 + r_y2);

        double gx = 0.0, gy = 0.0;
        for (size_t j = 0; j < N; ++j) {
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

    J /= (double)NN;
    std::transform(gradJ.begin(), gradJ.end(), gradJ.begin(),
                                    [=](double x){ return x/(double)NN; });
}

//-----------------------------------------------------------------------------
// Generates initial space distribution of sensor points.
//-----------------------------------------------------------------------------
virtual void InitialGuess(const Configuration & conf,
                          double_array_t      & x,
                          double_array_t      & y,
                          const point2d_t     & pos) const
{
    const int ny = conf.asInt("num_subdomains_y");
    std::mt19937_64 gen(RandomSeed() + (uint64_t)(pos.x * ny + pos.y));
    std::uniform_real_distribution<double> distrib(0.0, 1.0);
    std::generate(x.begin(), x.end(), [&](){ return distrib(gen); });
    std::generate(y.begin(), y.end(), [&](){ return distrib(gen); });
}

//-----------------------------------------------------------------------------
// Minimizes the objective function by gradient descent.
//-----------------------------------------------------------------------------
virtual void OptimizePointLocations(double_array_t & x,
                                    double_array_t & y) const
{
    const size_t N = x.size();          // short-hand alias
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
        for (size_t k = 0; (k < N) && is_inside; ++k) {
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

};  // class SensorGenerator

}   // namespace amdados

