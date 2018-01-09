//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
//             Fearghal O'Donncha, feardonn@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/core/io.h"
#include "allscale/utils/assert.h"
#include "../include/debugging.h"
#include "../include/geometry.h"
#include "../include/amdados_utils.h"
#include "../include/matrix.h"
#include "../include/configuration.h"
#include "../include/cholesky.h"
#include "../include/demo_visualizer.h"
#include "../include/lu.h"
#include "../include/kalman_filter.h"

namespace amdados {

using ::allscale::api::user::algorithm::pfor;
using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;
using ::allscale::api::user::data::Direction;
using ::allscale::api::user::data::Direction::Up;
using ::allscale::api::user::data::Direction::Down;
using ::allscale::api::user::data::Direction::Left;
using ::allscale::api::user::data::Direction::Right;

// amdados_utils.cpp:
point2d_t GetGridSize(const Configuration & conf);

// scenario_sensors.cpp:
void LoadSensorLocations(const Configuration   & conf,
                         Grid<point_array_t,2> & sensors);
void LoadSensorMeasurements(const Configuration         & conf,
                            const Grid<point_array_t,2> & sensors,
                            Grid<Matrix,2>              & observations);

namespace {

const int NSIDES = 4;             // number of sides any subdomain has
const int ACTIVE_LAYER = L_100m;  // should be template a parameter eventually

/**
 * Structure keeps information about 4-side boundary of a subdomain
 * (Up, Down, Left, Right).
 */
struct Boundary
{
    double_array_t myself; // temporary storage for this subdomain boundary
    double_array_t remote; // temporary storage for remote boundary values

    double rel_diff;       // relative difference across in-flow borders
    bool   inflow[NSIDES]; // true, if flow is coming in along a side
    bool   outer[NSIDES];  // true, if side belongs to domain's outer boundary
};

typedef std::pair<double,double> flow_t;    // flow components (flow_x, flow_y)

//#############################################################################
// @Utilities.
//#############################################################################

/**
 * Function computes: z = H * observations(t). Since H is a simple 0/1 matrix
 * that just picks up the observations at sensor locations, instead of
 * matrix-vector multiplication we get the observations directly.
 */
void GetObservations(Vector & z, const Matrix & observations, int timestep,
                     const Matrix & H, const point_array_t & sensors)
{
    const int n = static_cast<int>(observations.NCols());
    assert_true(z.Size() == n);
    for (int i = 0; i < n; ++i) { z(i) = observations(timestep, i); }

    (void)H;(void)sensors;
#if 1
#warning "Some extra validity test"
    Matrix subfield(SUBDOMAIN_X, SUBDOMAIN_Y);
    for (int i = 0; i < static_cast<int>(sensors.size()); ++i) {
        subfield(sensors[i].x, sensors[i].y) = observations(timestep, i);
    }
    Vector zz(n);
    MatVecMult(zz, H, subfield);    // z = H * observations(t)
    assert_true(std::equal(zz.begin(), zz.end(),
                            z.begin(),  z.end()));
#endif
}

/**
 * Function writes the entire state field into a file.
 */
void WriteFullField(const Configuration & conf,
                    const domain_t & state, const int timestamp)
{
    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;

    const int SUFFIX_LEN = 63;
    char suffix[SUFFIX_LEN + 1];
    snprintf(suffix, SUFFIX_LEN, "time%05d.txt", timestamp);

    // Open file manager and the output file for writing.
    std::string filename = MakeFileName(conf, "field", suffix);
    FileIOManager & manager = FileIOManager::getInstance();
    Entry e = manager.createEntry(filename, Mode::Text);
    auto out = manager.openOutputStream(e);

    pfor(point2d_t(0,0), state.size(), [&](const point2d_t & idx) {
        auto & f = state[idx].getLayer<ACTIVE_LAYER>();
        for (int y = 0; y < SUBDOMAIN_Y; ++y) { int j = Sub2GloY(idx, y);
        for (int x = 0; x < SUBDOMAIN_X; ++x) { int i = Sub2GloX(idx, x);
            double v = f[{x,y}];
            out.atomic([=](auto & file) {
                file << i << " " << j << " " << v << "\n";
            });
        }}
    });
    manager.close(out);
}

//#############################################################################
// @Initialization.
//#############################################################################

/**
 * Function initializes dependent parameters given the primary ones
 * specified by user.
 */
void InitDependentParams(Configuration & conf)
{
    // Check some global constants.
    static_assert((0 <= Up   ) && (Up    < NSIDES), "");
    static_assert((0 <= Down ) && (Down  < NSIDES), "");
    static_assert((0 <= Left ) && (Left  < NSIDES), "");
    static_assert((0 <= Right) && (Right < NSIDES), "");
    assert_true((SUBDOMAIN_X >= 3) && (SUBDOMAIN_Y >= 3))
                                    << "subdomain must be at least 3x3";

    // Ensure integer values for certain parameters.
    assert_true(conf.IsInteger("num_subdomains_x"));
    assert_true(conf.IsInteger("num_subdomains_y"));
    assert_true(conf.IsInteger("subdomain_x"));
    assert_true(conf.IsInteger("subdomain_y"));
    assert_true(conf.IsInteger("integration_nsteps"));

    // Check the subdomain size: hard-coded value must match the parameter.
    assert_true(conf.asInt("subdomain_x") == SUBDOMAIN_X)
                        << "subdomain_x mismatch" << std::endl;
    assert_true(conf.asInt("subdomain_y") == SUBDOMAIN_Y)
                        << "subdomain_y mismatch" << std::endl;

    const point2d_t grid_size = GetGridSize(conf);
    const int nx = grid_size.x * SUBDOMAIN_X;       // global X-size
    const int ny = grid_size.y * SUBDOMAIN_Y;       // global Y-size

    const double D = conf.asDouble("diffusion_coef");
    assert_true(D > 0.0);

    conf.SetInt("global_problem_size", nx * ny);
    const double dx = conf.asDouble("domain_size_x") / (nx - 1);
    const double dy = conf.asDouble("domain_size_y") / (ny - 1);
    assert_true((dx > 0) && (dy > 0));
    conf.SetDouble("dx", dx);
    conf.SetDouble("dy", dy);

    // Deduce the optimal time step from the stability criteria.
    const double dt_base = conf.asDouble("integration_period") /
                           conf.asDouble("integration_nsteps");
    const double max_vx = conf.asDouble("flow_model_max_vx");
    const double max_vy = conf.asDouble("flow_model_max_vy");
    const double dt = std::min(dt_base,
                        std::min( std::min(dx*dx, dy*dy)/(2.0*D + TINY),
                                  1.0/(std::fabs(max_vx)/dx +
                                       std::fabs(max_vy)/dy + TINY) ));
    assert_true(dt > 0);
    conf.SetDouble("dt", dt);
    conf.SetInt("Nt", static_cast<int>(
                        std::ceil(conf.asDouble("integration_period") / dt)));
}

/**
 * Function applies Dirichlet zero boundary condition at the outer border
 * of the domain.
 */
void ApplyBoundaryCondition(domain_t & state, const point2d_t & idx)
{
    const int Ox = 0, Nx = state.size()[_X_], Sx = SUBDOMAIN_X;
    const int Oy = 0, Ny = state.size()[_Y_], Sy = SUBDOMAIN_Y;

    auto & f = state[idx].getLayer<ACTIVE_LAYER>();

    // Set the leftmost and rightmost.
    if (idx[_X_] == Ox)   { for (int y = 0; y < Sy; ++y) f[{   0,y}] = 0.0; }
    if (idx[_X_] == Nx-1) { for (int y = 0; y < Sy; ++y) f[{Sx-1,y}] = 0.0; }

    // Set the bottommost and topmost.
    if (idx[_Y_] == Oy)   { for (int x = 0; x < Sx; ++x) f[{x,   0}] = 0.0; }
    if (idx[_Y_] == Ny-1) { for (int x = 0; x < Sx; ++x) f[{x,Sy-1}] = 0.0; }
}

/**
 * Function initializes and fills up initial field of density distribution.
 * It is either all zeros or a spike at some point and zeros elsewhere.
 * Note, the spike is not very sharp to make the field differentiable.
 * \param  state       multi-layered structure that keeps density fields of
 *                     all the subdomains.
 * \param  conf        configuration parameters.
 * \param  field_type  initial field type is one of: {"zero", "gauss"}.
 */
void InitialField(domain_t            & state,
                  const Configuration & conf,
                  const char          * field_type)
{
    if (std::strcmp(field_type, "zero") == 0)
    {
        pfor(point2d_t(0,0), state.size(), [&](const point2d_t & idx) {
            state[idx].setActiveLayer(ACTIVE_LAYER);
            state[idx].forAllActiveNodes([](double & value) { value = 0.0; });
            ApplyBoundaryCondition(state, idx);
        });
        MY_INFO("%s", "Initial state field: all zeros")
    }
    else if (std::strcmp(field_type, "gauss") == 0)
    {
        // Global coordinates of the density spot center.
        const int cx = Round(conf.asDouble("spot_x") / conf.asDouble("dx"));
        const int cy = Round(conf.asDouble("spot_y") / conf.asDouble("dy"));
        const point2d_t glo = GlobalSize(state);
        assert_true((0 <= cx) && (cx < glo.x) &&
                    (0 <= cy) && (cy < glo.y))
                    << "high-concentration spot is not inside the domain"
                    << std::endl;

        // Parameters of the global 2D Gaussian model of the spike.
        const double sigma = 1.0;       // in logical units (point indices)
        const double a = conf.asDouble("spot_density") /
                            (std::pow(sigma,2) * 2.0 * M_PI);
        const double b = 1.0 / (2.0 * std::pow(sigma,2));

        // For all the subdomains initialize the spike distribution by parts.
        pfor(point2d_t(0,0), state.size(), [&](const point2d_t & idx) {
            // Set active layer to zero at the beginning.
            state[idx].setActiveLayer(ACTIVE_LAYER);
            state[idx].forAllActiveNodes([](double & value) { value = 0.0; });

            auto & subfield = state[idx].getLayer<ACTIVE_LAYER>();
            for (int x = 0; x < SUBDOMAIN_X; ++x) {
            for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                int dx = Sub2GloX(idx,x) - cx;
                int dy = Sub2GloY(idx,y) - cy;
                if ((std::abs(dx) <= 4*sigma) && (std::abs(dy) <= 4*sigma)) {
                    subfield[{x,y}] += a * std::exp(-b * (dx*dx + dy*dy));
                }
            }}
            ApplyBoundaryCondition(state, idx);
        });
        MY_INFO("%s", "Initial state field: Gaussian peaked at some point")
    } else {
        assert_true(0) << "InitialField(): unknown field type" << std::endl;
    }
}

//#############################################################################
// @Kalman filter stuff.
//#############################################################################

/**
 * Function computes the initial process model covariance matrix P
 * based on exponential distance.
 * TODO: better initial correlation radius
 */
void InitialCovar(const Configuration & conf, Matrix & P)
{
    // Here we express correlation distances in logical coordinates
    // of nodal points.
    const double variance = conf.asDouble("model_ini_var");
    const double covar_radius = conf.asDouble("model_ini_covar_radius");
    const double sx = std::max(covar_radius / conf.asDouble("dx"), 1.0);
    const double sy = std::max(covar_radius / conf.asDouble("dy"), 1.0);
    const int Rx = Round(std::ceil(4.0 * sx));
    const int Ry = Round(std::ceil(4.0 * sy));
    const int Nx = SUBDOMAIN_X;                 // short-hand aliases
    const int Ny = SUBDOMAIN_Y;

    assert_true(P.IsSquare());
    Fill(P, 0.0);
    for (int u = 0; u < SUBDOMAIN_X; ++u) {
    for (int v = 0; v < SUBDOMAIN_Y; ++v) {
        int i = Sub2Ind(u,v);
        double dx = 0.0, dy = 0.0;
        for (int x = u-Rx; x <= u+Rx; ++x) { if ((0 <= x) && (x < Nx)) {
        for (int y = v-Ry; y <= v+Ry; ++y) { if ((0 <= y) && (y < Ny)) {
            int j = Sub2Ind(x,y);
            if (i <= j) {
                dx = (u - x) / sx;
                dy = (v - y) / sy;
                P(i,j) = P(j,i) = variance * std::exp(-0.5 * (dx*dx + dy*dy));
            }
        }}}}
    }}
}

/**
 * Function computes the process model noise covariance matrix.
 */
void ComputeQ(const Configuration & conf, Matrix & Q)
{
    std::uniform_real_distribution<double> distrib;
    std::mt19937_64 gen(RandomSeed());
    const double model_noise_Q = conf.asDouble("model_noise_Q");

    assert_true(Q.IsSquare());
    Fill(Q, 0.0);
    for (int k = 0; k < Q.NRows(); ++k) {
        Q(k,k) = 1.0 + model_noise_Q * distrib(gen);    // always >= 1
    }
}

/**
 * Function computes the measurement noise covariance matrix.
 */
void ComputeR(const Configuration & conf, Matrix & R)
{
    std::uniform_real_distribution<double> distrib;
    std::mt19937_64 gen(RandomSeed());
    const double model_noise_R = conf.asDouble("model_noise_R");

    assert_true(R.IsSquare());
    Fill(R, 0.0);
    for (int k = 0; k < R.NRows(); ++k) {
        R(k,k) = 1.0 + model_noise_R * distrib(gen);    // always >= 1
    }
}

/**
 * Function initializes the observation matrix H of size:
 * (number of sensors in subdomain) x (number of nodal points in subdomain).
 */
void ComputeH(const point_array_t & sensors, Matrix & H)
{
    assert_true(H.NRows() == static_cast<int>(sensors.size()));
    assert_true(H.NCols() == SUB_PROBLEM_SIZE);
    Fill(H, 0.0);
    for (size_t k = 0; k < sensors.size(); ++k) {
        int x = sensors[k].x;
        int y = sensors[k].y;
        H(static_cast<int>(k), Sub2Ind(x,y)) = 1.0;
    }
}

//#############################################################################
// @Advection-diffusion PDE stuff.
//#############################################################################

/**
 * Function computes flow components given a time.
 * \return a pair of flow components (flow_x, flow_y).
 */
flow_t Flow(const Configuration & conf, const double physical_time)
{
    const double max_vx = conf.asDouble("flow_model_max_vx");
    const double max_vy = conf.asDouble("flow_model_max_vy");
    const double T = conf.asDouble("integration_period");
    return flow_t( -max_vx * std::sin(0.1 * physical_time / T - M_PI),
                   -max_vy * std::sin(0.2 * physical_time / T - M_PI) );
}

/**
 * Function initializes inverse matrix of implicit Euler time-integrator:
 * B * x_{t+1} = x_{t}, where B = A^{-1} is the matrix returned by this
 * function. The matrix must be inverted while iterating forward in time:
 * x_{t+1} = A * x_{t}.
 * Note, the matrix we generate here is acting on a subdomain.
 * Note, the model matrix B is supposed to be a sparse one. For now, since we
 * do not have a fast utility for sparse matrix inversion, we define B as
 * a dense one with many zeros.
 */
void InverseModelMatrix(Matrix & B, const Configuration & conf,
                        const Boundary & boundary, const flow_t & flow)
{
    const int Nx = SUBDOMAIN_X;    // short-hand aliases
    const int Ny = SUBDOMAIN_Y;

    const double D  = conf.asDouble("diffusion_coef");
    const double dx = conf.asDouble("dx");
    const double dy = conf.asDouble("dy");
    const double dt = conf.asDouble("dt");

    const double rho_x = D * dt / std::pow(dx,2);
    const double rho_y = D * dt / std::pow(dy,2);

    const double v0x = 2.0 * dx / dt;
    const double v0y = 2.0 * dy / dt;

    const double vx = flow.first  / v0x;
    const double vy = flow.second / v0y;

    // This operation can be avoided in case of sparse matrix.
    Fill(B, 0.0);

    // Process the internal and boundary subdomain points separately.
    // At each border we treat outside points with finite difference scheme.
    // If there is no incoming flow, we replace missed point values (outside
    // a subdomain) by the values inside subdomain in such way that the field
    // derivative along the border normal is zero.
    for (int x = 0; x < Nx; ++x) {
    for (int y = 0; y < Ny; ++y) {
        int i = Sub2Ind(x,y);

        if ((x == 0) || (x+1 == Nx) || (y == 0) || (y+1 == Ny)) { // border
            B(i,i) += 1 + 2*(rho_x + rho_y);

            if (x == 0) {
                if (boundary.inflow[Left]) {
                    B(i,Sub2Ind(x  ,y)) += - 2*vx - rho_x;
                    B(i,Sub2Ind(x+1,y)) += + 2*vx - rho_x;
                } else {
                    B(i,Sub2Ind(x+1,y)) += - vx - rho_x;
                    B(i,Sub2Ind(x+1,y)) += + vx - rho_x;
                }
            } else if (x == Nx-1) {
                if (boundary.inflow[Right]) {
                    B(i,Sub2Ind(x-1,y)) += - 2*vx - rho_x;
                    B(i,Sub2Ind(x  ,y)) += + 2*vx - rho_x;
                } else {
                    B(i,Sub2Ind(x-1,y)) += - vx - rho_x;
                    B(i,Sub2Ind(x-1,y)) += + vx - rho_x;
                }
            } else {
                B(i,Sub2Ind(x-1,y)) += - vx - rho_x;
                B(i,Sub2Ind(x+1,y)) += + vx - rho_x;
            }

            if (y == 0) {
                if (boundary.inflow[Down]) {
                    B(i,Sub2Ind(x,y  )) += - 2*vy - rho_y;
                    B(i,Sub2Ind(x,y+1)) += + 2*vy - rho_y;
                } else {
                    B(i,Sub2Ind(x,y+1)) += - vy - rho_y;
                    B(i,Sub2Ind(x,y+1)) += + vy - rho_y;
                }
            } else if (y == Ny-1) {
                if (boundary.inflow[Up]) {
                    B(i,Sub2Ind(x,y-1)) += - 2*vy - rho_y;
                    B(i,Sub2Ind(x,y  )) += + 2*vy - rho_y;
                } else {
                    B(i,Sub2Ind(x,y-1)) += - vy - rho_y;
                    B(i,Sub2Ind(x,y-1)) += + vy - rho_y;
                }
            } else {
                B(i,Sub2Ind(x,y-1)) += - vy - rho_y;
                B(i,Sub2Ind(x,y+1)) += + vy - rho_y;
            }
        } else {                            // internal point
            B(i,i) = 1 + 2*(rho_x + rho_y);
            B(i,Sub2Ind(x-1,y)) = - vx - rho_x;
            B(i,Sub2Ind(x+1,y)) = + vx - rho_x;
            B(i,Sub2Ind(x,y-1)) = - vy - rho_y;
            B(i,Sub2Ind(x,y+1)) = + vy - rho_y;
        }
    }}
}

/**
 * Function implements the idea of Schwarz method where the boundary values
 * of subdomain are get updated depending on flow direction.
 */
double SchwarzUpdate(const Configuration &,
                     Boundary & border, const point2d_t & idx,
                     domain_t & domain, const flow_t & flow)
{
    // Origin and the number of subdomains in both directions.
    const int Ox = 0, Nx = domain.size()[_X_];
    const int Oy = 0, Ny = domain.size()[_Y_];

    // Index increments and corresponding directions
    // of a neighbour (remote) subdomain.
    Direction remote_dir[NSIDES];
    point2d_t remote_idx[NSIDES];

    // Given index and position of a subdomain boundary, we use lookup tables
    // to access the neighbour (remote) subdomain.
    remote_dir[Left ] = Right;  remote_idx[Left ] = point2d_t{-1,0};
    remote_dir[Right] = Left;   remote_idx[Right] = point2d_t{+1,0};
    remote_dir[Down ] = Up;     remote_idx[Down ] = point2d_t{0,-1};
    remote_dir[Up   ] = Down;   remote_idx[Up   ] = point2d_t{0,+1};

    // Function updates a boundary depending on flow direction and accumulates
    // the difference between this and neighbour (remote) subdomain according
    // to Schwarz method.
    auto UpdBoundary = [&](Direction dir, bool is_outer, const point2d_t & idx,
                           double & numer_sum, double & denom_sum)
    {
        border.outer[dir] = is_outer;
        border.inflow[dir] = false;
        if (is_outer) return;       // no incoming flow on outer border

        // Make a normal vector pointing outside the subdomain.
        int normal_x = (dir == Right) ? +1 : ((dir == Left) ? -1 : 0);
        int normal_y = (dir ==    Up) ? +1 : ((dir == Down) ? -1 : 0);

        // Update the boundary points if flow enters the subdomain.
        if (normal_x * flow.first + normal_y * flow.second < 0) {
            Direction peer_dir = remote_dir[dir];
            point2d_t peer_idx = remote_idx[dir] + idx;
            assert_true((0 <= peer_idx.x) && (peer_idx.x < Nx));
            assert_true((0 <= peer_idx.y) && (peer_idx.y < Ny));

            // Copy values from peer's boundary to the current boundary.
            border.inflow[dir] = true;                      // in-flow boundary
            border.myself = domain[idx].getBoundary(dir);   // save old state
            border.remote = domain[peer_idx].getBoundary(peer_dir);
            domain[idx].setBoundary(dir, border.remote);

            // Compute and accumulate aggregated difference between
            // subdomains' borders.
            double rsum = 0.0, msum = 0.0, diff = 0.0;
            assert_true(border.myself.size() == border.remote.size());
            for (size_t k = 0, n = border.remote.size(); k < n; ++k) {
                diff += std::fabs(border.remote[k] - border.myself[k]);
                rsum += std::fabs(border.remote[k]);
                msum += std::fabs(border.myself[k]);
            }
            numer_sum += diff;
            denom_sum += std::max(rsum, msum);
        }
    };

    // Update subdomain boundaries.
    double numer_sum = 0.0, denom_sum = 0.0;

    UpdBoundary(Left,  (idx[_X_] == Ox),     idx, numer_sum, denom_sum);
    UpdBoundary(Right, (idx[_X_] == Nx - 1), idx, numer_sum, denom_sum);
    UpdBoundary(Down,  (idx[_Y_] == Oy),     idx, numer_sum, denom_sum);
    UpdBoundary(Up,    (idx[_Y_] == Ny - 1), idx, numer_sum, denom_sum);

    border.rel_diff = numer_sum / std::max(denom_sum, TINY);
    return border.rel_diff;
}

/**
 * Using model matrix A, the function integrates advection-diffusion equation
 * forward in time inside individual subdomains and records all the solutions
 * as the state fields. Schwarz iterations are applied in order to make the
 * global solution seamless along subdomain boundaries. On top of that, the
 * Kalman filters (separate filter in each subdomain) drive the solution
 * towards the observations at sensor locations (data assimilation).
 */
void RunDataAssimilation(const Configuration         & conf,
                         const Grid<point_array_t,2> & sensors,
                         const Grid<Matrix,2>        & observations)
{
    MY_INFO("%s", "Running simulation with data assimilation")

	Visualizer	   vis(conf.asString("output_dir"));
    AverageProfile diff_profile(conf);

    const point2d_t GridSize = GetGridSize(conf);   // size in subdomains
    const int       Nschwarz = conf.asInt("schwarz_num_iters");

    domain_t         state(GridSize);       // grid of state fields
    Grid<Matrix,2>   field(GridSize);       // same state field but as matrices
    Grid<Boundary,2> boundaries(GridSize);  // grid of boundaries

    Grid<KalmanFilter,2> Kalman(GridSize);  // grid of Kalman filters
    Grid<Matrix,2>       B(GridSize);       // grid of inverse model matrices

    Grid<Matrix,2> P(GridSize);  // grid of process model covariances
    Grid<Matrix,2> Q(GridSize);  // grid of process noise covariances
    Grid<Matrix,2> H(GridSize);  // grid of observation matrices
    Grid<Matrix,2> R(GridSize);  // grid of observation noise covariances
    Grid<Vector,2> z(GridSize);  // grid of observation vectors

    // Set up the initial field of the physical entity distribution.
    InitialField(state, conf, "zero");      // "zero" or "gauss"

    // Initialize the observation and model covariance matrices.
    pfor(point2d_t(0,0), GridSize, [&](const point2d_t & idx) {
        const int Nsensors = static_cast<int>(sensors[idx].size());
        field[idx].Resize(SUBDOMAIN_X, SUBDOMAIN_Y);
        B[idx].Resize(SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE);
        P[idx].Resize(SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE);
        Q[idx].Resize(SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE);
        H[idx].Resize(Nsensors, SUB_PROBLEM_SIZE);
        R[idx].Resize(Nsensors, Nsensors);
        z[idx].Resize(Nsensors);
        ComputeH(sensors[idx], H[idx]);
        InitialCovar(conf, P[idx]);
    });

    // Visualization: image of the sensor positions in the whole domain.
    vis.WriteImageOfSensors(H);

    // Time integration forward in time.
    for (int t = 0, Nt = conf.asInt("Nt"); t < Nt; ++t) {
        MY_PROGRESS('+')

        const double physical_time = t * conf.asDouble("dt");
        const flow_t flow = Flow(conf, physical_time);

        pfor(point2d_t(0,0), GridSize, [&](const point2d_t & idx) {
            // Get the current sensor measurements.
            GetObservations(z[idx], observations[idx], t,
                            H[idx], sensors[idx]);

            // Covariance matrices are allowed to change over time.
            ComputeQ(conf, Q[idx]);
            ComputeR(conf, R[idx]);

            // Fixed number of Schwarz (sub)iterations.
            for (int iter_no = 0; iter_no < Nschwarz; ++iter_no) {
                if (idx == point2d_t(0,0)) MY_PROGRESS('.')

                // On the first Schwarz iteration we computer prior
                // estimation of the state field.
                if (iter_no == 0) {
                    // Copy state field into more convenient matrix object.
                    MatrixFromAllscale(field[idx],
                                       state[idx].getLayer<ACTIVE_LAYER>());
                    // Prior estimation.
                    InverseModelMatrix(B[idx], conf, boundaries[idx], flow);
                    Kalman[idx].PropagateStateInverse(field[idx], P[idx],
                                                          B[idx], Q[idx]);
                    // Put new estimation back to the state field.
                    AllscaleFromMatrix(state[idx].getLayer<ACTIVE_LAYER>(),
                                       field[idx]);
                    // Ensure boundary conditions on the outer border.
                    ApplyBoundaryCondition(state, idx);
                }

                // Update subdomain boundaries according to Schwarz method.
                double diff = SchwarzUpdate(conf, boundaries[idx],
                                                             idx, state, flow);
                // Ensure boundary conditions on the outer border.
                ApplyBoundaryCondition(state, idx);
                // Accumulate history of discrepancies for debugging,
                diff_profile.Accumulate(idx, diff);
                // Copy state field into more convenient matrix object.
                MatrixFromAllscale(field[idx],
                                   state[idx].getLayer<ACTIVE_LAYER>());
                // Filtering by Kalman filter.
                Kalman[idx].SolveFilter(field[idx], P[idx],
                                            H[idx], R[idx], z[idx]);
                // Put new estimation back to the state field.
                AllscaleFromMatrix(state[idx].getLayer<ACTIVE_LAYER>(),
                                   field[idx]);
                // Ensure boundary conditions on the outer border.
                ApplyBoundaryCondition(state, idx);
            }
        });

        // Write a few full fields as the simulation results.
        if ((t == 0) || (t+1 == Nt) || ((10*t)/Nt != (10*(t+1))/Nt)) {
            WriteFullField(conf, state, t);
        }

        // Visualization: print and/or write the image of the state field.
        vis.ShowImage(field, "field", t, true, true);
    }
    diff_profile.PrintProfile(conf, "schwarz_diff");
    MY_INFO("%s", "\n\n")
}

} // anonymous namespace

/**
 * The main function of this application runs simulation with data
 * assimilation using Schwarz method to handle domain subdivision.
 */
void ScenarioSimulation(const std::string & config_file)
{
    MY_INFO("%s", "***** Amdados2D application *****")
    MY_INFO("%s", "AMDADOS_DEBUGGING is enabled")

    // Read application parameters from configuration file,
    // prepare the output directory.
    Configuration conf;
    conf.ReadConfigFile(config_file);
    InitDependentParams(conf);
    conf.PrintParameters();
    //CleanOutputDir(conf.asString("output_dir"));

    // Load sensor data obtained from Python code.
    Grid<point_array_t,2> sensors(GetGridSize(conf));
    Grid<Matrix,2> observations(GetGridSize(conf));
    LoadSensorLocations(conf, sensors);
    LoadSensorMeasurements(conf, sensors, observations);

    // Run the simulation with data assimilation. Important: by this time
    // some parameters had been initialized in InitDependentParams(..), so
    // we can safely proceed to the main part of the simulation algorithm.
    RunDataAssimilation(conf, sensors, observations);

    MakeVideo(conf, "field");
}

} // namespace amdados
