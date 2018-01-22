//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
//             Fearghal O'Donncha, feardonn@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include "allscale/api/user/data/adaptive_grid.h"
#include "allscale/api/user/algorithm/pfor.h"
#include "allscale/api/user/algorithm/stencil.h"
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



//=============================================================================
// Variables associated with a sub-domain.
//=============================================================================
struct SubdomainData
{
    Matrix       field;         // sub-domain represented as a matrix
    Boundary     boundaries;    // sub-domain boundaries

    KalmanFilter Kalman;        // Kalman filter
    Matrix       B;             // inverse model matrix

    Matrix       P;             // process model covariance
    Matrix       Q;             // process noise covariance
    Matrix       H;             // observation matrix
    Matrix       R;             // observation noise covariance
    Vector       z;             // observation vector

    point2d_t    idx;           // subdomain's position on the grid
    size_t       Nt;            // number of time integration steps
    size_t       Nschwarz;      // number of Schwarz iterations
    flow_t       flow;          // current flow vector (vel_x, vel_y)

    SubdomainData()
        : field(), boundaries()
        , Kalman(), B()
        , P(), Q(), H(), R(), z()
        , idx(-1,-1), Nt(0), Nschwarz(0), flow(0.0, 0.0) {}
};

// The whole domain where instead of grid cells we place sub-domain data.
using data_domain_t = ::allscale::api::user::data::Grid<SubdomainData, 2>;

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
    assert_true(dt > TINY);
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
 * Function computes flow components given a discrete time.
 * \return a pair of flow components (flow_x, flow_y).
 */
flow_t Flow(const Configuration & conf, const size_t discrete_time)
{
    const double max_vx = conf.asDouble("flow_model_max_vx");
    const double max_vy = conf.asDouble("flow_model_max_vy");
    const double t = static_cast<double>(discrete_time) / conf.asDouble("Nt");
    return flow_t( -max_vx * std::sin(0.1 * t - M_PI),
                   -max_vy * std::sin(0.2 * t - M_PI) );
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
                     Boundary            & border,
                     const point2d_t     & idx,
                     const domain_t      & curr_domain,
                     domain_t            & next_domain,
                     const flow_t        & flow)
{
    // Origin and the number of subdomains in both directions.
    const int Ox = 0, Nx = curr_domain.size()[_X_];
    const int Oy = 0, Ny = curr_domain.size()[_Y_];
    assert_true(Nx == next_domain.size()[_X_]);
    assert_true(Ny == next_domain.size()[_Y_]);

    // Index increments and corresponding directions
    // of a neighbour subdomain (remote peer).
    Direction remote_dir[NSIDES];
    point2d_t remote_idx[NSIDES];

    // Given index and position of a subdomain boundary, we use lookup tables
    // to access the neighbour subdomain (remote peer).
    remote_dir[Left ] = Right;  remote_idx[Left ] = point2d_t{-1,0};
    remote_dir[Right] = Left;   remote_idx[Right] = point2d_t{+1,0};
    remote_dir[Down ] = Up;     remote_idx[Down ] = point2d_t{0,-1};
    remote_dir[Up   ] = Down;   remote_idx[Up   ] = point2d_t{0,+1};

    // Function updates a boundary depending on flow direction and accumulates
    // the difference between this and peer subdomain using Schwarz method.
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

            // Copy values from peer's boundary to this subdomain boundary.
            border.inflow[dir] = true;          // in-flow boundary
            border.myself = next_domain[idx].getBoundary(dir);
            border.remote = curr_domain[peer_idx].getBoundary(peer_dir);
            next_domain[idx].setBoundary(dir, border.remote);

            // Compute and accumulate aggregated difference between
            // boundary values of this subdomain before and after update.
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
 * Function is invoked for each sub-domain during the time integration.
 */
const subdomain_cell_t & SubdomainRoutine(
                            const Configuration & conf,
                            const point_array_t & sensors,
                            const Matrix        & observations,
                            const bool            external,
                            const size_t          timestamp,
                            const domain_t      & curr_state,
                            domain_t            & next_state,
                            SubdomainData       & vars,
                            AverageProfile      & diff_profile)
{
    // TODO !!!!!!!!!!!!!!!!!!!!! properly distinguish current and next states

    assert_true(&curr_state != &next_state);
    assert_true(curr_state.size() == next_state.size());

    // Index (position) of the current subdomain.
    const point2d_t idx = vars.idx;

    // Important: synchronize active layers, otherwise when we exchange data
    // at the borders of neighbour subdomains the size mismatch can happen.
    next_state[idx].setActiveLayer(curr_state[idx].getActiveLayer());

    // Get the discrete time (index of iteration) in the range [0..Nt) and
    // the index of Schwarz sub-iteration in the range [0..Nschwarz).
    assert_true(timestamp < vars.Nt * vars.Nschwarz);
    const size_t t_discrete = timestamp / vars.Nschwarz;
    const size_t sub_iter   = timestamp % vars.Nschwarz;

    // At the beginning of a regular iteration (i.e. at the first Schwarz
    // sub-iteration) do: (1) update the flow velocity, (2) get the new
    // observations; (3) obtain new noise covariance matrices; (4) compute
    // the prior estimation of the state field.
    if (sub_iter == 0) {
        if (idx == point2d_t(0,0)) MY_PROGRESS("+")

        // Compute flow velocity vector.
        vars.flow = Flow(conf, t_discrete);

        // Get the current sensor measurements.
        GetObservations(vars.z, observations, t_discrete, vars.H, sensors);

        // Covariance matrices are allowed to change over time.
        ComputeQ(conf, vars.Q);
        ComputeR(conf, vars.R);

        // Copy state field into the matrix object.
        MatrixFromAllscale(vars.field,
                           curr_state[idx].getLayer<ACTIVE_LAYER>());
        // Prior estimation.
        InverseModelMatrix(vars.B, conf, vars.boundaries, vars.flow);
        vars.Kalman.PropagateStateInverse(vars.field, vars.P,
                                              vars.B, vars.Q);
        // Put the prior estimation back to the state field.
        AllscaleFromMatrix(next_state[idx].getLayer<ACTIVE_LAYER>(),
                           vars.field);
        // Ensure boundary conditions on the outer border.
        if (external) {
            ApplyBoundaryCondition(next_state, idx);
        }
    } else {
        if (idx == point2d_t(0,0)) MY_PROGRESS(".")
    }

    // Update subdomain boundaries according to Schwarz method.
    double diff = SchwarzUpdate(conf, vars.boundaries, idx,
                                curr_state, next_state, vars.flow);

    // Ensure boundary conditions on the outer border.
    if (external) {
        ApplyBoundaryCondition(next_state, idx);
    }

    // Accumulate history of discrepancies for debugging,
    diff_profile.Accumulate(idx, diff);

    // Copy state field into more convenient matrix object.
    MatrixFromAllscale(vars.field,
                       next_state[idx].getLayer<ACTIVE_LAYER>());

    // Filtering by Kalman filter.
    vars.Kalman.SolveFilter(vars.field, vars.P, vars.H, vars.R, vars.z);

    // Put new estimation back to the state field.
    AllscaleFromMatrix(next_state[idx].getLayer<ACTIVE_LAYER>(),
                       vars.field);

    // Ensure boundary conditions on the outer border.
    if (external) {
        ApplyBoundaryCondition(next_state, idx);
    }

    return next_state[idx];
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
    using ::allscale::api::core::FileIOManager;
    using ::allscale::api::core::Entry;
    using ::allscale::api::core::Mode;

    MY_INFO("%s", "Running simulation with data assimilation")

	//Visualizer	   vis(conf.asString("output_dir"));
    AverageProfile diff_profile(conf);

    const point2d_t GridSize = GetGridSize(conf);   // size in subdomains
    const size_t    Nt = conf.asUInt("Nt");
    const size_t    Nschwarz = conf.asUInt("schwarz_num_iters");
    const size_t    Nwrite = conf.asUInt("write_num_fields");

    // Visualization: image of the sensor positions in the whole domain.
//XXX    vis.WriteImageOfSensors(H);

    data_domain_t state_vars(GridSize);         // variables of each sub-domain
    domain_t      temp_field(GridSize);         // grid if sub-domains
    domain_t      state_field(GridSize);        // grid if sub-domains

    // Set up the initial field of the physical entity distribution.
    InitialField(state_field, conf, "zero");    // "zero" or "gauss"

    // Initialize the observation and model covariance matrices.
    pfor(point2d_t(0,0), GridSize, [&](const point2d_t & idx) {
        const int Nsensors = static_cast<int>(sensors[idx].size());
        state_vars[idx].field.Resize(SUBDOMAIN_X, SUBDOMAIN_Y);
        state_vars[idx].B.Resize(SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE);
        state_vars[idx].P.Resize(SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE);
        state_vars[idx].Q.Resize(SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE);
        state_vars[idx].H.Resize(Nsensors, SUB_PROBLEM_SIZE);
        state_vars[idx].R.Resize(Nsensors, Nsensors);
        state_vars[idx].z.Resize(Nsensors);
        state_vars[idx].idx = idx;
        state_vars[idx].Nt = Nt;
        state_vars[idx].Nschwarz = Nschwarz;
        ComputeH(sensors[idx], state_vars[idx].H);
        InitialCovar(conf, state_vars[idx].P);
        state_field[idx].setActiveLayer(ACTIVE_LAYER);  // this is important
        temp_field[idx].setActiveLayer(ACTIVE_LAYER);   // this is important
    });

    // Open file manager and the output file for writing.
    std::string filename = MakeFileName(conf, "field");
    FileIOManager & file_manager = FileIOManager::getInstance();
    Entry stream_entry = file_manager.createEntry(filename, Mode::Binary);
    auto out_stream = file_manager.openOutputStream(stream_entry);

    // Time integration forward in time. We want to make Nt (normal) iterations
    // and Nschwarz Schwarz sub-iterations within each (normal) iteration.
    ::allscale::api::user::algorithm::stencil(
        state_field, Nt * Nschwarz,
        // Process the internal subdomains.
        [&](time_t t, const point2d_t & idx, const domain_t & state)
        -> const subdomain_cell_t &  // cell is not copy-constructible, so '&'
        {
            assert_true(t >= 0);
            return SubdomainRoutine(conf, sensors[idx], observations[idx],
                                    false, size_t(t), state, temp_field,
                                    state_vars[idx], diff_profile);
        },
        // Process the external subdomains.
        [&](time_t t, const point2d_t & idx, const domain_t & state)
        -> const subdomain_cell_t &  // cell is not copy-constructible, so '&'
        {
            return SubdomainRoutine(conf, sensors[idx], observations[idx],
                                    true, size_t(t), state, temp_field,
                                    state_vars[idx], diff_profile);
        },
        // Monitoring.
        ::allscale::api::user::algorithm::observer(
            // Time filter: choose few moments evenly distributed on time axis.
            [=](time_t t) {
                // Filter out the first Schwarz sub-iteration, skip the others.
                if ((t % time_t(Nschwarz)) != 0)
                    return false;
                t /= time_t(Nschwarz);
                return ((((Nwrite-1)*(t-1))/(Nt-1) != ((Nwrite-1)*t)/(Nt-1)));
            },
            // Space filter: no specific points.
            [](const point2d_t &) { return true; },
            // Append a full field to the file of simulation results.
            [&](time_t t, const point2d_t & p, const subdomain_cell_t & cell) {
                t /= time_t(Nschwarz);
                const auto & subdomain = cell.getLayer<ACTIVE_LAYER>();
                int gx = 0, gy = 0;     // global coordinates
                for (int y = 0; y < SUBDOMAIN_Y; ++y) { gy = Sub2GloY(p, y);
                for (int x = 0; x < SUBDOMAIN_X; ++x) { gx = Sub2GloX(p, x);
                    double val = subdomain[{x,y}];
                    out_stream.atomic([=](auto & file) {
                        file.write(static_cast<float>(t));
                        file.write(static_cast<float>(gx));
                        file.write(static_cast<float>(gy));
                        file.write(static_cast<float>(val));
                    });
                }}
            }
        )
    );
    file_manager.close(out_stream);
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
