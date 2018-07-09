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
#include "../include/lu.h"
#include "../include/kalman_filter.h"
#include "../include/demo_average_profile.h"

// There are two methods to implement multi-scaling.
#define MY_MULTISCALE_METHOD 2
#if (MY_MULTISCALE_METHOD != 1) && (MY_MULTISCALE_METHOD != 2)
#error wrong MY_MULTISCALE_METHOD macro value
#endif

namespace amdados {

using ::allscale::api::user::algorithm::pfor;
using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;
using ::allscale::api::user::data::Direction;

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

	friend std::ostream& operator<<(std::ostream& out, const Boundary& b) {
		out << "Boundary: [ ";
		for(const auto& e : b.myself) { out << " " << e; }
		out << ", ";
		for(const auto& e : b.remote) { out << " " << e; }
		out << ", " << b.rel_diff;
		for(int i = 0; i < NSIDES; ++i) { out << " " << b.inflow[i]; }
		out << ", ";
		for(int i = 0; i < NSIDES; ++i) { out << " " << b.outer[i]; }
		out << " ]" << std::endl;
		return out;
	}

};

typedef std::pair<double,double> flow_t;    // flow components (flow_x, flow_y)



//=============================================================================
// Variables and data associated with a sub-domain.
//=============================================================================
struct SubdomainContext
{
    Matrix        field;        // sub-domain represented as a matrix
    Boundary      boundaries;   // sub-domain boundaries

    KalmanFilter  Kalman;       // Kalman filter
    Matrix        B;            // inverse model matrix

    Matrix        P;            // process model covariance
    Matrix        Q;            // process noise covariance
    Matrix        H;            // observation matrix
    Matrix        R;            // observation noise covariance
    Vector        z;            // observation vector

    point_array_t   sensors;    // sensor locations at the finest resolution
    LUdecomposition LU;         // used for state propagation without sensors
    Matrix          tmp_field;  // used for state propagation without sensors

    flow_t        flow;         // current flow vector (vel_x, vel_y)

    SubdomainContext()
        : field(), boundaries()
        , Kalman(), B()
        , P(), Q(), H(), R(), z()
        , sensors(), LU(), tmp_field()
        , flow(0.0, 0.0)
    {}

	friend std::ostream& operator<<(std::ostream& out, const SubdomainContext& ctx) {
		out << "SubdomainContext: [ ";
		out << ctx.field << ", ";
		out << ctx.boundaries << ", ";
		out << ctx.Kalman << ", ";
		out << ctx.B << ", ";
		out << ctx.P << ", ";
		out << ctx.Q << ", ";
		out << ctx.H << ", ";
		out << ctx.R << ", ";
		out << ctx.z;
		for(const auto& e : ctx.sensors) { out << ", " << e; }
		out << ctx.LU << ", ";
		out << ctx.tmp_field << ", ";
		out << ctx.flow.first << ", ";
		out << ctx.flow.second;
		out << " ]" << std::endl;
		return out;
	}

};

// The whole domain where instead of grid cells we place sub-domain data.
using context_domain_t = ::allscale::api::user::data::Grid<SubdomainContext,2>;

/**
 * Function converts 2D point to a flat 1D index for extended (!!!) subdomain.
 * Index layout ('y' is faster than 'x') matches to row-major Matrix class.
 */
inline index_t sub2ind(index_t x, index_t y, const size2d_t & layer_size)
{
#ifndef NDEBUG
    if (!((static_cast<size_t>(x) < static_cast<size_t>(layer_size.x + 2)) &&
          (static_cast<size_t>(y) < static_cast<size_t>(layer_size.y + 2))))
        assert_true(0);
#endif
    return (x * (layer_size.y + 2) + y);
}

/**
 * Function returns "true" if the active layer and the field specified
 * have the matching sizes.
 */
inline bool CheckSizes(const subdomain_t & cell, const Matrix & field)
{
    // Mind the extended subdomain: one extra point layer on either side.
    const size2d_t sz = const_cast<subdomain_t&>(cell).getActiveLayerSize();
    return ((field.NRows() == sz.x + 2) && (field.NCols() == sz.y + 2));
}

#if MY_MULTISCALE_METHOD == 1
/**
 * Function copies the peer subdomain boundary (bin) to the current subdomain
 * boundary (bout) rescaling the array of boundary values accordingly. Three
 * situations are possible: (1) both boundaries have the same size; (2) this
 * boundary has twice less points than the other one; (3) this boundary has
 * twice more points then the other one.
 * While doing downscaling, we apply the following anti-aliasing filter:
 * [1 2 2 1]*(1/6). While doing upscaling, the end point values are copied
 * as is, but the internal input points are mapped between successive
 * couples of the output points. Then, the values at the output points are
 * linearly interpolated using the mapped positions.
 */
void AdjustBoundary(double_array_t & bout, const double_array_t & bin, long N)
{
    if (bin.size() == size_t(N)) {          // same size
        bout = bin;
    } else if (bin.size() == 2*size_t(N)) { // down-scaling by factor 2
        bout.resize(size_t(N));
        const double scale = 1.0/6.0;
        bout[0] = scale * (bin[0] + 2*(bin[0] + bin[1]) + bin[2]);
        for (long i = 1; i < N - 1; ++i) {
            long q = 2*i;
            bout[i] = scale * (bin[q-1] + 2*(bin[q] + bin[q+1]) + bin[q+2]);
        }
        long q = 2*(N-1);
        bout[N-1] = scale * (bin[q-1] + 2*(bin[q] + bin[q+1]) + bin[q+1]);
    } else if (2*bin.size() == size_t(N)) { // up-scaling by factor 2
        bout.resize(size_t(N));
        bout[0] = bin.front();
        bout[N-1] = bin.back();
        for (long i = 1; i < N - 1; ++i) {
            long q = (i-1)/2;
            bout[i] = bin[q] + (bin[q+1] - bin[q]) * (0.5*i - 0.25 - q);
        }
    } else {
        assert_true(false);
    }
}
#endif  // MY_MULTISCALE_METHOD

/**
 * Function computes: z = H * observations(t). Since H is a simple 0/1 matrix
 * that just picks up the observations at sensor locations, instead of
 * matrix-vector multiplication we get the observations directly.
 */
void GetObservations(Vector & z, const Matrix & observations, int timestep,
                     const Matrix & H, const point_array_t & sensors,
                     const size2d_t & layer_size)
{
    const index_t n = observations.NCols();
    assert_eq(z.Size(), n);
    for (index_t i = 0; i < n; ++i) { z(i) = observations(timestep, i); }

    (void) H; (void) sensors; (void) layer_size;
#ifdef AMDADOS_DEBUGGING
    #warning "Some extra validity test"
    // Mind the extended subdomain: one extra point layer on either side.
    assert_true(H.NCols() == (layer_size.x + 2) * (layer_size.y + 2));
    Matrix subfield(layer_size.x + 2, layer_size.y + 2);
    for (index_t i = 0; i < static_cast<index_t>(sensors.size()); ++i) {
        subfield(sensors[i].x + 1, sensors[i].y + 1) =
            observations(timestep, i);
    }
    Vector _z(n);
    MatVecMult(_z, H, subfield);    // _z = H * observations(t)
    assert_true(std::equal(_z.begin(), _z.end(), z.begin(),  z.end()));
#endif
}

/**
 * Function applies Dirichlet zero boundary condition at the outer border
 * of the domain.
 */
void ApplyBoundaryCondition(domain_t & state, const point2d_t & idx)
{
    subdomain_t &  subdom = state[idx];
    const size2d_t subdom_size = subdom.getActiveLayerSize();

    const index_t Nx = state.size().x, Sx = subdom_size.x;
    const index_t Ny = state.size().y, Sy = subdom_size.y;

    // Set the leftmost and rightmost.
    if (idx.x == 0)      subdom.setBoundary(Direction::Left,  double_array_t(Sy, 0.0));
    if (idx.x == Nx - 1) subdom.setBoundary(Direction::Right, double_array_t(Sy, 0.0));

    // Set the bottommost and topmost.
    if (idx.y == 0)      subdom.setBoundary(Direction::Down, double_array_t(Sx, 0.0));
    if (idx.y == Ny - 1) subdom.setBoundary(Direction::Up,   double_array_t(Sx, 0.0));
}

/**
 * Function computes the initial process model covariance matrix P
 * based on exponential distance.
 */
void InitialCovar(const Configuration & conf, Matrix & P)
{
    // Here we express correlation distances in logical coordinates
    // of nodal points.
    const double variance = conf.asDouble("model_ini_var");
    const double covar_radius = conf.asDouble("model_ini_covar_radius");
    const double sigma_x = std::max(covar_radius, 1.0);
    const double sigma_y = std::max(covar_radius, 1.0);
    const index_t Rx = Round(std::ceil(4.0 * sigma_x));
    const index_t Ry = Round(std::ceil(4.0 * sigma_y));
    const index_t Sx = conf.asInt("subdomain_x");
    const index_t Sy = conf.asInt("subdomain_y");
    const size2d_t layer_size(Sx, Sy);      // size at the finest resolution

    assert_true(P.IsSquare());
    assert_true(P.NRows() == (Sx + 2) * (Sy + 2));

    // Mind the extended subdomain: one extra point layer on either side.
    Fill(P, 0.0);
    for (index_t u = 0; u < Sx + 2; ++u) {
    for (index_t v = 0; v < Sy + 2; ++v) {
        index_t i = sub2ind(u, v, layer_size);
        double dx = 0.0, dy = 0.0;
        for (index_t x = u-Rx; x <= u+Rx; ++x) { if ((0 <= x) && (x < Sx + 2)) {
        for (index_t y = v-Ry; y <= v+Ry; ++y) { if ((0 <= y) && (y < Sy + 2)) {
            index_t j = sub2ind(x, y, layer_size);
            if (i <= j) {
                dx = (u - x) / sigma_x;
                dy = (v - y) / sigma_y;
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
    for (index_t k = 0; k < Q.NRows(); ++k) {
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
    for (index_t k = 0; k < R.NRows(); ++k) {
        R(k,k) = 1.0 + model_noise_R * distrib(gen);    // always >= 1
    }
}

/**
 * Function initializes the observation matrix H of size:
 * (number of sensors in subdomain) x (number of nodal points in subdomain).
 */
void ComputeH(const point_array_t & sensors, const size2d_t & layer_size,
              Matrix & H)
{
    if (sensors.empty()) {
        H.Clear();
        return;
    }
    // Mind the extended subdomain: one extra point layer on either side.
    assert_true(H.NRows() == static_cast<index_t>(sensors.size()));
    assert_true(H.NCols() == (layer_size.x + 2) * (layer_size.y + 2));
    Fill(H, 0.0);
    for (size_t k = 0; k < sensors.size(); ++k) {
        H(k, sub2ind(sensors[k].x + 1, sensors[k].y + 1, layer_size)) = 1.0;
    }
}

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
                        const flow_t & flow, const size2d_t & layer_size,
                        unsigned resolution, int dt_divisor = 0)
{
    const index_t Sx = layer_size.x;
    const index_t Sy = layer_size.y;

    assert_true((B.NRows() == B.NCols()) && (B.NCols() == (Sx + 2)*(Sy + 2)));

    // Important: scale the space steps according to resolution.
    const double resol_ratio =
        (resolution == LayerFine) ? 1.0 : conf.asDouble("resolution_ratio");

    const double D  = conf.asDouble("diffusion_coef");
    const double dx = conf.asDouble("dx") * resol_ratio;
    const double dy = conf.asDouble("dy") * resol_ratio;
    const double dt = conf.asDouble("dt")/((dt_divisor > 0) ? dt_divisor : 1);

    const double rho_x = D * dt / std::pow(dx,2);
    const double rho_y = D * dt / std::pow(dy,2);

    const double v0x = 2.0 * dx / dt;
    const double v0y = 2.0 * dy / dt;

    const double vx = flow.first  / v0x;
    const double vy = flow.second / v0y;

    // The internal and the boundary points of extended subdomain are treated
    // differently. The border values are passed through as is (B(i,i) = 1).
    // Mind the extended subdomain: one extra point layer on either side.
    MakeIdentityMatrix(B);
    for (index_t x = 1; x <= Sx; ++x) {
    for (index_t y = 1; y <= Sy; ++y) {
        index_t i = sub2ind(x, y, layer_size);
        B(i,i) = 1.0 + 2*(rho_x + rho_y);
        B(i,sub2ind(x-1, y, layer_size)) = - vx - rho_x;
        B(i,sub2ind(x+1, y, layer_size)) = + vx - rho_x;
        B(i,sub2ind(x, y-1, layer_size)) = - vy - rho_y;
        B(i,sub2ind(x, y+1, layer_size)) = + vy - rho_y;
    }}
}

/**
 * Function copies an Allscale subdomain to the matrix. The output matrix
 * represents so called "extended subdomain" where one extra point layer on
 * either side is added by copying values from the peer subdomains' boundaries.
 * When a subdomain is located at the outer boundary of the whole domain,
 * we initialize the extended points in a way that provides zero boundary
 * condition on the density derivative along the normal: du/dn = 0. For example,
 * at the left boundary we set (see the code): field(0,:) = field(2,:);
 * Here field(0,:) addresses the points outside the domain (because this is an
 * extended subdomain) and field(1,:) addresses the points on the global outer
 * boundary. Derivative at the left outer boundary reads:
            d(field(1,y))/dx = (field(2,y) - field(0,y))/2 = 0,
 * according to above condition.
 */
void MatrixFromAllscale(Matrix & field,
                        const domain_t & dom, const point2d_t & idx)
{
    const unsigned layer_no = dom[idx].getActiveLayer();
    const index_t Nx = dom.size().x;
    const index_t Ny = dom.size().y;
    const index_t Sx = const_cast<subdomain_t&>(dom[idx]).getActiveLayerSize().x;
    const index_t Sy = const_cast<subdomain_t&>(dom[idx]).getActiveLayerSize().y;

    // Copy the internal points of a subdomain to the (internal part of) output
    // field. Mind the extended subdomain: an extra point layer on either side.
    assert_true(CheckSizes(dom[idx], field));
    dom[idx].forAllActiveNodes([&](const point2d_t & pos, const double & val) {
        field(pos.x + 1, pos.y + 1) = val;
    });

    double_array_t boundary;    // placeholder of boundary points


    // Set up left-most points from the right boundary of the left peer.
    if (idx.x > 0) {
#if MY_MULTISCALE_METHOD == 1
        AdjustBoundary(boundary, dom[{idx.x-1, idx.y}].getBoundary(Direction::Right), Sy);
#else // method == 2
        boundary = dom[{idx.x-1, idx.y}].data.getBoundary(layer_no, Direction::Right);
#endif
        assert_true(boundary.size() == size_t(Sy));
        for (index_t y = 0; y < Sy; ++y) field(0, y+1) = boundary[y];
    } else {
        for (index_t y = 0; y < Sy; ++y) field(0, y+1) = field(2, y+1);
    }

    // Set up right-most points from the left boundary of the right peer.
    if (idx.x+1 < Nx) {
#if MY_MULTISCALE_METHOD == 1
        AdjustBoundary(boundary, dom[{idx.x+1, idx.y}].getBoundary(Direction::Left), Sy);
#else // method == 2
        boundary = dom[{idx.x+1, idx.y}].data.getBoundary(layer_no, Direction::Left);
#endif
        assert_true(boundary.size() == size_t(Sy));
        for (index_t y = 0; y < Sy; ++y) field(Sx+1, y+1) = boundary[y];
    } else {
        for (index_t y = 0; y < Sy; ++y) field(Sx+1, y+1) = field(Sx-1, y+1);
    }

    // Set up bottom-most points from the top boundary of the bottom peer.
    if (idx.y > 0) {
#if MY_MULTISCALE_METHOD == 1
        AdjustBoundary(boundary, dom[{idx.x, idx.y-1}].getBoundary(Direction::Up), Sx);
#else // method == 2
        boundary = dom[{idx.x, idx.y-1}].data.getBoundary(layer_no, Direction::Up);
#endif
        assert_true(boundary.size() == size_t(Sx));
        for (index_t x = 0; x < Sx; ++x) field(x+1, 0) = boundary[x];
    } else {
        for (index_t x = 0; x < Sx; ++x) field(x+1, 0) = field(x+1, 2);
    }

    // Set up top-most points from the bottom boundary of the top peer.
    if (idx.y+1 < Ny) {
#if MY_MULTISCALE_METHOD == 1
        AdjustBoundary(boundary, dom[{idx.x, idx.y+1}].getBoundary(Direction::Down), Sx);
#else // method == 2
        boundary = dom[{idx.x, idx.y+1}].data.getBoundary(layer_no, Direction::Down);
#endif
        assert_true(boundary.size() == size_t(Sx));
        for (index_t x = 0; x < Sx; ++x) field(x+1, Sy+1) = boundary[x];
    } else {
        for (index_t x = 0; x < Sx; ++x) field(x+1, Sy+1) = field(x+1, Sy-1);
    }

    // The corner points are not used in finite-difference scheme applied to
    // internal subdomain points, however, Kalman filter uses the entire
    // subdomain field (as currently implemented but could be easily avoided).
    // For the latter reason, we have to assign some feasible values to the
    // corner points. Note, since we operate on extended subdomains, their
    // corners belong to unreachable diagonal peers.
    field(0,0)       = (field(1,0)     + field(1,1)   + field(0,1)    ) / 3.0;
    field(0,Sy+1)    = (field(1,Sy+1)  + field(1,Sy)  + field(0,Sy)   ) / 3.0;
    field(Sx+1,0)    = (field(Sx,0)    + field(Sx,1)  + field(Sx+1,1) ) / 3.0;
    field(Sx+1,Sy+1) = (field(Sx,Sy+1) + field(Sx,Sy) + field(Sx+1,Sy)) / 3.0;
}

/**
 * Function copies a matrix, which represents an extended subdomain,
 * to Allscale subdomain structure.
 */
void AllscaleFromMatrix(subdomain_t & cell, const Matrix & field)
{
    // Mind the extended subdomain: one extra point layer on either side.
    assert_true(CheckSizes(cell, field));
    cell.forAllActiveNodes([&](const point2d_t & pos, double & val) {
        val = field(pos.x + 1, pos.y + 1);
    });
}

/**
 * Function is invoked for each sub-domain, which contains at least one sensor,
 * during the time integration. For such a subdomain the Kalman filter governs
 * the simulation by pulling it towards the observed ground-truth.
 */
const subdomain_t & SubdomainRoutineKalman(
                            const Configuration & conf,
                            const point_array_t & sensors,
                            const Matrix        & observations,
                            const size_t          timestamp,
                            const domain_t      & curr_state,
                            domain_t            & next_state,
                            SubdomainContext    & ctx, 
                            const point2d_t     & idx,
                            const size_t          Nsubiter,
                            const size_t          Nt)
{
    const unsigned resolution = static_cast<unsigned>(LayerFine);
    assert_true(&curr_state != &next_state);
    assert_true(curr_state.size() == next_state.size());

    // Important: synchronize active layers, otherwise when we exchange data
    // at the borders of neighbour subdomains the size mismatch can happen.
    assert_true(curr_state[idx].getActiveLayer() == resolution);
    next_state[idx].setActiveLayer(resolution);
    const size2d_t layer_size =
        const_cast<subdomain_t&>(curr_state[idx]).getActiveLayerSize();

    // Get the discrete time (index of iteration) in the range [0..Nt) and
    // the index of sub-iteration in the range [0..Nsubiter).
    assert_true(timestamp < Nt * Nsubiter);
    const size_t t_discrete = timestamp / Nsubiter;
    const size_t sub_iter   = timestamp % Nsubiter;

#ifdef AMDADOS_DEBUGGING    // printing progress
    if ((idx == point2d_t(0,0)) && (sub_iter == 0)) {
        if (t_discrete == 0) std::cout << "Nt: " << Nt << std::endl;
        std::cout << "\rtime: " << t_discrete << std::flush;
        if (t_discrete + 1 == Nt) std::cout << std::endl << std::flush;
    }
#endif

    // Compute flow velocity vector.
    ctx.flow = Flow(conf, t_discrete);

    // Copy state field into the matrix object.
    MatrixFromAllscale(ctx.field, curr_state, idx);

    // At the beginning of a regular iteration (i.e. at the first
    // sub-iteration) do: (1) get the new observations; (2) obtain new noise
    // covariance matrices; (3) compute the prior state estimation.
    if (sub_iter == 0) {
        // Get the current sensor measurements.
        GetObservations(ctx.z, observations, static_cast<int>(t_discrete),
                        ctx.H, sensors, layer_size);

        // Covariance matrices can change over time.
        ComputeQ(conf, ctx.Q);
        ComputeR(conf, ctx.R);

        // Prior estimation.
        InverseModelMatrix(ctx.B, conf, ctx.flow, layer_size, resolution);
        ctx.Kalman.PropagateStateInverse(ctx.field, ctx.P, ctx.B, ctx.Q);
    }

    // Filtering by Kalman filter.
    ctx.Kalman.SolveFilter(ctx.field, ctx.P, ctx.H, ctx.R, ctx.z);

    // Put new estimation back to the Allscale state field.
    AllscaleFromMatrix(next_state[idx], ctx.field);

    // Ensure boundary conditions on the outer border.
    ApplyBoundaryCondition(next_state, idx);

    // Ensure non-negative (physically plausible) density.
    next_state[idx].forAllActiveNodes([](double & v){ if (v < 0.0) v = 0.0; });

#if MY_MULTISCALE_METHOD == 2
    // Make up the coarse layer (so the peer subdomains can use either
    // fine or low resolution one), then go back to the default resolution.
    next_state[idx].coarsen([](const double & elem) { return elem; });
    next_state[idx].setActiveLayer(resolution);
#endif

    return next_state[idx];
}

/**
 * Function is invoked for each sub-domain without sensors therein
 * during the time integration.
 */
const subdomain_t & SubdomainRoutineNoSensors(
                            const Configuration & conf,
                            const size_t          timestamp,
                            const domain_t      & curr_state,
                            domain_t            & next_state,
                            SubdomainContext    & ctx,
                            const point2d_t     & idx,
                            const size_t          Nsubiter,
                            const size_t          Nt)
{
    const unsigned resolution = static_cast<unsigned>(LayerLow);
    assert_true(&curr_state != &next_state);
    assert_true(curr_state.size() == next_state.size());

    // Important: synchronize active layers, otherwise when we exchange data
    // at the borders of neighbour subdomains the size mismatch can happen.
    assert_true(curr_state[idx].getActiveLayer() == resolution);
    next_state[idx].setActiveLayer(resolution);
    const size2d_t layer_size =
        const_cast<subdomain_t&>(curr_state[idx]).getActiveLayerSize();

    // Get the discrete time (index of iteration) in the range [0..Nt) and
    // the index of sub-iteration in the range [0..Nsubiter).
    assert_true(timestamp < Nt * Nsubiter);
    const size_t t_discrete = timestamp / Nsubiter;
    const size_t sub_iter   = timestamp % Nsubiter;  (void) sub_iter;

#ifdef AMDADOS_DEBUGGING    // printing progress
    if ((idx == point2d_t(0,0)) && (sub_iter == 0)) {
        if (t_discrete == 0) std::cout << "Nt: " << Nt << std::endl;
        std::cout << "\rtime: " << t_discrete << std::flush;
        if (t_discrete + 1 == Nt) std::cout << std::endl << std::flush;
    }
#endif

    // Compute flow velocity vector.
    ctx.flow = Flow(conf, t_discrete);

    // Copy state field into the matrix object.
    MatrixFromAllscale(ctx.field, curr_state, idx);

    // Prior estimation.
    InverseModelMatrix(ctx.B, conf, ctx.flow, layer_size,
                       resolution, static_cast<int>(Nsubiter));
    ctx.tmp_field = ctx.field;              // copy state into a temporary one
    ctx.LU.Init(ctx.B);                     // decompose: B = L*U
    ctx.LU.Solve(ctx.field, ctx.tmp_field); // new_field = B^{-1}*old_field

    // Put the estimation back to the Allscale state field.
    AllscaleFromMatrix(next_state[idx], ctx.field);

    // Ensure boundary conditions on the outer border.
    ApplyBoundaryCondition(next_state, idx);

    // Ensure non-negative (physically plausible) density.
    next_state[idx].forAllActiveNodes([](double & v){ if (v < 0.0) v = 0.0; });

#if MY_MULTISCALE_METHOD == 2
    // Make up the fine layer (so the peer subdomains can use either
    // fine or low resolution one), then go back to the default resolution.
    next_state[idx].refine([](const double & elem) { return elem; });
    next_state[idx].setActiveLayer(resolution);
#endif

    return next_state[idx];
}

} // anonymous namespace

/**
 * Using model matrix A, the function integrates advection-diffusion equation
 * forward in time inside individual subdomains and records all the solutions
 * as the state fields. Sub-iterations are applied in order to make the
 * global solution seam-less along subdomain boundaries. On top of that, the
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

    const point2d_t GridSize = GetGridSize(conf);   // size in subdomains
    const size_t    Nt = conf.asUInt("Nt");
    const size_t    Nsubiter = conf.asUInt("num_sub_iter");
	//const size_t    Nwrite = std::min(Nt, conf.asUInt("write_num_fields"));

    context_domain_t contexts(GridSize);    // variables of each sub-domain
    domain_t         temp_field(GridSize);  // grid of sub-domains
    domain_t         state_field(GridSize); // grid of sub-domains

    // Initialize the observation and model covariance matrices.
    pfor(point2d_t(0,0), GridSize, [&,conf](const point2d_t & idx) {
        // Zero field at the beginning for all the resolutions.
        static_assert(LayerFine <= LayerLow, "");
        //for (int layer = LayerFine; layer <= LayerLow; ++layer) {
        //    state_field[idx].setActiveLayer(layer);
        //     temp_field[idx].setActiveLayer(layer);
        //    state_field[idx].forAllActiveNodes([](double & v) { v = 0.0; });
        //     temp_field[idx].forAllActiveNodes([](double & v) { v = 0.0; });
        //    ApplyBoundaryCondition(state_field, idx);
        //    ApplyBoundaryCondition( temp_field, idx);
        //}

        // If there is at least one sensor in a subdomain, then we operate at
        // the fine resolution, otherwise at the low resolution.
        const index_t Nsensors = static_cast<index_t>(sensors[idx].size());
        if (Nsensors > 0) {
            state_field[idx].setActiveLayer(LayerFine);
             temp_field[idx].setActiveLayer(LayerFine);
        } else {
            state_field[idx].setActiveLayer(LayerLow);
             temp_field[idx].setActiveLayer(LayerLow);
        }

        const size2d_t layer_size = state_field[idx].getActiveLayerSize();
        const index_t Sx = layer_size.x;
        const index_t Sy = layer_size.y;
        const index_t sub_prob_size = (Sx + 2) * (Sy + 2);

        // Note, we initialize the Kalman filter matrices only in the
        // presence of sensor(s), otherwise they are useless. Also,
        // mind the extended subdomain: one extra point layer on either side.
        SubdomainContext & ctx = contexts[idx];
        ctx.field.Resize(Sx + 2, Sy + 2);
        ctx.B.Resize(sub_prob_size, sub_prob_size);
        if (Nsensors > 0) {
            ctx.P.Resize(sub_prob_size, sub_prob_size);
            ctx.Q.Resize(sub_prob_size, sub_prob_size);
            ctx.H.Resize(Nsensors, sub_prob_size);
            ctx.R.Resize(Nsensors, Nsensors);
            ctx.z.Resize(Nsensors);
            ctx.sensors = sensors[idx];
            ComputeH(sensors[idx], layer_size, ctx.H);
            InitialCovar(conf, ctx.P);
        }
    });

    // Time integration forward in time. We want to make Nt (normal) iterations
    // and Nsubiter sub-iterations within each (normal) iteration.
    ::allscale::api::user::algorithm::stencil<allscale::api::user::algorithm::implementation::coarse_grained_iterative>(
        state_field, Nt * Nsubiter,
        [&,conf,Nsubiter,Nt](time_t t, const point2d_t & idx, const domain_t & state)
        -> const subdomain_t &  // cell is not copy-constructible, so '&'
        {
            if (contexts[idx].sensors.size() > 0) {
                return SubdomainRoutineKalman(conf, sensors[idx],
                            observations[idx], size_t(t),
                            state, temp_field, contexts[idx], idx, Nsubiter, Nt);
            } else {
                return SubdomainRoutineNoSensors(conf, size_t(t),
                            state, temp_field, contexts[idx], idx, Nsubiter, Nt);
            }
        },
        // Monitoring.
        ::allscale::api::user::algorithm::observer(
            // Time filter: choose time-slices evenly distributed on time axis.
            [=](time_t /*t*/) {
                // Filter out the first sub-iteration, skip the others.
                //if ((t % time_t(Nsubiter)) != 0) return false;
                //t /= time_t(Nsubiter);
                //return ((((Nwrite-1)*(t-1))/(Nt-1) != ((Nwrite-1)*t)/(Nt-1)));
				return false;
            },
            // Space filter: no specific points.
            [](const point2d_t &) { return true; },
            // Append a full field to the file of simulation results.
            [&,Nsubiter](time_t /*t*/, const point2d_t & /*idx*/, const subdomain_t & /*cell*/) {
//                t /= time_t(Nsubiter);
//                // Important: we save field at the finest resolution.
//                subdomain_t temp;
//                temp = cell;
//                while (temp.getActiveLayer() != LayerFine) {
//                    temp.refine([](const double & elem) { return elem; });
//                }
//                const size2d_t finest_layer_size = temp.getActiveLayerSize();
//                // Write the subdomain into file.
//                temp.forAllActiveNodes([&](const point2d_t & loc, double val) {
//                    point2d_t glo = Sub2Glo(loc, idx, finest_layer_size);
//                    out_stream.atomic([=](auto & file) {
//                        file << t << " " << glo.x << " " << glo.y << " " << val << "\n";
//                    });
//                });
            }
        )
    );

	std::string filename = MakeFileName(conf, "field");
	::allscale::api::user::algorithm::async([=,&state_field](){
		// Open file manager and the output file for writing.
		FileIOManager & file_manager = FileIOManager::getInstance();
		Entry stream_entry = file_manager.createEntry(filename, Mode::Text);
		auto out_stream = file_manager.openOutputStream(stream_entry);

		for(index_t i = 0; i < GridSize.x; ++i) {
			for(index_t j = 0; j < GridSize.y; ++j) {
				const point2d_t idx{ i,j };
				const size_t t = Nt - 1;
				subdomain_t temp;
				temp = state_field[idx];
				while(temp.getActiveLayer() != LayerFine) {
					temp.refine([](const double & elem) { return elem; });
				}
				const size2d_t finest_layer_size = temp.getActiveLayerSize();
				// Write the subdomain into file.
				temp.forAllActiveNodes([&](const point2d_t & loc, double val) {
					point2d_t glo = Sub2Glo(loc, idx, finest_layer_size);
					out_stream << t << " " << glo.x << " " << glo.y << " " << val << "\n";
				});
			}
		}
    		file_manager.close(out_stream);
		// need to output result file name for the CI system to pick it up
	}).wait();
	std::cout << "Wrote result data to " << filename << std::endl;

    MY_INFO("%s", "\n\n")
}

/**
 * Function initializes dependent parameters given the primary ones
 * specified by user.
 */
void InitDependentParams(Configuration & conf)
{
    // Get resolution at the finest level by creating a temporary subdomain.
    subdomain_t temp;
    temp.setActiveLayer(LayerFine);
    const size2d_t fine_size = temp.getActiveLayerSize();
    const index_t Sx = fine_size.x;   // subdomain size x
    const index_t Sy = fine_size.y;   // subdomain size y

    // Compute the size ratio between fine and low resolution layers.
    // Same ratio is expected in both dimensions.
    temp.setActiveLayer(LayerLow);
    const size2d_t low_size = temp.getActiveLayerSize();
    assert_true(fine_size.x * low_size.y == fine_size.y * low_size.x);
    conf.SetDouble("resolution_ratio", static_cast<double>(fine_size.x) /
                                       static_cast<double>( low_size.x));

    // Check some global constants.
    static_assert((0 <= Direction::Up   ) && (Direction::Up    < NSIDES), "");
    static_assert((0 <= Direction::Down ) && (Direction::Down  < NSIDES), "");
    static_assert((0 <= Direction::Left ) && (Direction::Left  < NSIDES), "");
    static_assert((0 <= Direction::Right) && (Direction::Right < NSIDES), "");
    assert_true((Sx >= 3) && (Sy >= 3)) << "subdomain must be at least 3x3";

    // Ensure integer values for certain parameters.
    assert_true(conf.IsInteger("num_subdomains_x"));
    assert_true(conf.IsInteger("num_subdomains_y"));
    assert_true(conf.IsInteger("subdomain_x"));
    assert_true(conf.IsInteger("subdomain_y"));
    assert_true(conf.IsInteger("integration_nsteps"));

    // Check the subdomain size: hard-coded value must match the parameter.
    assert_true(conf.asInt("subdomain_x") == Sx)
                        << "subdomain_x mismatch" << std::endl;
    assert_true(conf.asInt("subdomain_y") == Sy)
                        << "subdomain_y mismatch" << std::endl;

    const point2d_t grid_size = GetGridSize(conf);
    const index_t nx = grid_size.x * Sx;                // global X-size
    const index_t ny = grid_size.y * Sy;                // global Y-size

    const double D = conf.asDouble("diffusion_coef");
    assert_true(D > 0.0);

    conf.SetInt("global_problem_size", static_cast<int>(nx * ny));
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
 * The main function of this application runs simulation with data
 * assimilation using method to handle domain subdivision.
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

    // Load sensor data obtained from Python code.
    Grid<point_array_t,2> sensors(GetGridSize(conf));
    Grid<Matrix,2> observations(GetGridSize(conf));
    LoadSensorLocations(conf, sensors);
    LoadSensorMeasurements(conf, sensors, observations);

    // Run the simulation with data assimilation. Important: by this time
    // some parameters had been initialized in InitDependentParams(..), so
    // we can safely proceed to the main part of the simulation algorithm.
    MY_TIME_IT("Running the simulation with data assimilation ...")
    RunDataAssimilation(conf, sensors, observations);
}

} // namespace amdados
