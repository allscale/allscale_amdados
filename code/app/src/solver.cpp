//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include "allscale/api/user/operator/pfor.h"
#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/geometry.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/sparse_matrix.h"
#include "amdados/app/utils/configuration.h"
#include "amdados/app/utils/cholesky.h"
#include "amdados/app/model/i_model.h"
#include "amdados/app/utils/kalman_filter.h"
#include "amdados/app/model/euler_finite_diff.h"

using namespace amdados::app::utils;

namespace amdados {
namespace app {

using ::allscale::api::user::pfor;
using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;

// Global variables are declared here.
// TODO: how does it fit into MPI framework???
//namespace {

// Total number of elements (or nodal points) in the entire domain in each dimension.
const int GLOBAL_NELEMS_X = NELEMS_X * NUM_DOMAINS_X;
const int GLOBAL_NELEMS_Y = NELEMS_Y * NUM_DOMAINS_Y;

// Number of elements (or nodal points) in a sub-domain.
const size_t SUB_PROBLEM_SIZE = NELEMS_X * NELEMS_Y;

// Number of measurement locations in a sub-domain.
const size_t NUM_MEASUREMENTS = NELEMS_X * NELEMS_Y;

using subdomain_sub2ind_t = Sub2Ind<NELEMS_X, NELEMS_Y>;
using subdomain_ind2sub_t = Ind2Sub<NELEMS_X, NELEMS_Y>;

using point2d_t = GridPoint<2>;
using size2d_t = point2d_t;

using point3d_t = GridPoint<3>;
using size3d_t = point3d_t;

using cube_t = Grid<double,3>;
using DA_vector_t = Vector<SUB_PROBLEM_SIZE>;
using DA_matrix_t = Matrix<SUB_PROBLEM_SIZE, SUB_PROBLEM_SIZE>;
using DA_sp_matrix_t = SpMatrix<SUB_PROBLEM_SIZE,SUB_PROBLEM_SIZE>;

using sub_domain_t = Cell<double,sub_domain_config_t>;
using vec_grid_t = Grid<DA_vector_t,2>;
using mat_grid_t = Grid<DA_matrix_t,2>;

// Zero-coordinate point and global grid size counted in sub-domains.
const point2d_t Origin = 0;
const size2d_t  SubDomGridSize = {NUM_DOMAINS_X, NUM_DOMAINS_Y};

// Create the overall grid. Two similar grids (A and B) are created
// for convenience to alternate between successive states.
//XXX:???: A is of form A[{ndox,ndomy}].layer[{xElCount,yElcount}]
Grid<sub_domain_t,2> gridA(SubDomGridSize);
Grid<sub_domain_t,2> gridB(SubDomGridSize);

// Data structures for observation data has enough room for all the nodal
// points accross the whole domain at all timestamps.
unique_ptr<cube_t> obsv_glob;

// Create data structure for storing data assimilation matrices.
mat_grid_t P(SubDomGridSize);    // state covariance matrix
mat_grid_t R(SubDomGridSize);    // observation noise covariance matrix
mat_grid_t Q(SubDomGridSize);    // process noise covariance matrix
mat_grid_t H(SubDomGridSize);    // projector from model space to observation space

// Create data structure for storing data assimilation vectors.
vec_grid_t Obvs(SubDomGridSize);
vec_grid_t forecast(SubDomGridSize);
vec_grid_t posterior(SubDomGridSize);

unique_ptr< IModel<NELEMS_X,NELEMS_Y> > model;

// Create Kalman filters for each sub-domain.
using Kalman_t = amdados::app::utils::KalmanFilter<SUB_PROBLEM_SIZE, NUM_MEASUREMENTS>;
using Kalman_ptr_t = unique_ptr<Kalman_t>;
Grid<Kalman_ptr_t,2> kalman_filters(SubDomGridSize);

//-------------------------------------------------------------------------------------------------
// Function reads observations from a file.
// TODO: right now assume data available every timestep; need to update
//-------------------------------------------------------------------------------------------------
void ReadObservations(cube_t & obsver, const std::string & filename, int observint)
{
    assert_true(observint == 1) << "so far, the observation interval must be equal to 1" << endl;
    std::ifstream file(filename);
    assert_true(file.good()) << "failed to open observation file: " << filename << endl;
    int debug_t = 0;
    for (int t = 0; file >> t;) {       // header line contains a time stamp
        assert_true(debug_t == t) << "time step missed in observations: " << debug_t << endl;
        ++debug_t;
        for (int y = 0; y < GLOBAL_NELEMS_Y; ++y) {
        for (int x = 0; x < GLOBAL_NELEMS_X; ++x) {   // XXX: x is the fastest
            int i = 0, j = 0;
            double val = 0;
            file >> i >> j >> val;
            if (!((i == x) && (j == y))) {
                assert_true(0) << "mismatch between the grid and observation layouts" << endl;
            }
            //assert_true((i == x) && (j == y))     // compilation BUG
            //<< "mismatch between the grid and observation layouts" << endl;
            obsver[{i,j,t}] = val;
        }}
    }
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// TODO: introduce a function to read flowfield data
// Read of form I,J,U,V and map to grid of amdados application
void ReadFlows(cube_t& /*flowfield*/,
        const std::string /*filename_flow*/, int /*nflopts_x*/ ,int /*nflopts_y*/)
{}

//-------------------------------------------------------------------------------------------------
// Utility function saves the current state.
//-------------------------------------------------------------------------------------------------
void SaveGrid2D(const Configuration & conf, const std::string & title, int t,
                const Grid<sub_domain_t,2> & grid)
{
    std::stringstream ss;
    ss << conf.asString("output_dir") << "/" << title
       << std::setfill('0') << std::setw(5) << t << ".txt";
    std::fstream f(ss.str(), std::ios::out | std::ios::trunc);
    assert_true(f.good()) << "failed to open the file for writing: " << ss.str() << endl;
    f << "# Layout: [1] dimensionality, [2..dim+1] sizes per dimension, [dim+2...] values" << endl;
    const size_t dim = 2;
    f << dim << endl << GLOBAL_NELEMS_X << endl << GLOBAL_NELEMS_Y << endl;
    for (int x = 0; x < GLOBAL_NELEMS_X; ++x) {
    for (int y = 0; y < GLOBAL_NELEMS_Y; ++y) { // XXX: y the fastest
        double value = grid[{x / NELEMS_X, y / NELEMS_Y}].getLayer<L_100m>()
                           [{x % NELEMS_X, y % NELEMS_Y}];
        f << value << endl;
    }}
    f.flush();
};

//-------------------------------------------------------------------------------------------------
// Function computes the model covariance matrix as function of exponential distance.
// XXX: what does it mean: Follow same structure as that implemented for Runga Kutta solver
//-------------------------------------------------------------------------------------------------
void GetModelCovar(const Configuration & conf, DA_matrix_t & covar)
{
    // TODO: check matrix sizes.
    const double var = conf.asDouble("model_covar_var");
    const double Rx = conf.asDouble("model_covar_Rx");
    const double Ry = conf.asDouble("model_covar_Ry");

    subdomain_ind2sub_t ind2sub;
    int                 x1 = 0, y1 = 0, x2 = 0, y2 = 0;

    for (int i = 0; i < static_cast<int>(SUB_PROBLEM_SIZE); ++i) {
        ind2sub(i, x1, y1);
        covar[{i,i}] = var;
        for (int k = i + 1; k < static_cast<int>(SUB_PROBLEM_SIZE); ++k) {
            ind2sub(k, x2, y2);
            covar[{k,i}] = covar[{i,k}] = var * exp(- Square((x1 - x2)/Rx) - Square((y1 - y2)/Ry));
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function computes the matrix of observations.
//-------------------------------------------------------------------------------------------------
template<size_t MSIZE>
void ComputeH(Matrix<MSIZE,MSIZE> & H)
{
    MakeIdentityMatrix(H);
}

//-------------------------------------------------------------------------------------------------
// Function computes the observation noise covariance matrix.
//-------------------------------------------------------------------------------------------------
template<size_t MSIZE>
void ComputeR(Matrix<MSIZE,MSIZE> & R)
{
    // Random values between 1 and 2, otherwise singular matrices can happen down the road.
    std::mt19937                           gen(std::time(nullptr));
    std::uniform_real_distribution<double> distrib(1.0, 2.0);

    for (int i = 0; i < static_cast<int>(MSIZE); ++i) {
        R[{i,i}] = distrib(gen);
    }
}

//-------------------------------------------------------------------------------------------------
// Function gets the observations at certain time and returns them unrolled into a vector.
//-------------------------------------------------------------------------------------------------
void GetObservation(DA_vector_t & obsver_vec, const cube_t & obsver,
                    int timestep, const point2d_t & domain)
{
    subdomain_sub2ind_t sub2ind;
    for (int x = 0; x < NELEMS_X; ++x) {
        int xg = x + domain[0] * NELEMS_X;                          // global abscissa
        for (int y = 0; y < NELEMS_Y; ++y) {
            int yg = y + domain[1] * NELEMS_Y;                      // global ordinate
            obsver_vec[{sub2ind(x,y)}] = obsver[{xg,yg,timestep}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// TODO: what is this?
//-------------------------------------------------------------------------------------------------
double mu1(double timestep)
{
    return -0.6 * sin(timestep/10 - M_PI) * 0.2;
}

//-------------------------------------------------------------------------------------------------
// TODO: what is this?
//-------------------------------------------------------------------------------------------------
double mu2(double timestep)
{
    return -1.2 * sin(timestep/5 - M_PI) * 0.2;
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
void InitSolver(const Configuration & conf)
{
    // Check the geometry.
    assert_true(conf.asInt("num_domains_x") == NUM_DOMAINS_X) << "num_domains_x mismatch" << endl;
    assert_true(conf.asInt("num_domains_y") == NUM_DOMAINS_Y) << "num_domains_y mismatch" << endl;
    assert_true(conf.asInt("num_elems_x") == NELEMS_X) << "num_elems_x mismatch" << endl;
    assert_true(conf.asInt("num_elems_y") == NELEMS_Y) << "num_elems_y mismatch" << endl;

    // Get the rounded initial position of a spot of substance and check the bounds.
    const int spot_x = static_cast<int>(std::floor(conf.asDouble("spot_x") + 0.5));
    const int spot_y = static_cast<int>(std::floor(conf.asDouble("spot_y") + 0.5));
    assert_true(IsInsideRange(spot_x, 0, GLOBAL_NELEMS_X) &&
                IsInsideRange(spot_y, 0, GLOBAL_NELEMS_Y))
                    << "species spot is not inside the domain" << endl;

    // Set all the matrices and vectors to zero.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        FillMatrix(gridA[idx].getLayer<L_100m>(), 0.0);
        FillMatrix(gridB[idx].getLayer<L_100m>(), 0.0);
        FillMatrix(P[idx], 0.0);
        FillMatrix(R[idx], 0.0);
        FillMatrix(Q[idx], 0.0);
        FillMatrix(H[idx], 0.0);
        FillVector(Obvs[idx], 0.0);
        FillVector(forecast[idx], 0.0);
        FillVector(posterior[idx], 0.0);
    });

    // TODO: comment.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & pos) {
        // initialize all cells on the 100m resolution
        gridA[pos].setActiveLayer(L_100m);
        // initialize the concentration, rho = 0
        gridA[pos].forAllActiveNodes([](double & value) { value = 0.0; });
    });

    // Set initial covariance matrices in each sub-domain.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        GetModelCovar(conf, P[idx]);
    });

    // Allocate a placeholder for the entire observation data and read the observations.
    const int T = conf.asInt("integration_T");
    assert_true((0 < T) && (T < numeric_limits<int>::max()/10)) << "wrong T" << endl;
    obsv_glob.reset(new cube_t(size3d_t{GLOBAL_NELEMS_X, GLOBAL_NELEMS_Y, T + 1}));
    ReadObservations(*obsv_glob, conf.asString("observation_file"),
                                 conf.asInt("observation_interval"));

    // Bring in a high concentration of substance at some spot.
    gridA[{spot_x / NELEMS_X, spot_y / NELEMS_Y}].getLayer<L_100m>()
         [{spot_x % NELEMS_X, spot_y / NELEMS_Y}] = conf.asDouble("spot_density");

//  TODO: load the flows computed previously
//  utils::Grid<double,3> flows({ num_domains_x, num_domains_y, T });

    // Initialize the sub-domain filters.
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        kalman_filters[idx].reset(new Kalman_t());
    });

    // Initialize the model.
    model.reset(new EulerFiniteDifferenceModel<NELEMS_X,NELEMS_Y>(conf));
}

//-------------------------------------------------------------------------------------------------
// steps for the solver;
// 1) Need to initialize all domains to zero (& all layers I would think L100, L20 and L4)
// 2) Need to distribute work across processors for solution of advection diffusion
//      In this manner it would be relatively straightforward task to implement the solver within
//      each subdomain
// 3) Need to integrate data assimilation schemes from kalman_filter.h and filter.h
//      Data structures are in place and simply needs to be integrated with appropriate data from
//      advection-diffusion process
// 4) Need to integrate adaptive meshing capabilities
//      Right now this links with (3): If data is available then resolve solution to higher
//      fidelity grid (down to 20m) and assimilate at this scale. Then solution returns to 100m
//      grid and solution continues. For time being; advection-diffusion
//      always happens on the 100m grid and we resolve at higher resolution only for data
//      assimilation
//-------------------------------------------------------------------------------------------------
void Compute(const Configuration & conf)
{
    const int save_int = conf.asInt("output_every_nth_time_step");
    const double time_delta = conf.asInt("time_delta");
    const double space_delta = conf.asInt("space_delta");

    // Time integration.
    for (int t = 0, T = conf.asInt("integration_T"); t <= T; t++) {
        std::cout << "Time = " << t << endl << flush;
        if ((save_int > 0) && (t % save_int == 0)) SaveGrid2D(conf, "state", t, gridA);

        // Go from time t to t+1. Compute next time step => store it in B.
        pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
            // "curr" is defined as grid A[i,j] and distributed across shared memory.
            // "next" is the next state estimation.
            const auto & curr = gridA[idx];
            auto &       next = gridB[idx];

            // Initially two states coincide.
            next = curr;

            int    nx = 1;                      // TODO: description?
            int    ny = 1;                      // TODO: description?
            double timept = t * time_delta;     // TODO: description?
            double flowu = mu1(timept);         // TODO: description?
            double flowv = mu2(timept);         // TODO: description?

            for (Direction dir : { Up, Down, Left, Right }) {
                // Skip global boarder; for example, if direction == up and no neighbour to south.
                if (dir == Up    && idx[0] == 0)                   continue;
                if (dir == Down  && idx[0] == SubDomGridSize[0]-1) continue;
                if (dir == Left  && idx[1] == 0)                   continue;
                if (dir == Right && idx[1] == SubDomGridSize[1]-1) continue;

                // Obtain the local boundary.
                auto local_boundary = curr.getBoundary(dir);

                // Obtain the neighboring boundary.
                auto remote_boundary =
                    (dir == Up)   ? gridA[idx + point2d_t{-1,0}].getBoundary(Down)  : // remote boundary is bottom strip of neighbour
                    (dir == Down) ? gridA[idx + point2d_t{ 1,0}].getBoundary(Up)    : // remote boundary is top of neighbour
                    (dir == Left) ? gridA[idx + point2d_t{0,-1}].getBoundary(Right) : // remote boundary is left of domain
                                    gridA[idx + point2d_t{0, 1}].getBoundary(Left);

                // Compute local flow in domain to decide if flow in or out of domain.
                if (dir == Down) ny = -1; // flow into domain from below
                if (dir == Left) nx = -1; // flow into domain from left
                double flow_boundary =  nx*flowu + ny*flowv;
                // TODO: scale the boundary vectors to the same resolution

                // Compute updated boundary.
                assert_true(local_boundary.size() == remote_boundary.size());
                if (flow_boundary < 0) {  // then flow into domain need to update boundary with neighbour value
                    for (size_t i = 0; i < local_boundary.size(); i++) {
                        // for now, we just take the average
                        // need to update this to account for flow direction (Fearghal)
                        local_boundary[i] = remote_boundary[i];
                    }
                }
                //	std::cout << "DIRECTION: " << (dir == Up ? "Up" :
                //	(dir == Down ? "Down" : (dir == Left ? "Left" : "Right"))) << std::endl;

                // Update boundary in result.
                next.setBoundary(dir, local_boundary);
            }

            //// 2) run iterations of runge-kutta
            //double error = 1000;
            //while (error > 0.1) {           // runge kutta iteration on sub-domain
                //error = applyRungeKutta(res,flowu,flowv,delta,stepsize);
            //}

            // For now data assimilation is done at periodic timesteps.
            if (t % 1 == 0) {
                // Compute mapping matrix H to project from model to observation grid.
                ComputeH(H[idx]);
                // Estimate measurement noise covariance matrix.
                ComputeR(R[idx]);
                // Extract appropriate observation data for subdomain and time.
                GetObservation(Obvs[idx], *obsv_glob, t, idx);
                // Estimate process noise covariance matrix.
                Q[idx] = R[idx];    // TODO: this is stub: same noise covariances

                // Get current state.
                utils::Reshape2Dto1D<NELEMS_X, NELEMS_Y>(forecast[idx], next.getLayer<L_100m>());

                // Ligh-weight adapter helps to get the model matrix from within Kalman framework.
                //struct ModelAdapter {
                    //ModelAdapter()
                    //const SpMatrix<SUB_PROBLEM_SIZE,SUB_PROBLEM_SIZE> & operator()() {
                        //return model->ModelMatrix(flowu, flowv, time_delta, space_delta, t);
                    //}
                //}
                auto ModelAdapter = [&]() -> const DA_sp_matrix_t & {
                    const DA_sp_matrix_t & spmat =
                            model->ModelMatrix(flowu, flowv, time_delta, space_delta, t);
                    return spmat;
                };

                // Upadate the Kalman filter using the process model and current process matrices.
                Kalman_t * kf = kalman_filters[idx].get();
                kf->IterateWithModel(ModelAdapter, Q[idx], H[idx], R[idx], Obvs[idx],
                                     forecast[idx], P[idx]);
                //// Get a-posteriory estimationis for the state and covariance.
                //posterior[idx] = kf->GetStateVector();
                //P[idx] = kf->GetCovariance();

                // next state = posterior.
                utils::Reshape1Dto2D<NELEMS_X, NELEMS_Y>(next.getLayer<L_100m>(), forecast[idx]);
            }
        }); // end pfor

        // swap A and B
        //gridA.swap(gridB);
    }
}

//} // anonymous namespace

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
void AmdadosSolver(const Configuration & conf)
{
    InitSolver(conf);
    Compute(conf);
}

} // end namespace app
} // end namespace amdados

