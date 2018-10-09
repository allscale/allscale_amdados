//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

// This is pure MPI implementation without Allscale API at all.

#pragma once

// Ignore annoying warnings in MPI headers.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include <mpi.h>
#pragma GCC diagnostic pop

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <map>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <thread>

#ifndef AMDADOS_PLAIN_MPI
#define AMDADOS_PLAIN_MPI
#endif
#include "../include/amdados/app/debugging.h"
#include "../include/amdados/app/amdados_utils.h"
#include "../include/amdados/app/matrix.h"
#include "../include/amdados/app/configuration.h"
#include "../include/amdados/app/cholesky.h"
#include "../include/amdados/app/lu.h"
#include "../include/amdados/app/kalman_filter.h"
#include "mpi_basic.h"
#include "mpi_grid.h"
#include "mpi_subdomain.h"
#include "mpi_input_data.h"
#include "mpi_output.h"
#include "../include/amdados/app/sensors_generator.h"
// Include source files directly for simpler maintenance of the MPI project.
#include "../src/amdados_utils.cpp"
#include "../src/configuration.cpp"
#include "../src/matrix.cpp"

std::fstream gLogFile;            // global instance of log-file
bool gTestBoundExchange = false;  // 0/1: testing boundary exchange mechanism

namespace amdados {

void TestExchangeCorrectness(SubDomain * sd, long timestamp);

//-----------------------------------------------------------------------------
// Function parses command-line options. Returns 'true' when the help message
// was requested.
//-----------------------------------------------------------------------------
bool ParseCommandArgs(int argc, char ** argv,
                      std::string & scenario, std::string & config_file)
{
    scenario = "simulation";
    config_file = "amdados.conf";
    for (int a = 0; a < argc; ++a) {
        std::string token = argv[a];
        if (token == "--scenario") {
            if (++a < argc) {
                scenario = argv[a];
            }
        } else if (token == "--config") {
            if (++a < argc) {
                config_file = argv[a];
            }
        } else if (token == "--test") {
            if (++a < argc) {
                if (std::string(argv[a]) == "boundary_exchange") {
                    gTestBoundExchange = true;
                } else {
                    assert_true(0) << "unknown test: " << argv[a];
                }
            }
        } else if ((token == "--help") || (token == "-h")) {
            if (GetRank() == 0) {
                MY_LOG(INFO) << "\n\n***** Help:\n"
                    << "This program repeats functionality of the main\n"
                    << "Amdados application but uses pure MPI instead\n"
                    << "of Allscale API.\n"
                    << "Options:\n"
                    << "1. --config filename : path to configuration file;\n"
                    << "\t\t\t default is 'amdados.conf'.\n"
                    << "2. --scenario name : 'sensors' or 'simulation';\n"
                    << "\t\t\t default is 'simulation', in another scenario\n"
                    << "\t\t\t the new sensor locations are generated.\n"
                    << "3. --test boundary_exchange : test for subdomain\n"
                    << "\t\t\t boundary exchange mechanism\n\n\n";
            }
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------
// Function generates new sensor locations (given configuration) and stores
// them into a file in the output directory ('output_dir').
//-----------------------------------------------------------------------------
void ScenarioSensors(const Configuration & conf)
{
    assert_true(GetRank() == 0) << "sensors must be generated on rank 0 node";
    std::string filename = MakeFileName(conf, "sensors");
    std::fstream file(filename, std::ios::trunc | std::ios::out);
    assert_true(file.good()) << "failed to open file: " << filename;
    point_array_t sensors;
    SensorsGenerator().MakeSensors(conf, sensors);
    for (size_t k = 0; k < sensors.size(); ++k) {
        file << sensors[k].x << " " << sensors[k].y << std::endl;
    }
    file.flush();
}

//-----------------------------------------------------------------------------
// Function initializes dependent parameters given the primary ones
// specified by user.
//-----------------------------------------------------------------------------
void InitDependentParams(Configuration & conf)
{
    using ::amdados::TINY;

    const int Sx = conf.asInt("subdomain_x");   // subdomain size x
    const int Sy = conf.asInt("subdomain_y");   // subdomain size y
    conf.SetDouble("resolution_ratio", 1.0);    // MPI version does not use it

    // Check some global constants.
    static_assert((0 <= Up   ) && (Up    < NSides), "");
    static_assert((0 <= Down ) && (Down  < NSides), "");
    static_assert((0 <= Left ) && (Left  < NSides), "");
    static_assert((0 <= Right) && (Right < NSides), "");
    static_assert(NSides == 4, "");
    assert_true((Sx >= 3) && (Sy >= 3)) << "subdomain must be at least 3x3";

    // Ensure integer values for certain parameters.
    assert_true(conf.IsInteger("num_subdomains_x"));
    assert_true(conf.IsInteger("num_subdomains_y"));
    assert_true(conf.IsInteger("subdomain_x"));
    assert_true(conf.IsInteger("subdomain_y"));
    assert_true(conf.IsInteger("integration_nsteps"));

    // Check the subdomain size: hard-coded value must match the parameter.
    assert_true(conf.asInt("subdomain_x") == Sx) << "subdomain_x mismatch";
    assert_true(conf.asInt("subdomain_y") == Sy) << "subdomain_y mismatch";

    const int nx = conf.asInt("num_subdomains_x") * Sx;     // global X-size
    const int ny = conf.asInt("num_subdomains_y") * Sy;     // global Y-size

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
    // This time step is defined according to stability criteria.
//    const double dt = std::min(dt_base,
//                        std::min( std::min(dx*dx, dy*dy)/(2.0*D + TINY),
//                                  1.0/(std::fabs(max_vx)/dx +
//                                       std::fabs(max_vy)/dy + TINY) ));
std::cout << "A T T E N T I O N: hard-coded dt" << std::endl;
const double dt = dt_base;      // XXX
    assert_true(dt > TINY);
    conf.SetDouble("dt", dt);
    conf.SetInt("Nt", static_cast<int>(
                        std::ceil(conf.asDouble("integration_period") / dt)));
}

//-----------------------------------------------------------------------------
// Function computes the initial process model covariance matrix P
// based on exponential distance.
//-----------------------------------------------------------------------------
void InitialCovar(const Configuration & conf, SubDomain * sd)
{
    // Subdomains without observations have empty Kalman matrices.
    if (sd->m_P.Size() == 0) {
        assert_true(sd->m_sensors.empty());
        return;
    } else {
        assert_true(!sd->m_sensors.empty());
    }

    // Correlation distances are given in logical coordinates of nodal points.
    const double variance = conf.asDouble("model_ini_var");
    const double covar_radius = conf.asDouble("model_ini_covar_radius");
    const double sigma_x = std::max(covar_radius, 1.0);
    const double sigma_y = std::max(covar_radius, 1.0);
    const int    Rx = static_cast<int>(std::ceil(4.0 * sigma_x));
    const int    Ry = static_cast<int>(std::ceil(4.0 * sigma_y));
    const int    Sx = conf.asInt("subdomain_x");
    const int    Sy = conf.asInt("subdomain_y");

    assert_true(sd->m_P.IsSquare());
    assert_true(sd->m_P.NRows() == (Sx + 2) * (Sy + 2));

    // Mind the extended subdomain: one extra point layer on either side.
    Fill(sd->m_P, 0.0);
    for (int u = 0; u < Sx + 2; ++u) {
    for (int v = 0; v < Sy + 2; ++v) {
        int i = (int)sd->sub2ind_ex(u, v);
        double dx = 0.0, dy = 0.0;
        for (int x = u-Rx; x <= u+Rx; ++x) { if ((0 <= x) && (x < Sx + 2)) {
        for (int y = v-Ry; y <= v+Ry; ++y) { if ((0 <= y) && (y < Sy + 2)) {
            int j = (int)sd->sub2ind_ex(x, y);
            if (i <= j) {
                dx = (u - x) / sigma_x;
                dy = (v - y) / sigma_y;
                sd->m_P(i,j) = sd->m_P(j,i) =
                        variance * std::exp(-0.5 * (dx*dx + dy*dy));
            }
        }}}}
    }}
}

//-----------------------------------------------------------------------------
// Function computes flow components given a discrete time.
// @return a pair of flow components (flow_x, flow_y).
//-----------------------------------------------------------------------------
flow_t Flow(const Configuration & conf, const size_t discrete_time)
{
    const double max_vx = conf.asDouble("flow_model_max_vx");
    const double max_vy = conf.asDouble("flow_model_max_vy");
    const double t = static_cast<double>(discrete_time) / conf.asDouble("Nt");
    return flow_t( -max_vx * std::sin(0.1 * t - M_PI),
                   -max_vy * std::sin(0.2 * t - M_PI) );
}

//-----------------------------------------------------------------------------
// Function applies Dirichlet zero boundary condition to those subdomains that
// are located on the outer border of the whole domain. Note, we set to zero
// the boundary points of the normal subdomain along with the extra points of
// the extended one.
//-----------------------------------------------------------------------------
void ApplyBoundaryCondition(Matrix          & state,
                            const size2d_t  & subdom_extended_size,
                            const point2d_t & subdom_pos,
                            const size2d_t  & grid_size)
{
    assert_true(state.NRows() == subdom_extended_size.x);
    assert_true(state.NCols() == subdom_extended_size.y);

    const int nx = static_cast<int>(subdom_extended_size.x);
    const int ny = static_cast<int>(subdom_extended_size.y);

    if (subdom_pos.x == 0) {
        for (int y = 0; y < ny; ++y) {
            state(0, y) = 0.0;          // extended subdomain points
            state(1, y) = 0.0;          // border points of normal subdomain
        }
    }
    if (subdom_pos.x == grid_size.x - 1) {
        for (int y = 0; y < ny; ++y) {
            state(nx - 2, y) = 0.0;     // border points of normal subdomain
            state(nx - 1, y) = 0.0;     // extended subdomain points
        }
    }

    if (subdom_pos.y == 0) {
        for (int x = 0; x < nx; ++x) {
            state(x, 0) = 0.0;          // extended subdomain points
            state(x, 1) = 0.0;          // border points of normal subdomain
        }
    }
    if (subdom_pos.y == grid_size.y - 1) {
        for (int x = 0; x < nx; ++x) {
            state(x, ny - 2) = 0.0;     // border points of normal subdomain
            state(x, ny - 1) = 0.0;     // extended subdomain points
        }
    }
}

//-----------------------------------------------------------------------------
// Function ensures non-negative (physically plausible) density.
//-----------------------------------------------------------------------------
void ValidateField(Matrix & field)
{
    std::for_each(field.begin(), field.end(),
                  [](double & v){ if (v < 0.0) v = 0.0; });
}

//-----------------------------------------------------------------------------
// Function initializes inverse matrix of implicit Euler time-integrator:
// B * x_{t+1} = x_{t}, where B = A^{-1} is the matrix returned by this
// function. The matrix must be inverted while iterating forward in time:
// x_{t+1} = A * x_{t} = B^{-1} * x_{t}.
// Note, the matrix we generate here is acting on a subdomain.
// Note, the model matrix B is supposed to be a sparse one. For now, since we
// do not have a fast utility for sparse matrix inversion, we define B as
// a dense one with many zeros.
//-----------------------------------------------------------------------------
void InverseModelMatrix(Matrix & B, const Configuration & conf,
                        const flow_t & flow, const size2d_t & subdomain_size,
                        int Nsubiter = 0)
{
    const long Sx = subdomain_size.x;
    const long Sy = subdomain_size.y;

    assert_true((B.NRows() == B.NCols()) && (B.NCols() == (Sx + 2)*(Sy + 2)));

    const double D  = conf.asDouble("diffusion_coef");
    const double dx = conf.asDouble("dx");
    const double dy = conf.asDouble("dy");
    const double dt = conf.asDouble("dt") / std::max(Nsubiter,1);

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
    for (int x = 1; x <= Sx; ++x) {
    for (int y = 1; y <= Sy; ++y) {
        int i = (int) base_sub2ind(x, y, Sx + 2, Sy + 2);
        B(i, i) = 1.0 + 2 * (rho_x + rho_y);
        B(i, (int) base_sub2ind(x - 1, y, Sx + 2, Sy + 2)) = -vx - rho_x;
        B(i, (int) base_sub2ind(x + 1, y, Sx + 2, Sy + 2)) = +vx - rho_x;
        B(i, (int) base_sub2ind(x, y - 1, Sx + 2, Sy + 2)) = -vy - rho_y;
        B(i, (int) base_sub2ind(x, y + 1, Sx + 2, Sy + 2)) = +vy - rho_y;
    }}
}

//-----------------------------------------------------------------------------
// Function computes the process model noise covariance matrix.
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function computes the measurement noise covariance matrix.
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function computes: z = H * observations(t). Since H is a simple 0/1 matrix
// that just picks up the observations at sensor locations, instead of
// matrix-vector multiplication we get the observations directly.
//-----------------------------------------------------------------------------
void GetObservations(SubDomain * sd, long timestamp)
{
    const index_t n = sd->m_observations.NCols();
    assert_true(sd->m_z.Size() == n);
    for (index_t i = 0; i < n; ++i) {
        sd->m_z(i) = sd->m_observations((index_t)timestamp, i);
    }

#ifdef AMDADOS_DEBUGGING
    #warning "Some extra validity test"
    // Mind the extended subdomain: one extra point layer on either side.
    assert_true(sd->m_H.NCols() == (sd->m_size.x + 2) * (sd->m_size.y + 2));
    Matrix subfield(static_cast<int>(sd->m_size.x + 2),
                    static_cast<int>(sd->m_size.y + 2));
    for (index_t i = 0; i < static_cast<index_t>(sd->m_sensors.size()); ++i) {
        subfield(static_cast<int>(sd->m_sensors[(size_t)i].x + 1),
                 static_cast<int>(sd->m_sensors[(size_t)i].y + 1)) =
            sd->m_observations((index_t)timestamp, i);
    }
    Vector _z(n);
    MatVecMult(_z, sd->m_H, subfield);    // _z = H * observations(t)
    assert_true(std::equal(_z.begin(), _z.end(),
                           sd->m_z.begin(), sd->m_z.end()));
#endif
}

//-----------------------------------------------------------------------------
// Handy function prints progress if AMDADOS_DEBUGGING macro is defined.
//-----------------------------------------------------------------------------
void PrintProgress(long t, long Nt)
{
#ifdef AMDADOS_DEBUGGING    // printing progress
    if ((GetRank() == 0) && (t < Nt)) {
        if (t == 0) std::cout << "Nt: " << Nt << std::endl;
        std::cout << "\rtimestep: " << t << std::flush;
        if (t + 1 == Nt) std::cout << std::endl << std::flush;
    }
#else
    (void)t; (void)Nt;
#endif
}

//-----------------------------------------------------------------------------
// Function does a single iteration for each sub-domain without sensors
// during the time integration.
//-----------------------------------------------------------------------------
void SubdomainRoutineNoSensors(const Configuration & conf, SubDomain * sd,
                               long timestamp, long sub_iter)
{
    (void)sub_iter;

    // Compute flow velocity vector.
    flow_t flow = Flow(conf, static_cast<size_t>(timestamp));

    // Construct (inverse) model matrix.
    InverseModelMatrix(sd->m_B, conf, flow, sd->m_size,
                        static_cast<int>(sd->m_Nsubiter));

    // Decompose: B = L*U.
    sd->m_LU.Init(sd->m_B);

    // Propagate state: next_field = B^{-1} * current_field.
    sd->m_LU.Solve(sd->m_next_field, sd->m_curr_field);

    // Ensure boundary conditions on the outer border.
    ApplyBoundaryCondition(sd->m_next_field, sd->m_ex_size,
                           sd->m_pos, sd->m_grid_size);

    // Physically plausible density field must be non-negative everywhere.
    ValidateField(sd->m_next_field);
}

//-----------------------------------------------------------------------------
// Function does a single iteration for each sub-domain, which contains at
// least one sensor, during the time integration. For such a subdomain the
// Kalman filter governs the simulation by pulling the solution towards
// the ground-truth observed at sensor locations.
//-----------------------------------------------------------------------------
void SubdomainRoutineKalman(const Configuration & conf, SubDomain * sd,
                            long timestamp, long sub_iter)
{
    // Compute flow velocity vector.
    flow_t flow = Flow(conf, static_cast<size_t>(timestamp));

    // Start from the current field.
    sd->m_next_field = sd->m_curr_field;

    // At the beginning of a regular iteration (i.e. at the first
    // sub-iteration) do: (1) get the new observations; (2) obtain new noise
    // covariance matrices; (3) compute the prior state estimation.
    if (sub_iter == 0) {
        // Get the current sensor measurements (into sd->m_z vector).
        GetObservations(sd, timestamp);

        // Covariance matrices can change over time.
        ComputeQ(conf, sd->m_Q);
        ComputeR(conf, sd->m_R);

        // Prior estimation.
        InverseModelMatrix(sd->m_B, conf, flow, sd->m_size, 0);
        sd->m_Kalman.PropagateStateInverse(sd->m_next_field,
                                           sd->m_P, sd->m_B, sd->m_Q);
    }

    // Filtering by Kalman filter (on every sub-iteration!).
    sd->m_Kalman.SolveFilter(sd->m_next_field,
                             sd->m_P, sd->m_H, sd->m_R, sd->m_z);

    // Ensure boundary conditions on the outer border.
    ApplyBoundaryCondition(sd->m_next_field, sd->m_ex_size,
                           sd->m_pos, sd->m_grid_size);

    // Physically plausible density field must be non-negative everywhere.
    ValidateField(sd->m_next_field);
}

//-----------------------------------------------------------------------------
// Application's main function runs simulation with data assimilation.
//-----------------------------------------------------------------------------
void RunSimulation(const Configuration & conf)
{
    // Initialize grid and create subdomains attached to this process.
    MpiGrid grid(conf);
    grid.forAll([&grid, &conf](point2d_t pos, int rank) {
        if (rank == grid.myRank()) {
            grid.setSubdomain(pos, new SubDomain(conf, grid, rank, pos));
        }
    });

    // Load sensor locations, sensor observations and finalize
    // Initialization of subdomains attached to this process.
    InputData().Load(conf, grid);
    grid.forAllLocal([&conf](SubDomain * sd) {
        InitialCovar(conf, sd);
    });

    const long Nt = gTestBoundExchange ? 100 : conf.asInt("Nt");
    const long Nsubiter = conf.asInt("num_sub_iter");
    const long Nwrite = std::min(Nt, (long)conf.asInt("write_num_fields"));

    MpiOutputWriter writer(conf);

    // Initialization is done by this line.
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    // Time integration.
    for (long t = 0; t < Nt; ++t) {
        // Write a snapshot (of entire field) in the file.
        if ((t == 0) || ( ((Nwrite - 1) * (t - 1)) / (Nt - 1) !=
                          ((Nwrite - 1) * (t + 0)) / (Nt - 1) ) ) {
            writer.Write(t, grid);
        }

        // Few sub-iterations iron out discrepancies along subdomain boundaries.
        for (long sub_iter = 0; sub_iter < Nsubiter; ++sub_iter) {
            const long timestamp = t * Nsubiter + sub_iter;
            // Initiate data exchange.
            grid.forAllLocal([&](SubDomain * sd) {
                sd->SendBoundariesToNeighbours(grid, timestamp);
            });
            // Receive data for all the subdomains attached to this process.
            grid.forAllLocal([&](SubDomain * sd) {
                sd->ReceiveBoundariesFromNeighbours(timestamp);
            });
            // Upon arrival of the peer data, we can do the processing.
            grid.forAllLocal([&](SubDomain * sd) {
                if (gTestBoundExchange) {
                    TestExchangeCorrectness(sd, timestamp);
                } else {
                    if (sd->m_sensors.empty()) {
                        SubdomainRoutineNoSensors(conf, sd, t, sub_iter);
                    } else {
                        SubdomainRoutineKalman(conf, sd, t, sub_iter);
                    }
                }
            });
            // Wait until remote subdomains confirmed data arrival.
            grid.forAllLocal([](SubDomain * sd) {
                sd->WaitForExchangeCompletion();
            });
            // Replace the current state by the newly computed one.
            grid.forAllLocal([](SubDomain * sd) {
                sd->m_curr_field.swap(sd->m_next_field);
            });
        }
        PrintProgress(t, Nt);
    }

    // Simulation is done by this line.
    MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD));

    // Delete all subdomains attached to this process.
    grid.forAll([&grid](point2d_t pos, int rank) {
        if (rank == grid.myRank()) {
            SubDomain * sd = grid.getSubdomain(pos);
            delete sd;
        }
    });
}

//-----------------------------------------------------------------------------
// Function open a log-file per process rank.
//-----------------------------------------------------------------------------
bool OpenLogFile()
{
    int rank = -1;
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        std::cerr << "ERROR: MPI_Comm_rank() failed" << std::endl;
        return false;
    }
    std::stringstream fname;
    fname << "rank" << rank << ".log";
    gLogFile.open(fname.str(), std::ios::out | std::ios::trunc);
    if (!gLogFile.good()) {
        std::cerr << "ERROR: failed to open log-file" << std::endl;
        return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
// Function tests the correctness of boundary data exchange.
//-----------------------------------------------------------------------------
void TestExchangeCorrectness(SubDomain * sd, long timestamp)
{
    const double EPS = 5.0 * std::numeric_limits<double>::epsilon();
    const int Sx = (int)sd->m_size.x;
    const int Sy = (int)sd->m_size.y;
    const long domain_x = sd->m_grid_size.x * Sx;
    const long domain_y = sd->m_grid_size.y * Sy;
    const int rank = GetRank(); (void)rank;         // useful for debugging

//if (GetRank() == 0) {
//    std::cout << "pos: (" << sd->m_pos.x << "," << sd->m_pos.y << ")" << std::endl;
//    if (sd->m_pos.x == 0 && sd->m_pos.y <= 1) {
//        fprintf(stdout, "+");
//    }
//}

    // Function generates expected unique value given global point coordinates
    // and the timestamp.
    auto Value = [=](long x, long y, long t) -> double {
        if (t < 0) return 0.0;      // once at the very beginning
        x = x + sd->m_pos.x * Sx;   // x <-- global x
        y = y + sd->m_pos.y * Sy;   // y <-- global y
        return static_cast<double>(t + base_sub2ind(x, y, domain_x, domain_y));
    };

    // ***** STAGE 1: check boundary values were correctly exchanged for
    //                current field.

    // Check left-most points from the right boundary of the left neighbour.
    if (sd->m_pos.x > 0) {
        for (int y = 0; y < Sy; ++y) {
            double v1 = sd->m_curr_field(0, y + 1);
            double v2 = Value(-1, y, timestamp - 1);
            bool ok = (std::fabs(v1 - v2) < EPS);
            if (!ok) {
                fprintf(stderr, "e");   // put breakpoint here
            }
            assert_true(ok);
        }
    }
    // Check right-most points from the left boundary of the right neighbour.
    if (sd->m_pos.x + 1 < sd->m_grid_size.x) {
        for (int y = 0; y < Sy; ++y) {
            double v1 = sd->m_curr_field(Sx + 1, y + 1);
            double v2 = Value(Sx, y, timestamp - 1);
            bool ok = (std::fabs(v1 - v2) < EPS);
            if (!ok) {
                fprintf(stderr, "e");   // put breakpoint here
            }
            assert_true(ok);
        }
    }
    // Check bottom-most points from the top boundary of the bottom neighbour.
    if (sd->m_pos.y > 0) {
        for (int x = 0; x < Sx; ++x) {
            double v1 = sd->m_curr_field(x + 1, 0);
            double v2 = Value(x, -1, timestamp - 1);
            bool ok = (std::fabs(v1 - v2) < EPS);
            if (!ok) {
                fprintf(stderr, "e");   // put breakpoint here
            }
            assert_true(ok);

//if (!ok) {
//    for (int u = 0; u < Sx+2; ++u) {
//        for (int v = 0; v < Sy+2; ++v) {
//            std::cout << std::setw(6) << sd->m_curr_field(u,v);
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//    for (int u = 0; u < Sx; ++u) {
//        std::cout << std::setw(6) << Value(u, -1, timestamp - 1);
//    }
//    std::cout << std::endl;
//    MPI_Abort(MPI_COMM_WORLD, 1);
//}
//            assert_true(ok);
//std::cout << "v1=" << v1 << ", v2=" << v2 << std::endl;

        }
    }
    // Check top-most points from the bottom boundary of the top neighbour.
    if (sd->m_pos.y + 1 < sd->m_grid_size.y) {
        for (int x = 0; x < Sx; ++x) {
            double v1 = sd->m_curr_field(x + 1, Sy + 1);
            double v2 = Value(x, Sy, timestamp - 1);
            bool ok = (std::fabs(v1 - v2) < EPS);
            if (!ok) {
                fprintf(stderr, "e");   // put breakpoint here
            }
            assert_true(ok);
        }
    }

    // ***** STAGE 2: create a new field for the next iteration.

    Fill(sd->m_next_field, -777.0);
    for (int x = 0; x < Sx; ++x) {
    for (int y = 0; y < Sy; ++y) {
        sd->m_next_field(x + 1, y + 1) = Value(x, y, timestamp);
    }}

//    if (GetRank() == 0) std::cout << "." << std::flush;
}

}   // namespace amdados

//-----------------------------------------------------------------------------
// Application's main function runs simulation with data assimilation.
//-----------------------------------------------------------------------------
int mpi_main(int * argc, char *** argv)
{
    int ret_val = EXIT_FAILURE;

    // Very first MPI initialization steps.
    if (MPI_Init(argc, argv) != MPI_SUCCESS) {
        std::cerr << "ERROR: MPI_Init() failed" << std::endl;
        return EXIT_FAILURE;
    }

    // Open log-file for this process rank.
//    if (!amdados::OpenLogFile()) {
//        MPI_Finalize();
//        return EXIT_FAILURE;
//    }

    try {
        MY_LOG(INFO) << "***** MPI Amdados2D application *****";

        // Get command-line options.
        std::string scenario, config_file;
        if (amdados::ParseCommandArgs(*argc, *argv, scenario, config_file)) {
            gLogFile.flush();
            gLogFile.close();
            return MPI_Finalize();  // print help and exit
        }

//// XXX
//config_file = "output/tmp.conf";
//gTestBoundExchange = true;

        // Read application parameters from configuration file.
        amdados::Configuration conf;
        conf.ReadConfigFile(config_file);
        amdados::InitDependentParams(conf);
        conf.PrintParameters();

        // Compute new sensor locations and exit, or proceed with simulation.
        if (scenario == "sensors") {
            MY_LOG(INFO) << "SCENARIO: generating random sensors.";
            if (amdados::GetRank() == 0) {
                amdados::ScenarioSensors(conf);
            }
            gLogFile.flush();
            gLogFile.close();
            return MPI_Finalize();
        }
        MY_LOG(INFO) << "SCENARIO: simulation with MPI framework.";

        const auto Nx = conf.asInt("num_subdomains_x") ;
        const auto Ny = conf.asInt("num_subdomains_y") ;
        int steps = static_cast<int>(conf.asUInt("Nt"));

        if (amdados::GetRank() == 0) {
        	// print some status info for the user
        	std::cout << "Running MPI simulation on configuration file \"" << config_file;
        	std::cout << "\" with domain size " << Nx << "x" << Ny;
        	std::cout << " for " << steps << " time steps ...\n";
        	std::cout << "Running full scenario simulation ...\n";
        }
    	auto start = std::chrono::high_resolution_clock::now();
        // Now the main simulation.
        amdados::RunSimulation(conf);
        ret_val = EXIT_SUCCESS;
    	auto end = std::chrono::high_resolution_clock::now();

        if (amdados::GetRank() == 0) {
        	auto duration = end - start;

        	// --- summarize performance data ---
        	double time = (std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() / 1000.0f);
        	std::cout << "Simulation took " << time << "s\n";

        	double throughput = (Nx * Ny * steps) / time;
        	std::cout << "Throughput: " << throughput << " sub-domains/s\n";

        }
        if (gTestBoundExchange) {
            MY_LOG(INFO) << "Test for boundary values exchange has passed";
        } else {
            MY_LOG(INFO) << "Simulation has finished normally";
        }
    } catch (const std::domain_error & e) {
        MY_LOG(ERROR) << "domain error: " << e.what();
    } catch (const std::runtime_error & e) {
        MY_LOG(ERROR) << "runtime error: " << e.what();
    } catch (const std::exception & e) {
        MY_LOG(ERROR) << "exception: " << e.what();
    } catch (...) {
        MY_LOG(ERROR) << "Unsupported exception";
    }
    gLogFile.flush();
    gLogFile.close();
    return ((MPI_Finalize() == EXIT_SUCCESS) ? ret_val : EXIT_FAILURE);
}

//-----------------------------------------------------------------------------
// Function prints error message and terminates the MPI application.
//-----------------------------------------------------------------------------
void PrintErrorAndExit(const char * err_msg)
{
    if (err_msg != nullptr) { gLogFile << err_msg << std::endl; }
    gLogFile.flush();
    gLogFile.close();
    std::this_thread::sleep_for(std::chrono::seconds(2));
    MPI_Abort(MPI_COMM_WORLD, 1);
}
