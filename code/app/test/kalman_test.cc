//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#ifndef AMDADOS_PLAIN_MPI
// This implementation is based on Allscale API, no MPI at all.

#include <gtest/gtest.h>
#include "allscale/utils/assert.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include "../include/amdados_utils.h"
#include "../include/matrix.h"
#include "../include/cholesky.h"
#include "../include/lu.h"
#include "../include/kalman_filter.h"

//-----------------------------------------------------------------------------
// Function for testing a Kalman filter on a toy 2D problem.
//-----------------------------------------------------------------------------
TEST(KalmanFilter, Basic)
{
    using namespace ::amdados;

    // Open the output log-file.
    std::string fname = "kalman_test.log";
    std::fstream log_file(fname, std::ios::out | std::ios::trunc);
    assert_true(log_file.good()) << "failed to open: " << fname << std::endl;
    fname = "kalman_test.out";
    std::fstream track_file(fname, std::ios::out | std::ios::trunc);
    assert_true(track_file.good()) << "failed to open: " << fname << std::endl;

    const int    DIM = 3;
    const int    NUM_TIME_STEPS = 5000;
    const double MODEL_STD = 1.0;
    const double MEASUREMENT_STD = 1.0;
    const char   SPACE[] = " ";

    // Gaussian noise generator.
    std::mt19937_64            gen(RandomSeed());
    std::normal_distribution<> distrib_v(0.0, MEASUREMENT_STD);

    Vector x(DIM);      Fill(x, 0.0);
    Matrix P(DIM,DIM);  MakeIdentityMatrix(P);  ScalarMult(P, MODEL_STD);
    Matrix A(DIM,DIM);  MakeIdentityMatrix(A);
    Matrix Q(DIM,DIM);  MakeIdentityMatrix(Q);  ScalarMult(Q, MODEL_STD);
    Matrix H(DIM,DIM);  MakeIdentityMatrix(H);
    Matrix R(DIM,DIM);  MakeIdentityMatrix(R);  ScalarMult(R, MEASUREMENT_STD);

    Vector z(DIM), v(DIM);
    Vector true_x(DIM);
    Vector prev_x = x;
    Vector prev_true_x = x;
    double total_mean_dev = 0.0; // mean deviation estimated over entire period
    double mean_step = 0.0;      // mean step in space made by the particle
    double mean_inno = 0.0;      // mean innovation between successive states

    KalmanFilter kf;

    // Iterate the Kalman filter over time following a fancy spiral path.
    for(int k = 0; k < NUM_TIME_STEPS; ++k) {
        double t = k * (10.0 * M_PI) / (NUM_TIME_STEPS - 1);

        // Generate the next true state.
        double scale = (1 + 0.1*std::cos(5.0*t));
        true_x(0) = t * std::sqrt(std::fabs(t)) * scale * std::cos(t);
        true_x(1) = t * std::sqrt(std::fabs(t)) * scale * std::sin(t);
        true_x(2) = t                           * scale;

        // Some statistics for debugging.
        if (k > 0) { mean_step += NormDiff(true_x, prev_true_x); }
        prev_true_x = true_x;
        mean_inno += NormDiff(x, prev_x);
        prev_x = x;

        // Measurement noise.
        MakeRandom(v, 'n');
        ScalarMult(v, MEASUREMENT_STD);

        // Measurements: z = H*x + v.
        MatVecMult(z, H, true_x);
        AddVectors(z, z, v);

        // Make a Kalman filter iteration. Since A is identity matrix,
        // it does not matter whether we take A or its inverse.
        kf.PropagateStateInverse(x, P, A, Q);
        kf.SolveFilter(x, P, H, R, z);

        // Mean deviation is a square root of the mean diagonal element
        // of the matrix P.
        double mean_dev = Trace(P) / DIM;
        total_mean_dev += mean_dev;
        mean_dev = std::sqrt(std::fabs(mean_dev));

        // Print the current time followed the true state vector followed
        // by the noisy measurements followed by the state vector estimated
        // by the Kalman filter followed by the mean eigenvalue of the
        // estimated covariance matrix.
        track_file << t << SPACE;
        for (int j = 0; j < DIM; ++j) { track_file << true_x(j) << SPACE; }
        for (int j = 0; j < DIM; ++j) { track_file << z(j) << SPACE; }
        for (int j = 0; j < DIM; ++j) { track_file << x(j) << SPACE; }
        track_file << mean_dev << std::endl;
    }
    track_file.flush();
    mean_inno /= (NUM_TIME_STEPS - 1);   // first step was missed
    mean_step /= (NUM_TIME_STEPS - 1);   // first step was missed
    EXPECT_NEAR(mean_step, 0.47272, 1e-5);
    total_mean_dev = std::sqrt(std::fabs(total_mean_dev / NUM_TIME_STEPS));
    EXPECT_NEAR(total_mean_dev, 0.786159, 1e-5);

    log_file << "TestKalmanFilter():" << std::endl;
    log_file << "model noise deviation: " << MODEL_STD << std::endl;
    log_file << "measurement noise deviation: " << MEASUREMENT_STD << std::endl;
    log_file << "mean innovation: " << mean_inno << std::endl;
    log_file << "mean step: " << mean_step << std::endl;
    log_file << "mean estimated deviation: " << total_mean_dev << std::endl;
    log_file << std::endl << std::flush;
}

#endif  // AMDADOS_PLAIN_MPI
