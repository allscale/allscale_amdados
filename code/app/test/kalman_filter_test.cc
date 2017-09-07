//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/geometry.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/sparse_matrix.h"
#include "amdados/app/utils/cholesky.h"
#include "amdados/app/utils/lu.h"
#include "amdados/app/utils/configuration.h"
#include "amdados/app/model/i_model.h"
#include "amdados/app/utils/kalman_filter.h"

using namespace amdados::app;
using namespace amdados::app::utils;

//-------------------------------------------------------------------------------------------------
// Function for testing a Kalman filter on a toy 2D problem.
//-------------------------------------------------------------------------------------------------
TEST(KalmanFilter, Basic)
{
    // Read configuration settings.
    Configuration conf;         // a single global configuration
    conf.ReadConfigFile("../../amdados.conf");
    conf.PrintParameters();
    MakeDirectory(conf.asCString("test_output_dir"));

    // Open a file for printing summary of the test.
    std::string summaryFileName = conf.asString("test_output_dir") + "/kalman_summary_test.log";
    std::fstream summaryFile(summaryFileName, std::ios::out | std::ios::trunc);
    assert_true(summaryFile.good())
        << "failed to oped the summary file: " << summaryFileName << endl;

    // Open a file for storing the results of particle path tracking by the Kalman filter.
    std::string trackFileName = conf.asString("test_output_dir") + "/kalman_filter_test.out";
    std::fstream trackFile(trackFileName, std::ios::out | std::ios::trunc);
    assert_true(trackFile.good())
        << "failed to oped the track file: " << trackFileName << endl;

    const int    DIM = 3;
    const int    NUM_TIME_STEPS = 5000;
    const double DEVIATION_MODEL = 1.0;
    const double DEVIATION_MEASUREMENT = 1.0;
    const char   SPACE[] = " ";

    using matrix_t     = KalmanFilter<DIM,DIM>::matrix_t;
    using matrix_OxN_t = KalmanFilter<DIM,DIM>::matrix_OxN_t;
    using matrix_OxO_t = KalmanFilter<DIM,DIM>::matrix_OxO_t;
    using vector_t     = KalmanFilter<DIM,DIM>::vector_t;
    using vector_obs_t = KalmanFilter<DIM,DIM>::vector_obs_t;

    // Gaussian noise generator.
    std::mt19937               gen(std::time(nullptr));
    std::normal_distribution<> distrib_v(0.0, DEVIATION_MEASUREMENT);

    vector_t     x;     FillVector(x, 0.0);
    matrix_t     P;     MakeIdentityMatrix(P);      MatScalarMult(P, DEVIATION_MODEL);
    matrix_t     A;     MakeIdentityMatrix(A);
    matrix_t     Q;     MakeIdentityMatrix(Q);      MatScalarMult(Q, DEVIATION_MODEL);
    matrix_OxN_t H;     MakeIdentityMatrix(H);
    matrix_OxO_t R;     MakeIdentityMatrix(R);      MatScalarMult(R, DEVIATION_MEASUREMENT);
    vector_obs_t z, v;
    vector_t     true_x;
    vector_t     prev_x = x;
    vector_t     prev_true_x = x;
    double       total_mean_dev = 0.0;  // mean deviation estimated over the whole time period
    double       mean_step = 0.0;       // mean step in space made by the particle
    double       mean_inno = 0.0;       // mean innovation between successive states

    KalmanFilter<DIM,DIM> kf;

    // Iterate the Kalman filter over time following a fancy spiral path.
    for(int k = 0; k < NUM_TIME_STEPS; ++k) {
        double t = k * (10.0 * M_PI) / (NUM_TIME_STEPS - 1);

        // Generate the next true state.
        true_x(0) = t * std::sqrt(std::fabs(t)) * (1 + 0.1*std::cos(5.0*t)) * std::cos(t);
        true_x(1) = t * std::sqrt(std::fabs(t)) * (1 + 0.1*std::cos(5.0*t)) * std::sin(t);
        true_x(2) = t                           * (1 + 0.1*std::cos(5.0*t));

        // Some statistics for debugging.
        if (k > 0) { mean_step += NormVecDiff(true_x, prev_true_x); }
        prev_true_x = true_x;
        mean_inno += NormVecDiff(x, prev_x);
        prev_x = x;

        // Measurement noise.
        MakeRandomVector(v, 'n');
        VecScalarMult(v, DEVIATION_MEASUREMENT);

        // Measurements: z = H*x + v.
        MatVecMult(z, H, true_x);
        AddVectors(z, z, v);

        // Make a Kalman filter iteration.
        kf.Iterate(A, Q, H, R, z, x, P);

        // Mean deviation is a square root of the mean diagonal element of the matrix P.
        double mean_dev = Trace(P) / DIM;
        total_mean_dev += mean_dev;
        mean_dev = std::sqrt(std::fabs(mean_dev));

        // Print the current time followed the true state vector followed by the noisy
        // measurements followed by the state vector estimated by the Kalman filter
        // followed by the mean eigenvalue of the estimated covariance matrix.
        trackFile << t << SPACE;
        for (int j = 0; j < DIM; ++j) { trackFile << true_x(j) << SPACE; }
        for (int j = 0; j < DIM; ++j) { trackFile << z(j) << SPACE; }
        for (int j = 0; j < DIM; ++j) { trackFile << x(j) << SPACE; }
        trackFile << mean_dev << endl;
    }
    trackFile.flush();
    mean_inno /= (NUM_TIME_STEPS - 1);   // first step was missed
    mean_step /= (NUM_TIME_STEPS - 1);   // first step was missed
    EXPECT_NEAR(mean_step, 0.47272, 1e-5);
    total_mean_dev = std::sqrt(std::fabs(total_mean_dev / NUM_TIME_STEPS));
    EXPECT_NEAR(total_mean_dev, 0.786159, 1e-5);

    summaryFile << "TestKalmanFilter():" << endl;
    summaryFile << "model noise deviation: " << DEVIATION_MODEL << endl;
    summaryFile << "measurement noise deviation: " << DEVIATION_MEASUREMENT << endl;
    summaryFile << "mean innovation: " << mean_inno << endl;
    summaryFile << "mean step: " << mean_step << endl;
    summaryFile << "mean estimated deviation: " << total_mean_dev << endl;
    summaryFile << endl << flush;
}

