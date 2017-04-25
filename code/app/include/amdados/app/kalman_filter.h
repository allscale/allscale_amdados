#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

//=================================================================================================
// Class implements a Kalman filter using AllScale API grid for matrices and vectors.
//=================================================================================================
template<size_t PROBLEM_SIZE, size_t NUM_MEASUREMENTS>
class KalmanFilter
{
public:
    using matrix_t     = allscale::utils::grid<double, PROBLEM_SIZE, PROBLEM_SIZE>;
    using matrix_MxN_t = allscale::utils::grid<double, NUM_MEASUREMENTS, PROBLEM_SIZE>;
    using matrix_NxM_t = allscale::utils::grid<double, PROBLEM_SIZE, NUM_MEASUREMENTS>;
    using matrix_MxM_t = allscale::utils::grid<double, NUM_MEASUREMENTS, NUM_MEASUREMENTS>;
    using vector_t     = allscale::utils::grid<double, PROBLEM_SIZE>;
    using vector_obs_t = allscale::utils::grid<double, NUM_MEASUREMENTS>;

private:
    Cholesky<PROBLEM_SIZE> m_chol;  // object computes inverse matrix by Cholesky decomposition
    vector_t               m_x;     // vector of state variables
    matrix_t               m_P;     // covariance matrix

    vector_t     m_x_prior;     // placeholder for the vector x_{k|k-1} = A * x
    vector_obs_t m_y;           // placeholder vector of observations
    vector_obs_t m_invSy;       // placeholder vector for S^{-1} * y
    matrix_MxM_t m_S;           // placeholder for the matrix S = H * P_{k|k-1} * H^t + R
    matrix_t     m_P_prior;     // placeholder for the matrix P_{k|k-1}
    matrix_t     m_PAt;         // placeholder for the matrix P * A^t
    matrix_NxM_t m_PHt;         // placeholder for the matrix P_{k|k-1} * H^t
    matrix_MxN_t m_HP;          // placeholder for the matrix H * P_{k|k-1}
    matrix_MxN_t m_invSHP;      // placeholder for the matrix S^{-1} * H * P_{k|k-1}

public:
//-------------------------------------------------------------------------------------------------
// Constructor initializes the vector of state variables and the covariance matrix
// by the original values at the very first timestamp.
//-------------------------------------------------------------------------------------------------
KalmanFilter(const vector_t & x0, const matrix_t & P0) : m_x(x0), m_P(P0)
{
}

//-------------------------------------------------------------------------------------------------
// Function returns the current vector of state variables "x".
//-------------------------------------------------------------------------------------------------
const vector_t & GetStateVector() const
{
    return m_x;
}

//-------------------------------------------------------------------------------------------------
// Function returns the current covariance matrix "P".
//-------------------------------------------------------------------------------------------------
const matrix_t & GetCovariance() const
{
    return m_P;
}

//-------------------------------------------------------------------------------------------------
// Function makes an iteration of Kalman filter which includes prediction and correction phases.
// \param  A  model transition matrix: x_k = A_k * x_{k-1} + w_{k-1}.
// \param  Q  process noise (w_k) covariance.
// \param  H  observation model: z_k = H_k * x_k + v_k.
// \param  R  measurement noise (v_k) covariance.
// \param  z  vector of observations.
//-------------------------------------------------------------------------------------------------
void Iterate(const matrix_t & A, const matrix_t & Q,
             const matrix_MxN_t & H, const matrix_MxM_t & R, const vector_obs_t & z)
{
    // x_prior = A * x
    MatVecMult(m_x_prior, A, m_x);

    // P_prior = A * P * A^t + Q
    MatMultTransposed(m_PAt, m_P, A);
    MatMult(m_P_prior, A, m_PAt);
    AddMatrices(m_P_prior, m_P_prior, Q);

    // y = z - H * x_prior
    MatVecMult(m_y, H, m_x_prior);
    SubtractVectors(m_y, z, m_y);

    // S = H * P_prior * H^t + R
    MatMultTransposed(m_PHt, m_P_prior, H);
    MatMult(m_S, H, m_PHt);
    AddMatrices(m_S, m_S, R);

    // Correct symmetry loss due to round-off errors.
    Symmetrize(m_S);

    // Compute Cholesky decomposition  S = L * L^t  to facilitate matrix inversion.
    m_chol.ComputeDecomposition(m_S);

    // m_invSy = S^{-1} * y
    m_chol.Solve(m_invSy, m_y);

    // x  =  x_prior + K * y  =  x_prior + P_prior * H^t * S^{-1} * y
    MatVecMult(m_x, m_PHt, m_invSy);
    AddVectors(m_x, m_x, m_x_prior);

    // m_invSHP = S^{-1} * H * P_prior
    GetTransposed(m_HP, m_PHt);
    m_chol.BatchSolve(m_invSHP, m_HP);

    // P  =  (I - K * H) * P_prior  =  P_prior - P_prior * H^t * S^{-1} * H * P_prior.
    MatMult(m_P, m_PHt, m_invSHP);
    SubtractMatrices(m_P, m_P_prior, m_P);

    // Correct symmetry loss due to round-off errors.
    Symmetrize(m_P);
}

}; // class KalmanFilter

//-------------------------------------------------------------------------------------------------
// Function for testing a Kalman filter on a toy 2D problem.
//-------------------------------------------------------------------------------------------------
void TestKalmanFilter()
{
    std::cout << "Running TestKalmanFilter() ..." << std::endl;

    const size_t DIM = 3;
    const int    NUM_TIME_STEPS = 5000;
    const double DEVIATION_MODEL = 1.0;
    const double DEVIATION_MEASUREMENT = 1.0;
    const char   SPACE[] = " ";

    using matrix_t     = KalmanFilter<DIM,DIM>::matrix_t;
    using matrix_MxN_t = KalmanFilter<DIM,DIM>::matrix_MxN_t;
    using matrix_MxM_t = KalmanFilter<DIM,DIM>::matrix_MxM_t;
    using vector_t     = KalmanFilter<DIM,DIM>::vector_t;
    using vector_obs_t = KalmanFilter<DIM,DIM>::vector_obs_t;

    // Gaussian noise generator.
    std::random_device         rd;
    std::mt19937               gen(rd());
    std::normal_distribution<> distrib_v(0.0, DEVIATION_MEASUREMENT);

    vector_t     x0;    FillVector(x0, 0.0);
    matrix_t     P0;    MakeIdentityMatrix(P0);     MatScalarMult(P0, DEVIATION_MODEL);
    matrix_t     A;     MakeIdentityMatrix(A);
    matrix_t     Q;     MakeIdentityMatrix(Q);      MatScalarMult(Q, DEVIATION_MODEL);
    matrix_MxN_t H;     MakeIdentityMatrix(H);
    matrix_MxM_t R;     MakeIdentityMatrix(R);      MatScalarMult(R, DEVIATION_MEASUREMENT);
    vector_t     x;
    vector_obs_t z, v;
    vector_t     prev_x;
    double       total_mean_dev = 0.0;  // mean estimated deviation over the whole time period
    double       mean_step = 0.0;       // mean step in space made by the particle

    // Initialize the Kalman filter.
    KalmanFilter<DIM,DIM> kf(x0, P0);

    // Open a file for storing the results of particle path tracking by the Kalman filter.
    std::cout.flush();
    std::fstream file;
    std::string filename = "/tmp/kalman_filter_test.out";   // on Linux or MacOS
    file.exceptions(std::fstream::failbit | std::fstream::badbit);
    try {
        file.open(filename, std::ios::trunc | std::ios::out);
    } catch (std::fstream::failure e) {
        filename = "kalman_filter_test.out";                // on Windows
        file.open(filename, std::ios::trunc | std::ios::out);
    }

    // Iterate the Kalman filter over time following a fancy spiral path.
    prev_x = x0;
    for (int k = 0; k < NUM_TIME_STEPS; ++k) {
        double t = k * (10.0 * M_PI) / (NUM_TIME_STEPS - 1);

        // Generate the next true state. Important: x == x0 at the beginning.
        x[{0}] = t * std::sqrt(std::fabs(t)) * (1 + 0.1*std::cos(5.0*t)) * std::cos(t);
        x[{1}] = t * std::sqrt(std::fabs(t)) * (1 + 0.1*std::cos(5.0*t)) * std::sin(t);
        x[{2}] = t                           * (1 + 0.1*std::cos(5.0*t));
        if (k == 0) assert(NormVecDiff(x, x0) < 3*std::numeric_limits<double>::epsilon());
        if (k > 0) mean_step += NormVecDiff(x, prev_x);

        // Measurement noise.
        for (int j = 0; j < static_cast<int>(DIM); ++j) {
            v[{j}] = distrib_v(gen);
        }

        // Measurements: z = H*x + v.
        MatVecMult(z, H, x);
        AddVectors(z, z, v);

        // Make a Kalman filter iteration.
        kf.Iterate(A, Q, H, R, z);
        const vector_t & x_est = kf.GetStateVector();
        const matrix_t & P_est = kf.GetCovariance();

        // Mean deviation is a square root of the mean diagonal element of the matrix P_est.
        double mean_dev = 0.0;
        for (int j = 0; j < static_cast<int>(DIM); ++j) mean_dev += P_est[{j,j}];
        mean_dev /= static_cast<double>(DIM);
        total_mean_dev += mean_dev;
        mean_dev = std::sqrt(std::fabs(mean_dev));

        // Print the current time followed the true state vector followed by the noisy
        // measurements followed by the state vector estimated by the Kalman filter
        // followed by the mean eigenvalue of the estimated covariance matrix.
        file << t << SPACE;
        for (int j = 0; j < static_cast<int>(DIM); ++j) { file <<     x[{j}] << SPACE; }
        for (int j = 0; j < static_cast<int>(DIM); ++j) { file <<     z[{j}] << SPACE; }
        for (int j = 0; j < static_cast<int>(DIM); ++j) { file << x_est[{j}] << SPACE; }
        file << mean_dev << std::endl;

        prev_x = x;
    }
    file.flush();
    mean_step /= static_cast<double>(NUM_TIME_STEPS - 1);   // first step was missed
    total_mean_dev = std::sqrt(std::fabs(total_mean_dev / static_cast<double>(NUM_TIME_STEPS)));
    std::cout << "model noise deviation: " << DEVIATION_MODEL << std::endl;
    std::cout << "measurement noise deviation: " << DEVIATION_MEASUREMENT << std::endl;
    std::cout << "mean step: " << mean_step << std::endl;
    std::cout << "mean estimated deviation: " << total_mean_dev << std::endl;
    std::cout << "TestKalmanFilter() is done" << std::endl << std::endl;
}

