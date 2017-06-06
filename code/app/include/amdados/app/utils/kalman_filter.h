//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once
#include "amdados/app/static_grid.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/cholesky.h"

namespace amdados {
namespace app {
namespace utils {

//=================================================================================================
// Class implements a Kalman filter using AllScale API grid for matrices and vectors.
//=================================================================================================
template<size_t PROBLEM_SIZE, size_t NUM_MEASUREMENTS>
class KalmanFilter
{
public:
    using matrix_t     = Matrix<PROBLEM_SIZE, PROBLEM_SIZE>;
    using matrix_MxN_t = Matrix<NUM_MEASUREMENTS, PROBLEM_SIZE>;
    using matrix_NxM_t = Matrix<PROBLEM_SIZE, NUM_MEASUREMENTS>;
    using matrix_MxM_t = Matrix<NUM_MEASUREMENTS, NUM_MEASUREMENTS>;
    using vector_t     = Vector<PROBLEM_SIZE>;
    using vector_obs_t = Vector<NUM_MEASUREMENTS>;

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

    bool         m_ready;       // true, if initial values had been set to m_x and m_P

public:
//-------------------------------------------------------------------------------------------------
// Default constructor sets all the variables to zero.
//-------------------------------------------------------------------------------------------------
KalmanFilter() : m_chol()
{
    FillVector(m_x);
    FillMatrix(m_P);
    FillVector(m_x_prior);
    FillVector(m_y);
    FillVector(m_invSy);
    FillMatrix(m_S);
    FillMatrix(m_P_prior);
    FillMatrix(m_PAt);
    FillMatrix(m_PHt);
    FillMatrix(m_HP);
    FillMatrix(m_invSHP);
    m_ready = false;
}

//-------------------------------------------------------------------------------------------------
// Constructor initializes the vector of state variables and the covariance matrix
// by the original values at the very first timestamp.
//-------------------------------------------------------------------------------------------------
void Init(const vector_t & x0, const matrix_t & P0)
{
    m_x = x0;
    m_P = P0;
    m_ready = true;
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
    // Estimate posterior state and covariance.
    PosteriorEstimation(H, R, z);
}

//-------------------------------------------------------------------------------------------------
// Function makes an iteration of Kalman filter which includes prediction and correction phases.
// It is assumed that the model adapter can compute the prior estimations of state and covariance,
// e.g. for the linear transition operator A: x_prior = A * x, P_prior = A * P * A^t.
// \param  model  model adapter for computing the prior estimations of state and covariance.
// \param  Q      process noise (w_k) covariance.
// \param  H      observation model: z_k = H_k * x_k + v_k.
// \param  R      measurement noise (v_k) covariance.
// \param  z      vector of observations.
//-------------------------------------------------------------------------------------------------
template<typename ModelAdapter>
void IterateWithModel(ModelAdapter & model, const matrix_t & Q,
                      const matrix_MxN_t & H, const matrix_MxM_t & R, const vector_obs_t & z)
{
    // x_prior = A * x, P_prior = A * P * A^t, where A is some linear operator (matrix here)
    m_x_prior = m_x;
    m_P_prior = m_P;
    model(m_x_prior, m_P_prior);    // model updates the state and covariance
    // P_prior += Q
    AddMatrices(m_P_prior, m_P_prior, Q);
    // Estimate posterior state and covariance.
    PosteriorEstimation(H, R, z);
}

private:
//-------------------------------------------------------------------------------------------------
// Function makes an iteration of Kalman filter given already estimated "x_prior" and "P_prior".
// \param  H  observation model: z_k = H_k * x_k + v_k.
// \param  R  measurement noise (v_k) covariance.
// \param  z  vector of observations.
//-------------------------------------------------------------------------------------------------
void PosteriorEstimation(const matrix_MxN_t & H, const matrix_MxM_t & R, const vector_obs_t & z)
{
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



//>>>>> TODO: temporary function <<<<<
public:
void Iterate(const vector_t & x_prior, const matrix_t & Q,
             const matrix_MxN_t & H, const matrix_MxM_t & R, const vector_obs_t & z)
{
    assert(m_ready);

    // In conventional Kalman filter we would do: x_prior = A * x, but here we just copy.
    m_x_prior = x_prior;

    // In conventional Kalman filter we would do: P_prior = A * P * A^t + Q, but here A = I.
    AddMatrices(m_P_prior, m_P, Q);

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

} // end namespace utils
} // end namespace app
} // end namespace amdados

