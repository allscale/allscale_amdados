//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Class implements a Kalman filter.
//=============================================================================
class KalmanFilter
{
private:
    Cholesky        m_chol;     // Cholesky decomposition solver
    LUdecomposition m_lu;       // LU decomposition solver

    Vector m_x_tmp;     // placeholder for the vector x_{k|k-1} = A*x
    Vector m_y;         // placeholder vector of observations
    Vector m_invSy;     // placeholder vector for S^{-1}*y
    Matrix m_S;         // placeholder for the matrix S = H*P_{k|k-1}*H^t + R
    Matrix m_P_tmp;     // placeholder for the matrix P_{k|k-1}
    Matrix m_PHt;       // placeholder for the matrix P_{k|k-1}*H^t
    Matrix m_HP;        // placeholder for the matrix H*P_{k|k-1}
    Matrix m_invSHP;    // placeholder for the matrix S^{-1}*H*P_{k|k-1}

public:
//-----------------------------------------------------------------------------
// Constructor initializes all vector/matrix variables by zeros.
//-----------------------------------------------------------------------------
KalmanFilter()
{
}

friend std::ostream& operator<<(std::ostream& out, const KalmanFilter& kf) {
	out << "KalmanFilter: [ ";
	out << kf.m_chol << ", ";
	out << kf.m_lu << ", ";
	out << kf.m_x_tmp << ", ";
	out << kf.m_y << ", ";
	out << kf.m_invSy << ", ";
	out << kf.m_S << ", ";
	out << kf.m_P_tmp << ", ";
	out << kf.m_PHt << ", ";
	out << kf.m_HP << ", ";
	out << kf.m_invSHP << ", ";
	out << " ]" << std::endl;
	return out;
}

//-----------------------------------------------------------------------------
// Function propagates state and covariance one timestep ahead and obtains
// prior estimations: x_prior = A*x, P_prior = A*P*A^t + Q, where A is
// the process model matrix, which is available by its inversion B = A^{-1}.
// @param  x  in: current state; out: prior state estimation.
// @param  P  in: current covariance; out: prior state covariance estimation.
// @param  B  inverse model matrix: B = A^{-1}.
// @param  Q  process noise covariance.
//-----------------------------------------------------------------------------
void PropagateStateInverse(VectorView & x, Matrix & P,
                           const Matrix & B, const Matrix & Q)
{
    const int N = x.Size();     // problem size

    assert_true((P.NRows() == N) && (P.NCols() == N));
    assert_true(B.SameSize(P) && Q.SameSize(P));

    m_x_tmp = x;                    // copy state and covariance into
    m_P_tmp = P;                    // the separate temporary objects

    m_lu.Init(B);                   // decompose: B = L*U
    m_lu.Solve(x, m_x_tmp);         // x_prior = B^{-1}*x

    m_lu.BatchSolve(m_P_tmp, P);    // P_tmp = B^{-1}*P, where P is symmetric
    m_lu.BatchSolveTr(P, m_P_tmp);  // P_prior = B^{-1}*(B^{-1}*P)^t = A*P*A^t

    AddMatrices(P, P, Q);           // P_prior = A*P*A^t + Q
    Symmetrize(P);                  // correct the loss of symmetry
}

//-----------------------------------------------------------------------------
// Function makes an iteration of Kalman filter given already estimated
// (prior) state and its covariance.
// @param  x  in: prior state estimation;
//            out: posterior state estimation.
// @param  P  in: prior state covariance estimation;
//            out: posterior state covariance estimation.
// @param  H  observation model: z = H*x + v.
// @param  R  measurement noise (v) covariance.
// @param  z  vector of observations.
//-----------------------------------------------------------------------------
void SolveFilter(VectorView & x, Matrix & P,
                 const Matrix & H, const Matrix & R, const VectorView & z)
{
    const int N = x.Size();     // problem size
    const int O = z.Size();     // number of observations

    assert_true((P.NRows() == N) && (P.NCols() == N));
    assert_true((H.NRows() == O) && (H.NCols() == N));
    assert_true((R.NRows() == O) && (R.NCols() == O));

    // Resize temporary buffer without initialization.
    m_y.Resize(O, false);
    m_invSy.Resize(O, false);
    m_S.Resize(O, O, false);
    m_PHt.Resize(N, O, false);
    m_HP.Resize(O, N, false);
    m_invSHP.Resize(O, N, false);

    m_x_tmp = x;                        // copy state and covariance into
    m_P_tmp = P;                        // the separate temporary objects

    const auto & x_prior = m_x_tmp;     // original state
    const auto & P_prior = m_P_tmp;     // original covariance

    // y = z - H*x_prior
    MatVecMult(m_y, H, x_prior);
    SubtractVectors(m_y, z, m_y);

    // S = H*P_prior*H^t + R
    MatMultTr(m_PHt, P_prior, H);
    MatMult(m_S, H, m_PHt);
    AddMatrices(m_S, m_S, R);

    // Correct symmetry loss due to round-off errors.
    Symmetrize(m_S);

    // Compute Cholesky decomposition S = L*L^t to facilitate matrix inversion.
    m_chol.Init(m_S);

    // m_invSy = S^{-1}*y
    m_chol.Solve(m_invSy, m_y);

    // x = x_prior + K*y = x_prior + P_prior*H^t*S^{-1}*y
    MatVecMult(x, m_PHt, m_invSy);
    AddVectors(x, x, x_prior);

    // m_invSHP = S^{-1}*H*P_prior
    GetTransposed(m_HP, m_PHt);
    m_chol.BatchSolve(m_invSHP, m_HP);

    // P = (I - K*H)*P_prior = P_prior - P_prior*H^t*S^{-1}*H*P_prior.
    MatMult(P, m_PHt, m_invSHP);
    SubtractMatrices(P, P_prior, P);

    // Correct symmetry loss due to round-off errors.
    Symmetrize(P);
}

}; // class KalmanFilter

} // namespace amdados
