//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Class for computing Cholesky decomposition of a symmetric, squared,
// positive-definite matrix. Once decomposition is done in constructor, the
// class instance can be used to solve linear systems A*x = b and A*X = B,
// where x, b are vectors and X, B are matrices respectively.
// The implementation was adopted from "Numerical Recipes" book,
// 3rd edition, chapter 2.9.
//=============================================================================
class Cholesky
{
private:
    Matrix m_L;     // lower triangular matrix of decomposition

public:
//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
Cholesky() : m_L()
{
}

//-----------------------------------------------------------------------------
// Function computes and stores Cholesky decomposition
// of a positive-definite symmetric matrix: A = L * L^t.
//-----------------------------------------------------------------------------
void Init(const Matrix & A)
{
    const double TINY = std::numeric_limits<double>::min() /
		       std::pow(std::numeric_limits<double>::epsilon(),3);

    assert_true(A.IsSquare());
    const int N = A.NRows();    // problem size, A is square

    m_L = A;                // copy the input matrix, then do decomposition
    Matrix & L = m_L;       // short-hand alias

    // Compute the lower triangular matrix of Cholesky decomposition.
    for (int i = 0; i < N; i++) {
    for (int j = i; j < N; j++) {
        if (L(j, i) != L(i, j))
            assert_true(0) << "Cholesky expects a symmetric matrix";

        double sum = L(i,j);
        for (int k = i - 1; k >= 0; k--) { sum -= L(i,k) * L(j,k); }

        if (i == j) {
            assert_gt(sum, TINY) << "Cholesky failed, sum: " << sum;
            L(i,i) = std::sqrt(sum);
        } else {
            L(j,i) = sum / L(i,i);
        }
    }}

    // Put the upper triangular matrix to zero.
    for (int i = 0; i < N; i++) {
    for (int j = 0; j < i; j++) { L(j,i) = 0.0; }}
}

//-----------------------------------------------------------------------------
// Function solves a linear system A*x = b, where A is the matrix whose
// Cholesky decomposition was computed by the Init() function.
//-----------------------------------------------------------------------------
void Solve(Vector & x, const Vector & b) const
{
    const Matrix & L = m_L;         // short-hand alias
    const int      N = L.NRows();   // problem size, L is symmetric

    bool ok = ((x.Size() == N) && (b.Size() == N));
    assert_true(ok);

    for (int i = 0; i < N; i++) {
        double sum = b(i);
        for (int k = i - 1; k >= 0; k--) { sum -= L(i,k) * x(k); }
        x(i) = sum / L(i,i);
    }

    for (int i = N - 1; i >= 0; i--) {
        double sum = x(i);
        for (int k = i + 1; k < N; k++) { sum -= L(k,i) * x(k); }
        x(i) = sum / L(i,i);
    }
}

//-----------------------------------------------------------------------------
// Function solves a collection of linear systems A*X = B, where A is the
// matrix whose Cholesky decomposition was computed by the Init() function,
// X and B are the matrices of the same size.
//-----------------------------------------------------------------------------
void BatchSolve(Matrix & X, const Matrix & B) const
{
    const Matrix & L = m_L;         // short-hand alias
    const int      N = L.NRows();   // problem size, L is square symmetric
    const int      K = X.NCols();   // number of linear systems to solve

    assert_true((N == X.NRows()) && X.SameSize(B));

    for (int c = 0; c < K; c++) {
        for (int i = 0; i < N; i++) {
            double sum = B(i,c);
            for (int k = i - 1; k >= 0; k--) { sum -= L(i,k) * X(k,c); }
            X(i,c) = sum / L(i,i);
        }

        for (int i = N - 1; i >= 0; i--) {
            double sum = X(i,c);
            for (int k = i + 1; k < N; k++) { sum -= L(k,i) * X(k,c); }
            X(i,c) = sum / L(i,i);
        }
    }
}

//-----------------------------------------------------------------------------
// Function solves a collection of linear systems A*X = B with transposed
// right-hand size B^t, where A is the matrix whose Cholesky decomposition was
// computed by the Init() function, X and B are the same size matrices.
//-----------------------------------------------------------------------------
void BatchSolveTr(Matrix & X, const Matrix & Bt) const
{
    const Matrix & L = m_L;         // short-hand alias
    const int      N = L.NRows();   // problem size, L is square symmetric
    const int      K = X.NCols();   // number of linear systems to solve

    assert_true((N == X.NRows()) && X.SameSizeTr(Bt));

    for (int c = 0; c < K; c++) {
        for (int i = 0; i < N; i++) {
            double sum = Bt(c,i);       // transposed B
            for (int k = i - 1; k >= 0; k--) { sum -= L(i,k) * X(k,c); }
            X(i,c) = sum / L(i,i);
        }

        for (int i = N - 1; i >= 0; i--) {
            double sum = X(i,c);
            for (int k = i + 1; k < N; k++) { sum -= L(k,i) * X(k,c); }
            X(i,c) = sum / L(i,i);
        }
    }
}

}; // class Cholesky

} // namespace amdados

