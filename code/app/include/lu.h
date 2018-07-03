//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//=============================================================================
// Class for computing LU decomposition of a non-singular matrix. Once
// decomposition is done in Init() function, the class instance can be used
// to solve linear systems A*x = b and A*X = B, where x, b are vectors and X,
// B are matrices respectively.
//=============================================================================
class LUdecomposition
{
private:
    Matrix           m_LU;  // lower/upper triangular matrices of decomposition
    std::vector<int> m_Perm;// row index permutation for partial pivoting

public:
//-----------------------------------------------------------------------------
// Constructor.
//-----------------------------------------------------------------------------
LUdecomposition() : m_LU(), m_Perm()
{
}

friend std::ostream& operator<<(std::ostream& out, const LUdecomposition& lud) {
	out << "LUdecomposition: [ ";
	out << lud.m_LU;
	for(const auto& e : lud.m_Perm) { out << ", " << e; }
	out << " ]" << std::endl;
	return out;
}

//-----------------------------------------------------------------------------
// Function computes and stores LU decomposition: M = L*U.
// The decomposition P*M = L*U is based on so called partial pivoting of matrix
// rows (permutation matrix P), which is often sufficient in practice.
//-----------------------------------------------------------------------------
void Init(const Matrix & M)
{
    const double TINY = std::numeric_limits<double>::min() /
		       std::pow(std::numeric_limits<double>::epsilon(),3);

    assert_true(M.IsSquare());
    const int N = M.NRows();     // problem size; M is square

	m_LU = M;                    // copy input matrix, then do decomposition
	m_Perm.resize(N);            // resize the permutation vector accordingly

	auto & A = m_LU;             // short-hand alias for LU
	int  * P = m_Perm.data();    // short-hand alias for permutation

	// There is no permutation at the beginning.
    for (int i = 0; i < N; i++) { P[i] = i; }

	// Process column by column. The last column is trivial, so we skip it.
    for (int i = 0; i < N - 1; ++i) {
        double maxA = 0.0;
        int    imax = i;

        // Find the largest by module element from
        // the main diagonal and all way down.
        for (int k = i; k < N; ++k) {
			double absA = std::fabs(A(P[k],i));
            if (maxA < absA) {
                maxA = absA;
                imax = k;
            }
        }

        // Check the matrix is not singular.
        assert_gt(maxA, TINY) << "LU failed, max. element: " << maxA;

        // Pivoting P: the diagonal element A[P[i]][i] takes
        // the maximum value in i-th column.
        if (i != imax) std::swap(P[i], P[imax]);

		// Eliminate the elements underneath the diagonal element A[P[i]][i].
		// Everything to the left is already zero, everything to the right
		// will be combined with P[i]-th row according to the classic
        // Gauss method.
        for (int j = i + 1; j < N; ++j) {
			const int    Pj = P[j];
			const int    Pi = P[i];
            const double Aji = (A(Pj,i) /= A(Pi,i));

            for (int k = i + 1; k < N; ++k) { A(Pj,k) -= Aji * A(Pi,k); }
        }
    }
}

//-----------------------------------------------------------------------------
// Function solves a linear system A*x = b through the back-substitution, where
// A is the matrix whose LU decomposition was computed by Init() function.
//-----------------------------------------------------------------------------
void Solve(VectorView & x, const VectorView & b) const
{
	const auto & A = m_LU;              // short-hand alias
	const int  * P = m_Perm.data();     // permutation
	const int    N = A.NRows();         // problem size; A is square

	bool ok = ((x.Size() == N) && (b.Size() == N));
    assert_true(ok);

    for (int i = 0; i < N; ++i) {
        const int Pi = P[i];
        x(i) = b(Pi);
        for (int k = 0; k < i; ++k) { x(i) -= A(Pi,k) * x(k); }
    }

    for (int i = N - 1; i >= 0; --i) {
        const int Pi = P[i];
        for (int k = i + 1; k < N; ++k) { x(i) -= A(Pi,k) * x(k); }
        x(i) /= A(Pi,i);
    }
}

//-----------------------------------------------------------------------------
// Function solves a collection of linear systems A*X = B, where A is the
// matrix whose LU decomposition was computed by the Init() function,
// X and B are the matrices of the same size.
//-----------------------------------------------------------------------------
void BatchSolve(Matrix & X, const Matrix & B) const
{
    const auto & A = m_LU;              // short-hand alias
    const int  * P = m_Perm.data();     // permutation
    const int    N = A.NRows();         // problem size; A is square
    const int    K = X.NCols();         // number of linear systems to solve

    assert_true((N == X.NRows()) && X.SameSize(B));

    for (int c = 0; c < K; ++c) {
        for (int i = 0; i < N; ++i) {
            const int Pi = P[i];
            X(i,c) = B(Pi,c);
            for (int k = 0; k < i; ++k) { X(i,c) -= A(Pi,k) * X(k,c); }
        }

        for (int i = N - 1; i >= 0; --i) {
            const int Pi = P[i];
            for (int k = i + 1; k < N; ++k) { X(i,c) -= A(Pi,k) * X(k,c); }
            X(i,c) /= A(Pi,i);
        }
    }
}

//-----------------------------------------------------------------------------
// Function solves a collection of linear systems A*X = B with transposed
// right-hand size B^t, where A is the matrix whose LU decomposition was
// computed by the Init() function, X and B are the matrices of the same size.
//-----------------------------------------------------------------------------
void BatchSolveTr(Matrix & X, const Matrix & Bt) const
{
    const auto & A = m_LU;              // short-hand alias
    const int  * P = m_Perm.data();     // permutation
    const int    N = A.NRows();         // problem size; A is square
    const int    K = X.NCols();         // number of linear systems to solve

    assert_true((N == X.NRows()) && X.SameSizeTr(Bt));

    for (int c = 0; c < K; ++c) {
        for (int i = 0; i < N; ++i) {
            const int Pi = P[i];
            X(i,c) = Bt(c,Pi);      // transposed B
            for (int k = 0; k < i; ++k) { X(i,c) -= A(Pi,k) * X(k,c); }
        }

        for (int i = N - 1; i >= 0; --i) {
            const int Pi = P[i];
            for (int k = i + 1; k < N; ++k) { X(i,c) -= A(Pi,k) * X(k,c); }
            X(i,c) /= A(Pi,i);
        }
    }
}

}; // class LUdecomposition

} // namespace amdados


