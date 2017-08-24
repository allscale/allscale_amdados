//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {
namespace utils {

//=================================================================================================
// Class for computing LU decomposition of a non-singular matrix. Once decomposition is done
// in Init() function, the class instance can be used to solve linear systems A*x = b and
// A*X = B, where x, b are vectors and X, B are matrices respectively.
//=================================================================================================
template<int MSIZE>
class LUdecomposition
{
private:
    Matrix<MSIZE,MSIZE> m_LU;       // lower and upper triangular matrices of decomposition
    std::vector<int>    m_Perm;     // row index permutation for partial pivoting

public:
//-------------------------------------------------------------------------------------------------
// Constructor.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
LUdecomposition() : m_LU()
{
    m_Perm.resize(MSIZE, 0);
}

//-------------------------------------------------------------------------------------------------
// Function computes and stores LU decomposition: M = L*U.
// The decomposition P*M = L*U is based on so called partial pivoting of matrix row (permutation
// matrix P), which is often sufficient in practice.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
void Init(const Matrix<MSIZE,MSIZE> & M)
{
    const double TINY = numeric_limits<double>::min() /
		       std::pow(numeric_limits<double>::epsilon(),3);

	m_LU = M;                   // copy input matrix, then carry out in-place decomposition

	auto & A = m_LU;            // short-hand aliases
	int  * P = m_Perm.data();

	// There is no permutation at the beginning.
    for (int i = 0; i < MSIZE; i++) { P[i] = i; }

	// Process column by column. The last column is trivial, so we skip it.
    for (int i = 0; i < MSIZE - 1; ++i) {
        double maxA = 0.0;
        int    imax = i;

        // Find the largest by module element from the main diagonal and all way down.
        for (int k = i; k < MSIZE; ++k) {
			double absA = fabs(A(P[k],i));
            if (maxA < absA) {
                maxA = absA;
                imax = k;
            }
        }

        // Check the matrix is not singular.
        assert_gt(maxA, TINY) << "LU failed, max. element: " << maxA;

        // Pivoting P: the diagonal element A[P[i]][i] takes the maximum value in i-th column.
        if (i != imax) std::swap(P[i], P[imax]);

		// Eliminate the elements undernearth the diagonal element A[P[i]][i].
		// Everything to the left is already zero, everything to the right
		// will be combined with P[i]-th row according to the classic Gauss method.
        for (int j = i + 1; j < MSIZE; ++j) {
			const int    Pj = P[j];
			const int    Pi = P[i];
            const double Aji = (A(Pj,i) /= A(Pi,i));

            for (int k = i + 1; k < MSIZE; ++k) { A(Pj,k) -= Aji * A(Pi,k); }
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function solves a linear system A*x = b through the backsubstitution, where A is
// the matrix whose LU decomposition was computed by Init() function.
//-------------------------------------------------------------------------------------------------
IBM_NOINLINE
void Solve(VectorView<MSIZE> & x, const VectorView<MSIZE> & b) const
{
	const auto & A = m_LU;                  // short-hand aliases
	const int  * P = m_Perm.data();

    for (int i = 0; i < MSIZE; ++i) {
        const int Pi = P[i];
        x(i) = b(Pi);
        for (int k = 0; k < i; ++k) { x(i) -= A(Pi,k) * x(k); }
    }

    for (int i = MSIZE - 1; i >= 0; --i) {
        const int Pi = P[i];
        for (int k = i + 1; k < MSIZE; ++k) { x(i) -= A(Pi,k) * x(k); }
        x(i) /= A(Pi,i);
    }
}

//-------------------------------------------------------------------------------------------------
// Function solves a collection of linear systems A*X = B, where A is the matrix whose LU
// decomposition was computed by the Init() function, X and B are the matrices of the same size.
//-------------------------------------------------------------------------------------------------
template<int NCOLS>
IBM_NOINLINE
void BatchSolve(Matrix<MSIZE,NCOLS> & X, const Matrix<MSIZE,NCOLS> & B) const
{
	const auto & A = m_LU;                  // short-hand aliases
	const int  * P = m_Perm.data();

    for (int c = 0; c < NCOLS; ++c) {
        for (int i = 0; i < MSIZE; ++i) {
            const int Pi = P[i];
            X(i,c) = B(Pi,c);
            for (int k = 0; k < i; ++k) { X(i,c) -= A(Pi,k) * X(k,c); }
        }

        for (int i = MSIZE - 1; i >= 0; --i) {
            const int Pi = P[i];
            for (int k = i + 1; k < MSIZE; ++k) { X(i,c) -= A(Pi,k) * X(k,c); }
            X(i,c) /= A(Pi,i);
        }
    }
}

}; // class LUdecomposition

} // end namespace utils
} // end namespace app
} // end namespace amdados


