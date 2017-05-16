#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

//=================================================================================================
// Class for computing Cholesky decomposition of a symmetric, squared, positive-definite matrix.
// Once decomposition is done in constructor, the class instance can be used to solve linear
// systems A*x = b and A*X = B, where x, b are vectors and X, B are matrices respectively.
// The implementation closely follows "Numerical Recipes" book, 3rd edition, chapter 2.9.
// The matrices and vector are assumed to be the instances of AllScale API grid.
//=================================================================================================
#include "amdados/app/utils/matrix.h"

namespace amdados {
namespace app {
namespace utils {

	template<size_t MSIZE>
	class Cholesky
	{
	public:
		using vector_t = allscale::utils::grid<double, MSIZE>;
		using matrix_t = allscale::utils::grid<double, MSIZE, MSIZE>;

	private:
		matrix_t m_L;       // lower triangular matrix of decomposition

	public:
	//-------------------------------------------------------------------------------------------------
	// Constructor.
	//-------------------------------------------------------------------------------------------------
	Cholesky()
	{
		FillMatrix(m_L, 0.0);
	}

	//-------------------------------------------------------------------------------------------------
	// Constructor computes and stores Cholesky decomposition
	// of a positive-definite symmetric matrix: A = L * L^t.
	//-------------------------------------------------------------------------------------------------
	void ComputeDecomposition(const matrix_t & A)
	{
		const int N = static_cast<int>(MSIZE);
		const double EPS = std::numeric_limits<double>::epsilon();
		const double TINY = std::numeric_limits<double>::min() / (EPS * EPS);

		// Compute the lower triangular matrix of Cholesky decomposition.
		m_L = A;
		for (int i = 0; i < N; i++) {
			for (int j = i; j < N; j++) {
				if (m_L[{j,i}] != m_L[{i,j}]) {
					throw std::runtime_error("Cholesky expects a symmetric matrix");
				}

				double sum = m_L[{i,j}];
				for (int k = i - 1; k >= 0; k--) {
					sum -= m_L[{i,k}] * m_L[{j,k}];
				}

				if (i == j) {
					if (sum <= TINY) {
						std::cout << "sum: " << sum << std::endl;
						throw std::runtime_error("Cholesky failed");
					}
					m_L[{i,i}] = std::sqrt(sum);
				} else {
					m_L[{j,i}] = sum / m_L[{i,i}];
				}
			}
		}

		// Put the upper triangular matrix to zero.
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < i; j++) {
				m_L[{j,i}] = 0.0;
			}
		}
	}

	//-------------------------------------------------------------------------------------------------
	// Function solves a linear system A*x = b, where A is the matrix whose Cholesky
	// decomposition was computed in the constructor.
	//-------------------------------------------------------------------------------------------------
	void Solve(vector_t & x, const vector_t & b)
	{
		const int N = static_cast<int>(MSIZE);

		for (int i = 0; i < N; i++) {
			double sum = b[{i}];
			for (int k = i - 1; k >= 0; k--) {
				sum -= m_L[{i,k}] * x[{k}];
			}
			x[{i}] = sum / m_L[{i,i}];
		}

		for (int i = N - 1; i >= 0; i--) {
			double sum = x[{i}];
			for (int k = i + 1; k < N; k++) {
				sum -= m_L[{k,i}] * x[{k}];
			}
			x[{i}] = sum / m_L[{i,i}];
		}
	}

	//-------------------------------------------------------------------------------------------------
	// Function solves a collection of linear systems A*X = B, where A is the matrix whose Cholesky
	// decomposition was computed in the constructor, X and B are the matrices of the same size.
	//-------------------------------------------------------------------------------------------------
	template<size_t NCOLS>
	void BatchSolve(allscale::utils::grid<double, MSIZE, NCOLS> & X,
			  const allscale::utils::grid<double, MSIZE, NCOLS> & B)
	{
		const int N = static_cast<int>(MSIZE);

		for (int c = 0; c < static_cast<int>(NCOLS); c++) {
			for (int i = 0; i < N; i++) {
				double sum = B[{i,c}];
				for (int k = i - 1; k >= 0; k--) {
					sum -= m_L[{i,k}] * X[{k,c}];
				}
				X[{i,c}] = sum / m_L[{i,i}];
			}

			for (int i = N - 1; i >= 0; i--) {
				double sum = X[{i,c}];
				for (int k = i + 1; k < N; k++) {
					sum -= m_L[{k,i}] * X[{k,c}];
				}
				X[{i,c}] = sum / m_L[{i,i}];
			}
		}
	}

	}; // class Cholesky

} // end namespace utils
} // end namespace app
} // end namespace amdados
