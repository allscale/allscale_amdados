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
#include "amdados/app/matrix.h"

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

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers given problem size.
//-------------------------------------------------------------------------------------------------
template<int PROBLEM_SIZE>
void TestCholeskyGivenProblemSize(int numTests)
{
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> distrib(0.0,1.0);

    const int NCOLS      = PROBLEM_SIZE % 37 + 7;
    using matrix_t       = allscale::utils::grid<double, PROBLEM_SIZE, PROBLEM_SIZE>;
    using matrix_other_t = allscale::utils::grid<double, PROBLEM_SIZE, NCOLS>;
    using vector_t       = allscale::utils::grid<double, PROBLEM_SIZE>;

    matrix_t       A, tmpA;
    matrix_other_t B, X;
    vector_t       b, x;

    std::cout << "Testing Cholesky decomposition, problem size: " << PROBLEM_SIZE
              << std::endl << std::endl;
    for (int testNo = 0; testNo < numTests; ++testNo) {
        double norm_b = 0.0, norm_B = 0.0, diff_norm = 0.0;

        // Create a random, square, symmetric, positive-definite matrix A: A = tmpA * tmpA.
        for (int i = 0; i < PROBLEM_SIZE; ++i) {
            for (int j = 0; j < PROBLEM_SIZE; ++j) {
                tmpA[{i,j}] = distrib(gen);
            }
        }
        Symmetrize(tmpA);
        for (int r = 0; r < PROBLEM_SIZE; ++r) {
            for (int c = 0; c < PROBLEM_SIZE; ++c) {
                double v = 0.0;
                for (int i = 0; i < PROBLEM_SIZE; ++i) {
                    v += tmpA[{r,i}] * tmpA[{i,c}];
                }
                A[{r,c}] = v;
            }
        }
        Symmetrize(A);

        // Create the right-hand side matrix and vector.
        for (int r = 0; r < PROBLEM_SIZE; ++r) {
            for (int c = 0; c < NCOLS; ++c) {
                double v = distrib(gen);
                B[{r,c}] = v;
                norm_B += v * v;
            }
            double v = distrib(gen);
            b[{r}] = v;
            norm_b += v * v;
        }
        norm_b = std::sqrt(norm_b);
        norm_B = std::sqrt(norm_B);

        // Solve linear systems.
        Cholesky<PROBLEM_SIZE> chol;
        chol.ComputeDecomposition(A);
        chol.Solve(x, b);
        chol.BatchSolve(X, B);

        // Compute |A*x - b| and print the relative error.
        diff_norm = 0.0;
        for (int r = 0; r < PROBLEM_SIZE; ++r) {
            double v = 0.0;
            for (int i = 0; i < PROBLEM_SIZE; ++i) {
                v += A[{r,i}] * x[{i}];
            }
            v -= b[{r}];
            diff_norm += v * v;
        }
        diff_norm = std::sqrt(diff_norm);
        std::cout << "Relative error |A*x - b|/|b|: " << (diff_norm/norm_b) << std::endl;

        // Compute |A*X - B| and print the relative error.
        diff_norm = 0.0;
        for (int r = 0; r < PROBLEM_SIZE; ++r) {
            for (int c = 0; c < NCOLS; ++c) {
                double v = 0.0;
                for (int i = 0; i < PROBLEM_SIZE; ++i) {
                    v += A[{r,i}] * X[{i,c}];
                }
                v -= B[{r,c}];
                diff_norm += v * v;
            }
        }
        diff_norm = std::sqrt(diff_norm);
        std::cout << "Relative error |A*X - B|/|B|: " << (diff_norm/norm_B) << std::endl;
    }
    std::cout << std::endl;
}

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers
// for different matrix sizes.
//-------------------------------------------------------------------------------------------------
void TestCholesky()
{
    const int numTests = 10;
    TestCholeskyGivenProblemSize<11>(numTests);
    TestCholeskyGivenProblemSize<17>(numTests);
    TestCholeskyGivenProblemSize<23>(numTests);
    TestCholeskyGivenProblemSize<37>(numTests);
    TestCholeskyGivenProblemSize<41>(numTests);
    TestCholeskyGivenProblemSize<53>(numTests);
    TestCholeskyGivenProblemSize<67>(numTests);
    TestCholeskyGivenProblemSize<73>(numTests);
    TestCholeskyGivenProblemSize<89>(numTests);
    TestCholeskyGivenProblemSize<97>(numTests);
}

