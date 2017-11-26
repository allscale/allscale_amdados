//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "allscale/utils/assert.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/cholesky.h"

using namespace ::amdados::app;
using namespace ::amdados::app::utils;

// Tolerance on relative error.
const double TOL = std::pow(std::numeric_limits<double>::epsilon(), 0.25);

//-----------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers
// given problem size.
//-----------------------------------------------------------------------------
void TestCholeskyGivenProblemSize(double & max_rel_err, const int N)
{
    const int K = N % 37 + 7;

    // Create a random, square, symmetric, positive-definite matrix A:
    // A = tmpA * tmpA^T.
    Matrix A(N, N);
    {
        Matrix tmpA(N, N);
        MakeRandom(tmpA, 'u');
        MatMultTr(A, tmpA, tmpA);
        Symmetrize(A);
        double dval = Trace(A) * TOL / N;
        for (int i = 0; i < N; ++i) {
            A(i, i) += dval;
        }
    }

    // Solve linear systems.
    Cholesky chol;
    chol.Init(A);

    // Compute |A*x - b| and print the relative error.
    {
        Vector b(N), x(N), Ax(N);
        MakeRandom(b, 'u');
        chol.Solve(x, b);
        MatVecMult(Ax, A, x);
        double diff_norm = NormDiff(Ax, b);
        double norm_b = Norm(b);
        max_rel_err = std::max(max_rel_err, diff_norm);
        EXPECT_LT(diff_norm / norm_b, TOL)
                << "Relative error |A*x - b|/|b| exceeded tolerance";
    }

    // Compute |A*X - B| and print the relative error.
    {
        Matrix B(N, K), X(N, K), AX(N, K);
        MakeRandom(B, 'u');
        chol.BatchSolve(X, B);
        MatMult(AX, A, X);
        double diff_norm = NormDiff(AX, B);
        double norm_B = Norm(B);
        max_rel_err = std::max(max_rel_err, diff_norm);
        EXPECT_LT(diff_norm / norm_B, TOL)
                << "Relative error |A*X - B|/|B| exceeded tolerance";
    }

    // Compute |A*X - B^t| and print the relative error.
    {
        Matrix Bt(K, N), B(N, K), X(N, K), AX(N, K);
        MakeRandom(B, 'u');
        GetTransposed(Bt, B);
        chol.BatchSolveTr(X, Bt);
        MatMult(AX, A, X);
        double diff_norm = NormDiff(AX, B);
        double norm_B = Norm(Bt);
        max_rel_err = std::max(max_rel_err, diff_norm);
        EXPECT_LT(diff_norm / norm_B, TOL)
                << "Relative error |A*X - Bt|/|Bt| exceeded tolerance";
    }
}

//-----------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers
// for different matrix sizes.
//-----------------------------------------------------------------------------
TEST(Cholesky, Basic)
{
    // Open the output log-file.
    std::fstream log_file;
    OpenTextFileForUnitTest(log_file, "cholesky_test.log");

    double max_rel_err = 0.0;
    for (int n = 1; n < 127; ++n) {
        TestCholeskyGivenProblemSize(max_rel_err, n);
    }
    log_file << "TestCholesky(): max. relative error: "
             << max_rel_err << std::endl << std::endl << std::flush;
}

