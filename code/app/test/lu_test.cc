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
#include "../include/amdados_utils.h"
#include "../include/matrix.h"
#include "../include/lu.h"

// Tolerance on relative error.
const double TOL = std::sqrt(std::sqrt(std::numeric_limits<double>::epsilon()));

//-----------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers
// given problem size.
//-----------------------------------------------------------------------------
void TestLUGivenProblemSize(double & max_rel_err, const int N)
{
    using namespace ::amdados;
    const int K = N % 37 + 7;

    // Create a random non-singular matrix A.
    Matrix A(N, N);
    {
        MakeRandom(A, 'u');
        double dval = Trace(A) * TOL / N;
        for (int i = 0; i < N; ++i) {
            A(i, i) += dval;
        }
    }

    // Solve linear systems.
    LUdecomposition lu;
    lu.Init(A);

    // Compute |A*x - b| and print the relative error.
    {
        Vector b(N), x(N), Ax(N);
        MakeRandom(b, 'u');
        lu.Solve(x, b);
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
        lu.BatchSolve(X, B);
        MatMult(AX, A, X);
        double diff_norm = NormDiff(AX, B);
        double norm_B = Norm(B);
        max_rel_err = std::max(max_rel_err, diff_norm);
        EXPECT_LT(diff_norm / norm_B, TOL)
                << "Relative error |A*X - B|/|B| exceeded tolerance";
    }

    // Compute |A*X - Bt| and print the relative error.
    {
        Matrix Bt(K, N), B(N, K), X(N, K), AX(N, K);
        MakeRandom(B, 'u');
        GetTransposed(Bt, B);
        lu.BatchSolveTr(X, Bt);
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
TEST(LU, Basic)
{
    // Open the output log-file.
    std::string fname = "lu_test.log";
    std::fstream log_file(fname, std::ios::out | std::ios::trunc);
    assert_true(log_file.good()) << "failed to open: " << fname << std::endl;

    double max_rel_err = 0.0;
    //TestLUGivenProblemSize<3600>(max_rel_err);
    for (int n = 1; n < 127; ++n) {
        TestLUGivenProblemSize(max_rel_err, n);
    }
    log_file << "TestLU(): max. relative error: "
             << max_rel_err << std::endl << std::endl << std::flush;
}
