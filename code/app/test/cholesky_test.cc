//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "amdados/app/amdados_grid.h"
#include "allscale/utils/assert.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/cholesky.h"
#include "amdados/app/utils/configuration.h"

using namespace amdados::app::utils;

static Configuration gConf;         // a single global configuration

// Tolerance on relative error.
const double TOL = std::sqrt(std::sqrt(std::numeric_limits<double>::epsilon()));

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers given problem size.
//-------------------------------------------------------------------------------------------------
template<int PROBLEM_SIZE>
void TestCholeskyGivenProblemSize(double & max_rel_err)
{
    const int NCOLS = PROBLEM_SIZE % 37 + 7;
    using matrix_t = allscale::utils::grid<double, PROBLEM_SIZE, PROBLEM_SIZE>;
    using matrix_other_t = allscale::utils::grid<double, PROBLEM_SIZE, NCOLS>;
    using vector_t = allscale::utils::grid<double, PROBLEM_SIZE>;

    // Create a random, square, symmetric, positive-definite matrix A: A = tmpA * tmpA^T.
    matrix_t A;
    {
        matrix_t tmpA;
        MakeRandomMatrix(tmpA);
        MatMultTransposed(A, tmpA, tmpA);
        Symmetrize(A);
        double dval = Trace(A) * TOL / PROBLEM_SIZE;
        for (int i = 0; i < PROBLEM_SIZE; ++i) { A[{i,i}] += dval; }
    }

    // Solve linear systems.
    Cholesky<PROBLEM_SIZE> chol;
    chol.ComputeDecomposition(A);

    // Compute |A*x - b| and print the relative error.
    {
        vector_t b, x, Ax;
        MakeRandomVector(b);
        chol.Solve(x, b);
        MatVecMult(Ax, A, x);
        double diff_norm = NormVecDiff(Ax, b);
        double norm_b = NormVec(b);
        max_rel_err = std::max(max_rel_err, diff_norm);
        EXPECT_LT(diff_norm / norm_b, TOL) << "Relative error |A*x - b|/|b| exceeded tolerance";
    }

    // Compute |A*X - B| and print the relative error.
    {
        matrix_other_t B, X, AX;
        MakeRandomMatrix(B);
        chol.BatchSolve(X, B);
        MatMult(AX, A, X);
        double diff_norm = NormMatDiff(AX, B);
        double norm_B = NormMat(B);
        max_rel_err = std::max(max_rel_err, diff_norm);
        EXPECT_LT(diff_norm / norm_B, TOL) << "Relative error |A*X - B|/|B| exceeded tolerance";
    }
}

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers
// for different matrix sizes.
//-------------------------------------------------------------------------------------------------
TEST(Cholesky, Basic)
{
    // Read configuration settings.
    gConf.ReadConfigFile("../../amdados_unittest.conf");
    gConf.PrintParameters();
    MakeDirectory(gConf.asCString("output_dir"));

    std::string filename = gConf.asString("output_dir") + "/cholesky_test.log";
    std::fstream fout(filename, std::ios::out | std::ios::trunc);
    assert_true(fout.good()) << "failed to oped the summary file: " << filename << std::endl;

    double max_rel_err = 0.0;
    for (int testNo = 0; testNo < 10; ++testNo) {
        TestCholeskyGivenProblemSize<11>(max_rel_err);
        TestCholeskyGivenProblemSize<17>(max_rel_err);
        TestCholeskyGivenProblemSize<23>(max_rel_err);
        TestCholeskyGivenProblemSize<37>(max_rel_err);
        TestCholeskyGivenProblemSize<41>(max_rel_err);
        TestCholeskyGivenProblemSize<53>(max_rel_err);
        TestCholeskyGivenProblemSize<67>(max_rel_err);
        TestCholeskyGivenProblemSize<73>(max_rel_err);
        TestCholeskyGivenProblemSize<89>(max_rel_err);
        TestCholeskyGivenProblemSize<97>(max_rel_err);
    }
    fout << "TestCholesky(): max. relative error: " << max_rel_err << std::endl;
    fout << std::endl << std::flush;
}

