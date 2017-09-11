//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "allscale/api/user/data/adaptive_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/geometry.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/cholesky.h"
#include "amdados/app/utils/configuration.h"

using namespace amdados::app;
using namespace amdados::app::utils;

// Tolerance on relative error.
const double TOL = std::pow(numeric_limits<double>::epsilon(), 0.25);

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers given problem size.
//-------------------------------------------------------------------------------------------------
template<int PROBLEM_SIZE>
void TestCholeskyGivenProblemSize(double & max_rel_err)
{
    const int NCOLS = PROBLEM_SIZE % 37 + 7;
    using matrix_t = Matrix<PROBLEM_SIZE, PROBLEM_SIZE>;
    using matrix_other_t = Matrix<PROBLEM_SIZE, NCOLS>;
    using vector_t = Vector<PROBLEM_SIZE>;

    // Create a random, square, symmetric, positive-definite matrix A: A = tmpA * tmpA^T.
    matrix_t A;
    {
        matrix_t tmpA;
        MakeRandomMatrix(tmpA);
        MatMultTr(A, tmpA, tmpA);
        Symmetrize(A);
        double dval = Trace(A) * TOL / PROBLEM_SIZE;
        for (int i = 0; i < PROBLEM_SIZE; ++i) { A(i,i) += dval; }
    }

    // Solve linear systems.
    Cholesky<PROBLEM_SIZE> chol;
    chol.Init(A);

    // Compute |A*x - b| and print the relative error.
    {
        vector_t b, x, Ax;
        MakeRandomVector(b, 'u');
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
    Configuration conf;
    conf.ReadConfigFile("../../amdados.conf");
    conf.PrintParameters();
    MakeDirectory(conf.asCString("test_output_dir"));

    std::string filename = conf.asString("test_output_dir") + "/cholesky_test.log";
    std::fstream fout(filename, std::ios::out | std::ios::trunc);
    assert_true(fout.good()) << "failed to oped the summary file: " << filename << std::endl;

    double max_rel_err = 0.0;
    for (int testNo = 0; testNo < 3; ++testNo) {
        TestCholeskyGivenProblemSize< 1>(max_rel_err);
        TestCholeskyGivenProblemSize< 2>(max_rel_err);
        TestCholeskyGivenProblemSize< 3>(max_rel_err);
        TestCholeskyGivenProblemSize< 4>(max_rel_err);
        TestCholeskyGivenProblemSize< 5>(max_rel_err);
        TestCholeskyGivenProblemSize< 6>(max_rel_err);
        TestCholeskyGivenProblemSize< 7>(max_rel_err);
        TestCholeskyGivenProblemSize< 8>(max_rel_err);
        TestCholeskyGivenProblemSize< 9>(max_rel_err);
        TestCholeskyGivenProblemSize<10>(max_rel_err);
        TestCholeskyGivenProblemSize<11>(max_rel_err);
        TestCholeskyGivenProblemSize<12>(max_rel_err);
        TestCholeskyGivenProblemSize<13>(max_rel_err);
        TestCholeskyGivenProblemSize<14>(max_rel_err);
        TestCholeskyGivenProblemSize<15>(max_rel_err);
        TestCholeskyGivenProblemSize<16>(max_rel_err);
        TestCholeskyGivenProblemSize<17>(max_rel_err);
        TestCholeskyGivenProblemSize<18>(max_rel_err);
        TestCholeskyGivenProblemSize<19>(max_rel_err);
        TestCholeskyGivenProblemSize<20>(max_rel_err);
        TestCholeskyGivenProblemSize<21>(max_rel_err);
        TestCholeskyGivenProblemSize<22>(max_rel_err);
        TestCholeskyGivenProblemSize<23>(max_rel_err);
        TestCholeskyGivenProblemSize<24>(max_rel_err);
        TestCholeskyGivenProblemSize<25>(max_rel_err);
        TestCholeskyGivenProblemSize<26>(max_rel_err);
        TestCholeskyGivenProblemSize<27>(max_rel_err);
        TestCholeskyGivenProblemSize<28>(max_rel_err);
        TestCholeskyGivenProblemSize<29>(max_rel_err);
        TestCholeskyGivenProblemSize<30>(max_rel_err);
        TestCholeskyGivenProblemSize<31>(max_rel_err);
        TestCholeskyGivenProblemSize<32>(max_rel_err);
        TestCholeskyGivenProblemSize<33>(max_rel_err);
        TestCholeskyGivenProblemSize<34>(max_rel_err);
        TestCholeskyGivenProblemSize<35>(max_rel_err);
        TestCholeskyGivenProblemSize<36>(max_rel_err);
        TestCholeskyGivenProblemSize<37>(max_rel_err);
        TestCholeskyGivenProblemSize<38>(max_rel_err);
        TestCholeskyGivenProblemSize<39>(max_rel_err);
        TestCholeskyGivenProblemSize<40>(max_rel_err);
        TestCholeskyGivenProblemSize<41>(max_rel_err);
        TestCholeskyGivenProblemSize<42>(max_rel_err);
        TestCholeskyGivenProblemSize<43>(max_rel_err);
        TestCholeskyGivenProblemSize<44>(max_rel_err);
        TestCholeskyGivenProblemSize<45>(max_rel_err);
        TestCholeskyGivenProblemSize<46>(max_rel_err);
        TestCholeskyGivenProblemSize<47>(max_rel_err);
        TestCholeskyGivenProblemSize<48>(max_rel_err);
        TestCholeskyGivenProblemSize<49>(max_rel_err);
        TestCholeskyGivenProblemSize<50>(max_rel_err);
        TestCholeskyGivenProblemSize<51>(max_rel_err);
        TestCholeskyGivenProblemSize<52>(max_rel_err);
        TestCholeskyGivenProblemSize<53>(max_rel_err);
        TestCholeskyGivenProblemSize<54>(max_rel_err);
        TestCholeskyGivenProblemSize<55>(max_rel_err);
        TestCholeskyGivenProblemSize<56>(max_rel_err);
        TestCholeskyGivenProblemSize<57>(max_rel_err);
        TestCholeskyGivenProblemSize<58>(max_rel_err);
        TestCholeskyGivenProblemSize<59>(max_rel_err);
        TestCholeskyGivenProblemSize<60>(max_rel_err);
        TestCholeskyGivenProblemSize<61>(max_rel_err);
        TestCholeskyGivenProblemSize<62>(max_rel_err);
        TestCholeskyGivenProblemSize<63>(max_rel_err);
        TestCholeskyGivenProblemSize<64>(max_rel_err);
        TestCholeskyGivenProblemSize<65>(max_rel_err);
        TestCholeskyGivenProblemSize<66>(max_rel_err);
        TestCholeskyGivenProblemSize<67>(max_rel_err);
        TestCholeskyGivenProblemSize<68>(max_rel_err);
        TestCholeskyGivenProblemSize<69>(max_rel_err);
        TestCholeskyGivenProblemSize<70>(max_rel_err);
        TestCholeskyGivenProblemSize<71>(max_rel_err);
        TestCholeskyGivenProblemSize<72>(max_rel_err);
        TestCholeskyGivenProblemSize<73>(max_rel_err);
        TestCholeskyGivenProblemSize<74>(max_rel_err);
        TestCholeskyGivenProblemSize<75>(max_rel_err);
        TestCholeskyGivenProblemSize<76>(max_rel_err);
        TestCholeskyGivenProblemSize<77>(max_rel_err);
        TestCholeskyGivenProblemSize<78>(max_rel_err);
        TestCholeskyGivenProblemSize<79>(max_rel_err);
        TestCholeskyGivenProblemSize<80>(max_rel_err);
        TestCholeskyGivenProblemSize<81>(max_rel_err);
        TestCholeskyGivenProblemSize<82>(max_rel_err);
        TestCholeskyGivenProblemSize<83>(max_rel_err);
        TestCholeskyGivenProblemSize<84>(max_rel_err);
        TestCholeskyGivenProblemSize<85>(max_rel_err);
        TestCholeskyGivenProblemSize<86>(max_rel_err);
        TestCholeskyGivenProblemSize<87>(max_rel_err);
        TestCholeskyGivenProblemSize<88>(max_rel_err);
        TestCholeskyGivenProblemSize<89>(max_rel_err);
        TestCholeskyGivenProblemSize<90>(max_rel_err);
        TestCholeskyGivenProblemSize<91>(max_rel_err);
        TestCholeskyGivenProblemSize<92>(max_rel_err);
        TestCholeskyGivenProblemSize<93>(max_rel_err);
        TestCholeskyGivenProblemSize<94>(max_rel_err);
        TestCholeskyGivenProblemSize<95>(max_rel_err);
        TestCholeskyGivenProblemSize<96>(max_rel_err);
        TestCholeskyGivenProblemSize<97>(max_rel_err);
        TestCholeskyGivenProblemSize<98>(max_rel_err);
        TestCholeskyGivenProblemSize<99>(max_rel_err);
    }
    fout << "TestCholesky(): max. relative error: " << max_rel_err << std::endl;
    fout << std::endl << std::flush;
}

