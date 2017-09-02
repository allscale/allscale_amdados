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
#include "amdados/app/utils/lu.h"
#include "amdados/app/utils/configuration.h"

using namespace amdados::app;
using namespace amdados::app::utils;

// Tolerance on relative error.
const double TOL = std::sqrt(std::sqrt(std::numeric_limits<double>::epsilon()));

//-------------------------------------------------------------------------------------------------
// Function tests the Cholesky decomposition and the linear systems solvers given problem size.
//-------------------------------------------------------------------------------------------------
template<int PROBLEM_SIZE>
void TestLUGivenProblemSize(double & max_rel_err)
{
    const int NCOLS = PROBLEM_SIZE % 37 + 7;
    using matrix_t = Matrix<PROBLEM_SIZE, PROBLEM_SIZE>;
    using matrix_other_t = Matrix<PROBLEM_SIZE, NCOLS>;
    using vector_t = Vector<PROBLEM_SIZE>;

    // Create a random non-singular matrix A.
    matrix_t A;
    {
        MakeRandomMatrix(A);
        double dval = Trace(A) * TOL / PROBLEM_SIZE;
        for (int i = 0; i < PROBLEM_SIZE; ++i) { A(i,i) += dval; }
    }

    // Solve linear systems.
    LUdecomposition<PROBLEM_SIZE> lu;
    lu.Init(A);

    // Compute |A*x - b| and print the relative error.
    {
        vector_t b, x, Ax;
        MakeRandomVector(b, 'u');
        lu.Solve(x, b);
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
        lu.BatchSolve(X, B);
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
TEST(LU, Basic)
{
    // Read configuration settings.
    Configuration conf;
    conf.ReadConfigFile("../../amdados.conf");
    conf.PrintParameters();
    MakeDirectory(conf.asCString("test_output_dir"));

    std::string filename = conf.asString("test_output_dir") + "/lu_test.log";
    std::fstream fout(filename, std::ios::out | std::ios::trunc);
    assert_true(fout.good()) << "failed to oped the summary file: " << filename << std::endl;

    double max_rel_err = 0.0;
    //TestLUGivenProblemSize<3600>(max_rel_err);
    for (int testNo = 0; testNo < 3; ++testNo) {
        TestLUGivenProblemSize< 1>(max_rel_err);
        TestLUGivenProblemSize< 2>(max_rel_err);
        TestLUGivenProblemSize< 3>(max_rel_err);
        TestLUGivenProblemSize< 4>(max_rel_err);
        TestLUGivenProblemSize< 5>(max_rel_err);
        TestLUGivenProblemSize< 6>(max_rel_err);
        TestLUGivenProblemSize< 7>(max_rel_err);
        TestLUGivenProblemSize< 8>(max_rel_err);
        TestLUGivenProblemSize< 9>(max_rel_err);
        TestLUGivenProblemSize<10>(max_rel_err);
        TestLUGivenProblemSize<11>(max_rel_err);
        TestLUGivenProblemSize<12>(max_rel_err);
        TestLUGivenProblemSize<13>(max_rel_err);
        TestLUGivenProblemSize<14>(max_rel_err);
        TestLUGivenProblemSize<15>(max_rel_err);
        TestLUGivenProblemSize<16>(max_rel_err);
        TestLUGivenProblemSize<17>(max_rel_err);
        TestLUGivenProblemSize<18>(max_rel_err);
        TestLUGivenProblemSize<19>(max_rel_err);
        TestLUGivenProblemSize<20>(max_rel_err);
        TestLUGivenProblemSize<21>(max_rel_err);
        TestLUGivenProblemSize<22>(max_rel_err);
        TestLUGivenProblemSize<23>(max_rel_err);
        TestLUGivenProblemSize<24>(max_rel_err);
        TestLUGivenProblemSize<25>(max_rel_err);
        TestLUGivenProblemSize<26>(max_rel_err);
        TestLUGivenProblemSize<27>(max_rel_err);
        TestLUGivenProblemSize<28>(max_rel_err);
        TestLUGivenProblemSize<29>(max_rel_err);
        TestLUGivenProblemSize<30>(max_rel_err);
        TestLUGivenProblemSize<31>(max_rel_err);
        TestLUGivenProblemSize<32>(max_rel_err);
        TestLUGivenProblemSize<33>(max_rel_err);
        TestLUGivenProblemSize<34>(max_rel_err);
        TestLUGivenProblemSize<35>(max_rel_err);
        TestLUGivenProblemSize<36>(max_rel_err);
        TestLUGivenProblemSize<37>(max_rel_err);
        TestLUGivenProblemSize<38>(max_rel_err);
        TestLUGivenProblemSize<39>(max_rel_err);
        TestLUGivenProblemSize<40>(max_rel_err);
        TestLUGivenProblemSize<41>(max_rel_err);
        TestLUGivenProblemSize<42>(max_rel_err);
        TestLUGivenProblemSize<43>(max_rel_err);
        TestLUGivenProblemSize<44>(max_rel_err);
        TestLUGivenProblemSize<45>(max_rel_err);
        TestLUGivenProblemSize<46>(max_rel_err);
        TestLUGivenProblemSize<47>(max_rel_err);
        TestLUGivenProblemSize<48>(max_rel_err);
        TestLUGivenProblemSize<49>(max_rel_err);
        TestLUGivenProblemSize<50>(max_rel_err);
        TestLUGivenProblemSize<51>(max_rel_err);
        TestLUGivenProblemSize<52>(max_rel_err);
        TestLUGivenProblemSize<53>(max_rel_err);
        TestLUGivenProblemSize<54>(max_rel_err);
        TestLUGivenProblemSize<55>(max_rel_err);
        TestLUGivenProblemSize<56>(max_rel_err);
        TestLUGivenProblemSize<57>(max_rel_err);
        TestLUGivenProblemSize<58>(max_rel_err);
        TestLUGivenProblemSize<59>(max_rel_err);
        TestLUGivenProblemSize<60>(max_rel_err);
        TestLUGivenProblemSize<61>(max_rel_err);
        TestLUGivenProblemSize<62>(max_rel_err);
        TestLUGivenProblemSize<63>(max_rel_err);
        TestLUGivenProblemSize<64>(max_rel_err);
        TestLUGivenProblemSize<65>(max_rel_err);
        TestLUGivenProblemSize<66>(max_rel_err);
        TestLUGivenProblemSize<67>(max_rel_err);
        TestLUGivenProblemSize<68>(max_rel_err);
        TestLUGivenProblemSize<69>(max_rel_err);
        TestLUGivenProblemSize<70>(max_rel_err);
        TestLUGivenProblemSize<71>(max_rel_err);
        TestLUGivenProblemSize<72>(max_rel_err);
        TestLUGivenProblemSize<73>(max_rel_err);
        TestLUGivenProblemSize<74>(max_rel_err);
        TestLUGivenProblemSize<75>(max_rel_err);
        TestLUGivenProblemSize<76>(max_rel_err);
        TestLUGivenProblemSize<77>(max_rel_err);
        TestLUGivenProblemSize<78>(max_rel_err);
        TestLUGivenProblemSize<79>(max_rel_err);
        TestLUGivenProblemSize<80>(max_rel_err);
        TestLUGivenProblemSize<81>(max_rel_err);
        TestLUGivenProblemSize<82>(max_rel_err);
        TestLUGivenProblemSize<83>(max_rel_err);
        TestLUGivenProblemSize<84>(max_rel_err);
        TestLUGivenProblemSize<85>(max_rel_err);
        TestLUGivenProblemSize<86>(max_rel_err);
        TestLUGivenProblemSize<87>(max_rel_err);
        TestLUGivenProblemSize<88>(max_rel_err);
        TestLUGivenProblemSize<89>(max_rel_err);
        TestLUGivenProblemSize<90>(max_rel_err);
        TestLUGivenProblemSize<91>(max_rel_err);
        TestLUGivenProblemSize<92>(max_rel_err);
        TestLUGivenProblemSize<93>(max_rel_err);
        TestLUGivenProblemSize<94>(max_rel_err);
        TestLUGivenProblemSize<95>(max_rel_err);
        TestLUGivenProblemSize<96>(max_rel_err);
        TestLUGivenProblemSize<97>(max_rel_err);
        TestLUGivenProblemSize<98>(max_rel_err);
        TestLUGivenProblemSize<99>(max_rel_err);
    }
    fout << "TestLU(): max. relative error: " << max_rel_err << std::endl;
    fout << std::endl << std::flush;
}

