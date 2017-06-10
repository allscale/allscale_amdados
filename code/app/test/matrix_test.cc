//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>

#define ARMA_USE_CXX11
#define ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_LAPACK
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_NEWARP
#define ARMA_DONT_USE_ARPACK
#define ARMA_DONT_USE_SUPERLU
#define ARMA_DONT_USE_HDF5
#define ARMA_DONT_USE_OPENMP
#include "../../../3rd-party/include/armadillo"

#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"
#include "amdados/app/utils/common.h"
#include "amdados/app/utils/amdados_utils.h"
#include "amdados/app/utils/configuration.h"
#include "amdados/app/utils/matrix.h"
#include "amdados/app/utils/sparse_matrix.h"

using namespace amdados::app::utils;

static Configuration gConf;         // a single global configuration

// Tolerance on relative error.
const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());

//-------------------------------------------------------------------------------------------------
// Function copies an "allscale" matrix to "armadillo" one.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void CopyToArma(arma::mat & dest, const Matrix<NROWS,NCOLS> & source)
{
    dest.zeros(NROWS, NCOLS);
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            dest(r,c) = source[{r,c}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function copies an "allscale" vector to "armadillo" one.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
void CopyToArma(arma::vec & dest, const Vector<LENGTH> & source)
{
    dest.zeros(LENGTH);
    for (int k = 0; k < static_cast<int>(LENGTH); ++k) {
        dest(k) = source[{k}];
    }
}

//-------------------------------------------------------------------------------------------------
// Function initializes a random "allscale" matrix along with corresponding "armadillo" one.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void InitRandomSparse(SpMatrix<NROWS,NCOLS> & A, arma::sp_mat & B)
{
    // Make "allscale" sparse matrix.
    std::vector<Triplet> triplets = MakeRandomSpMatrix(A, 0.1);

    // Make "armadillo" sparse matrix.
    arma::umat locations(2, triplets.size());
    arma::vec  values(triplets.size());
    for (arma::uword k = 0; k < triplets.size(); ++k) {
        locations(0,k) = triplets[k].row();
        locations(1,k) = triplets[k].col();
        values(k)      = triplets[k].value();
    }
    bool add_values = true;
    B = arma::sp_mat(add_values, locations, values, static_cast<arma::uword>(A.NRows()),
                                                    static_cast<arma::uword>(A.NCols()));
}

//-------------------------------------------------------------------------------------------------
// Function compares matrices for equality and updates the max. relative error.
//-------------------------------------------------------------------------------------------------
double RelativeErrorCheckAndUpdate(const arma::mat & A, const arma::mat & B,
                                   double max_rel_err, const char * what)
{
    double rel_err = arma::norm(A - B, "fro") / arma::norm(A, "fro");
    max_rel_err = std::max(max_rel_err, rel_err);
    EXPECT_LT(rel_err, TOL) << "Relative error of " << what << " exceeded tolerance";
    return max_rel_err;
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double TestDenseMulDense(double max_rel_err)
{
    const size_t MSIZE = 167;

    Matrix<NROWS,MSIZE> A1;
    Matrix<MSIZE,NCOLS> A2;
    Matrix<NROWS,NCOLS> A_mul;
    arma::mat           B1, B2, B_mul, A_tmp;

    // Create random matrices.
    MakeRandomMatrix(A1);
    MakeRandomMatrix(A2);
    CopyToArma(B1, A1);
    CopyToArma(B2, A2);

    // Check multiplication.
    MatMult(A_mul, A1, A2);         // A_mul = A1 * A2
    B_mul = B1 * B2;                // B_mul = B1 * B2
    CopyToArma(A_tmp, A_mul);       // A_tmp = A_mul

    return RelativeErrorCheckAndUpdate(B_mul, A_tmp, max_rel_err, "dense*dense multiplication");
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double TestSparseMulDense(double max_rel_err)
{
    const size_t MSIZE = 111;

    // Make "allscale" and "armadillo" sparse matrices.
    SpMatrix<NROWS,MSIZE> A;
    arma::sp_mat B;
    InitRandomSparse(A, B);

    // Auxiliary matrices.
    Matrix<MSIZE,NCOLS> A_right;
    Matrix<NROWS,NCOLS> A_mul;
    arma::mat           B_right, B_mul, A_tmp;

    // Make random, dense matrices.
    MakeRandomMatrix(A_right);              // A_right <-- random dense
    CopyToArma(B_right, A_right);           // B_right = A_right

    // Check right-hand side sparse to dense multiplication.
    SparseMulDense(A_mul, A, A_right);      // A_mul = A * A_right
    B_mul = B * B_right;                    // B_mul = B * B_right
    CopyToArma(A_tmp, A_mul);               // A_tmp = A_mul

    return RelativeErrorCheckAndUpdate(B_mul, A_tmp, max_rel_err, "sparse*dense multiplication");
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double TestDenseMulSparse(double max_rel_err)
{
    const size_t MSIZE = 143;

    // Make "allscale" and "armadillo" sparse matrices.
    SpMatrix<MSIZE,NCOLS> A;
    arma::sp_mat B;
    InitRandomSparse(A, B);

    // Auxiliary matrices.
    Matrix<NROWS,MSIZE> A_left;
    Matrix<NROWS,NCOLS> A_mul;
    arma::mat           B_left, B_mul, A_tmp;

    // Make random, dense matrices.
    MakeRandomMatrix(A_left);               // A_left <-- random dense
    CopyToArma(B_left, A_left);             // B_left = A_left

    // Check left-hand side sparse to dense multiplication.
    DenseMulSparse(A_mul, A_left, A);       // A_mul = A_left * A
    B_mul = B_left * B;                     // B_mul = B_left * B
    CopyToArma(A_tmp, A_mul);               // A_tmp = A_mul

    return RelativeErrorCheckAndUpdate(B_mul, A_tmp, max_rel_err, "dense*sparse multiplication");
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double TestDenseMulSparseTransposed(double max_rel_err)
{
    const size_t MSIZE = 119;

    // Make "allscale" and "armadillo" sparse matrices.
    SpMatrix<NCOLS,MSIZE> A;
    arma::sp_mat B;
    InitRandomSparse(A, B);

    // Auxiliary matrices.
    Matrix<NROWS,MSIZE> A_left;
    Matrix<NROWS,NCOLS> A_mul;
    arma::mat           B_left, B_mul, A_tmp;

    // Make random, dense matrices.
    MakeRandomMatrix(A_left);               // A_left <-- random dense
    CopyToArma(B_left, A_left);             // B_left = A_left

    // Check left-hand side sparse to dense multiplication.
    DenseMulSparseTr(A_mul, A_left, A);     // A_mul = A_left * A^T
    B_mul = B_left * B.t();                 // B_mul = B_left * B^T
    CopyToArma(A_tmp, A_mul);               // A_tmp = A_mul

    return RelativeErrorCheckAndUpdate(B_mul, A_tmp, max_rel_err, "dense*sparse^T multiplication");
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double TestMatrixMulVector(double max_rel_err)
{
    Vector<NROWS> result;
    Vector<NCOLS> source;
    arma::vec res1, res2, src;
    MakeRandomVector(source);
    CopyToArma(src, source);

    {
        SpMatrix<NROWS,NCOLS> A;
        arma::sp_mat B;
        InitRandomSparse(A, B);

        SparseMulVector(result, A, source);
        CopyToArma(res1, result);
        res2 = B * src;
        max_rel_err = RelativeErrorCheckAndUpdate(res1, res2, max_rel_err,
                                                    "sparse*vector multiplication");
    }

    {
        Matrix<NROWS,NCOLS> A;
        arma::mat B;
        MakeRandomMatrix(A);
        CopyToArma(B, A);

        MatVecMult(result, A, source);
        CopyToArma(res1, result);
        res2 = B * src;
        max_rel_err = RelativeErrorCheckAndUpdate(res1, res2, max_rel_err,
                                                    "dense*vector multiplication");
    }
    return max_rel_err;
}

//-------------------------------------------------------------------------------------------------
// Function tests the matrix library written by means of Allscale API.
//-------------------------------------------------------------------------------------------------
TEST(MatrixTests, Basic)
{
    // Read configuration settings.
    gConf.ReadConfigFile("../../amdados.conf");
    gConf.PrintParameters();
    MakeDirectory(gConf.asCString("test_output_dir"));

    std::string filename = gConf.asString("test_output_dir") + "/matrix_test.log";
    std::fstream fout(filename, std::ios::out | std::ios::trunc);
    assert_true(fout.good()) << "failed to oped the summary file: " << filename << std::endl;

    {
        double max_rel_err = 0.0;
        for (int testNo = 0; testNo < 10; ++testNo) {
            max_rel_err = TestSparseMulDense<37,67>(max_rel_err);
            max_rel_err = TestSparseMulDense<43,53>(max_rel_err);
            max_rel_err = TestSparseMulDense<73,11>(max_rel_err);
            max_rel_err = TestSparseMulDense<57,67>(max_rel_err);
            max_rel_err = TestSparseMulDense<89,97>(max_rel_err);
            max_rel_err = TestSparseMulDense<97,17>(max_rel_err);
            max_rel_err = TestSparseMulDense<37,67>(max_rel_err);
            max_rel_err = TestSparseMulDense<33,41>(max_rel_err);
            max_rel_err = TestSparseMulDense<89,47>(max_rel_err);
            max_rel_err = TestSparseMulDense<97,41>(max_rel_err);
        }
        fout << "TestSparseMulDense(): max. relative error: " << max_rel_err << std::endl;
    }

    {
        double max_rel_err = 0.0;
        for (int testNo = 0; testNo < 10; ++testNo) {
            max_rel_err = TestDenseMulSparse<37,67>(max_rel_err);
            max_rel_err = TestDenseMulSparse<43,53>(max_rel_err);
            max_rel_err = TestDenseMulSparse<73,11>(max_rel_err);
            max_rel_err = TestDenseMulSparse<57,67>(max_rel_err);
            max_rel_err = TestDenseMulSparse<89,97>(max_rel_err);
            max_rel_err = TestDenseMulSparse<97,17>(max_rel_err);
            max_rel_err = TestDenseMulSparse<37,67>(max_rel_err);
            max_rel_err = TestDenseMulSparse<33,41>(max_rel_err);
            max_rel_err = TestDenseMulSparse<89,47>(max_rel_err);
            max_rel_err = TestDenseMulSparse<97,41>(max_rel_err);
        }
        fout << "TestDenseMulSparse(): max. relative error: " << max_rel_err << std::endl;
    }

    {
        double max_rel_err = 0.0;
        for (int testNo = 0; testNo < 10; ++testNo) {
            max_rel_err = TestDenseMulSparseTransposed<37,67>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<43,53>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<73,11>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<57,67>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<89,97>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<97,17>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<37,67>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<33,41>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<89,47>(max_rel_err);
            max_rel_err = TestDenseMulSparseTransposed<97,41>(max_rel_err);
        }
        fout << "TestDenseMulSparseTransposed(): max. relative error: " << max_rel_err << std::endl;
    }

    {
        double max_rel_err = 0.0;
        for (int testNo = 0; testNo < 10; ++testNo) {
            max_rel_err = TestDenseMulDense<37,67>(max_rel_err);
            max_rel_err = TestDenseMulDense<43,53>(max_rel_err);
            max_rel_err = TestDenseMulDense<73,11>(max_rel_err);
            max_rel_err = TestDenseMulDense<57,67>(max_rel_err);
            max_rel_err = TestDenseMulDense<89,97>(max_rel_err);
            max_rel_err = TestDenseMulDense<97,17>(max_rel_err);
            max_rel_err = TestDenseMulDense<37,67>(max_rel_err);
            max_rel_err = TestDenseMulDense<33,41>(max_rel_err);
            max_rel_err = TestDenseMulDense<89,47>(max_rel_err);
            max_rel_err = TestDenseMulDense<97,41>(max_rel_err);
        }
        fout << "TestDenseMulDense(): max. relative error: " << max_rel_err << std::endl;
    }

    {
        double max_rel_err = 0.0;
        for (int testNo = 0; testNo < 50; ++testNo) {
            max_rel_err = TestMatrixMulVector<37,67>(max_rel_err);
            max_rel_err = TestMatrixMulVector<43,53>(max_rel_err);
            max_rel_err = TestMatrixMulVector<73,11>(max_rel_err);
            max_rel_err = TestMatrixMulVector<57,67>(max_rel_err);
            max_rel_err = TestMatrixMulVector<89,97>(max_rel_err);
            max_rel_err = TestMatrixMulVector<97,17>(max_rel_err);
            max_rel_err = TestMatrixMulVector<37,67>(max_rel_err);
            max_rel_err = TestMatrixMulVector<33,41>(max_rel_err);
            max_rel_err = TestMatrixMulVector<89,47>(max_rel_err);
            max_rel_err = TestMatrixMulVector<97,41>(max_rel_err);
        }
        fout << "TestMatrixMulVector(): max. relative error: " << max_rel_err << std::endl;
    }
    fout << std::endl << std::flush;
}

