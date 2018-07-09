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
#include "amdados/app/amdados_utils.h"
#include "amdados/app/matrix.h"
#define ARMA_USE_CXX11
#define ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_LAPACK
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_NEWARP
#define ARMA_DONT_USE_ARPACK
#define ARMA_DONT_USE_SUPERLU
#define ARMA_DONT_USE_HDF5
#define ARMA_DONT_USE_OPENMP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
	#include "../../../api/armadillo/include/armadillo"
#pragma GCC diagnostic pop

using namespace ::amdados;

// Tolerance on relative error.
const double TOL = std::sqrt(std::numeric_limits<double>::epsilon());

const int PRIMES[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                      47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103};
const int NumPrimes = static_cast<int>(sizeof(PRIMES)/sizeof(PRIMES[0]));

//-----------------------------------------------------------------------------
// Function copies an "allscale" matrix to "armadillo" one.
//-----------------------------------------------------------------------------
void CopyToArma(arma::mat & dest, const Matrix & source)
{
    dest.zeros(source.NRows(), source.NCols());
    for (int r = 0; r < source.NRows(); ++r) {
    for (int c = 0; c < source.NCols(); ++c) { dest(r,c) = source(r,c); }}
}

//-----------------------------------------------------------------------------
// Function copies an "allscale" vector to "armadillo" one.
//-----------------------------------------------------------------------------
void CopyToArma(arma::vec & dest, const Vector & source)
{
    dest.zeros(source.Size());
    for (int k = 0; k < source.Size(); ++k) { dest(k) = source(k); }
}

//-----------------------------------------------------------------------------
// Function compares matrices for equality and updates the max. relative error.
//-----------------------------------------------------------------------------
void RelativeErrorCheckAndUpdate(const arma::mat & A, const arma::mat & B,
                                 double & max_rel_err, const char * what)
{
    double rel_err = arma::norm(A - B, "fro") / arma::norm(A, "fro");
    max_rel_err = std::max(max_rel_err, rel_err);
    EXPECT_LT(rel_err, TOL) << "Relative error of "
                            << what << " exceeded tolerance";
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void TestDenseMulDense(double & max_rel_err, int nrows, int ncols)
{
    const int msize = 167;

    Matrix    A1(nrows, msize);
    Matrix    A2(msize, ncols);
    Matrix    A_mul(nrows, ncols);
    arma::mat B1, B2, B_mul, A_tmp;

    // Create random matrices.
    MakeRandom(A1, 'u');
    MakeRandom(A2, 'u');
    CopyToArma(B1, A1);
    CopyToArma(B2, A2);

    // Check multiplication.
    MatMult(A_mul, A1, A2);         // A_mul = A1 * A2
    B_mul = B1 * B2;                // B_mul = B1 * B2
    CopyToArma(A_tmp, A_mul);       // A_tmp = A_mul

    RelativeErrorCheckAndUpdate(B_mul, A_tmp, max_rel_err,
            "dense*dense multiplication");
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void TestMatrixMulVector(double & max_rel_err, int nrows, int ncols)
{
    Vector result(nrows);
    Vector source(ncols);
    arma::vec res1, res2, src;
    MakeRandom(source, 'u');
    CopyToArma(src, source);

    Matrix A(nrows, ncols);
    arma::mat B;
    MakeRandom(A, 'u');
    CopyToArma(B, A);

    MatVecMult(result, A, source);
    CopyToArma(res1, result);
    res2 = B * src;
    RelativeErrorCheckAndUpdate(res1, res2, max_rel_err,
            "dense*vector multiplication");
}

//-----------------------------------------------------------------------------
// Function tests the matrix library written by means of Allscale API.
//-----------------------------------------------------------------------------
TEST(MatrixTests, Basic)
{
    // Open the output log-file.
    std::string fname = "matrix_test.log";
    std::fstream log_file(fname, std::ios::out | std::ios::trunc);
    assert_true(log_file.good()) << "failed to open: " << fname << std::endl;

    // Test dense matrix by dense matrix multiplication.
    double max_rel_err = 0.0;
    for (int i = 0; i < NumPrimes; ++i) {
    for (int k = 0; k < NumPrimes; ++k) {
        TestDenseMulDense(max_rel_err, PRIMES[i], PRIMES[k]);
    }}
    log_file << "TestDenseMulDense(): max. relative error: "
             << max_rel_err << std::endl;

    // Test dense matrix by vector multiplication.
    max_rel_err = 0.0;
    for (int i = 0; i < NumPrimes; ++i) {
    for (int k = 0; k < NumPrimes; ++k) {
        TestMatrixMulVector(max_rel_err, PRIMES[i], PRIMES[k]);
    }}
    log_file << "TestMatrixMulVector(): max. relative error: "
             << max_rel_err << std::endl << std::endl << std::flush;
}

