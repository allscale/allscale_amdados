//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#include "allscale/utils/assert.h"
#include <algorithm>
#include <random>
#include "../include/amdados_utils.h"
#include "../include/matrix.h"

namespace amdados {

//-----------------------------------------------------------------------------
// Matrix multiplication: result = A * B.
// @param  result  out: nrows-x-ncols matrix.
// @param  A       nrows-x-msize matrix.
// @param  B       msize-x-ncols matrix.
//-----------------------------------------------------------------------------
void MatMult(Matrix & result, const Matrix & A, const Matrix & B)
{
    const int nrows = A.NRows();
    const int msize = A.NCols();
    const int ncols = B.NCols();
    assert_true((result.NRows() == nrows) &&
                (result.NCols() == ncols) && (msize == B.NRows()));
    for (int r = 0; r < nrows; ++r) {
    for (int c = 0; c < ncols; ++c) {
        double sum = 0.0;
        for (int k = 0; k < msize; ++k) { sum += A(r,k) * B(k,c); }
        result(r,c) = sum;
    }}
}

//-----------------------------------------------------------------------------
// Matrix multiplication with transposition: result = A * B^t.
// Note, the matrix B is only logically but not explicitly transposed.
// @param  result  out: nrows-x-ncols matrix.
// @param  A       nrows-x-msize matrix.
// @param  B       ncols-x-msize matrix.
//-----------------------------------------------------------------------------
void MatMultTr(Matrix & result, const Matrix & A, const Matrix & B)
{
    const int nrows = A.NRows();
    const int ncols = B.NRows();
    const int msize = B.NCols();
    assert_true((result.NRows() == nrows) &&
                (result.NCols() == ncols) && (A.NCols() == msize));
    for (int r = 0; r < nrows; ++r) {
    for (int c = 0; c < ncols; ++c) {     // get B as transposed
        double sum = 0.0;
        for (int k = 0; k < msize; ++k) { sum += A(r,k) * B(c,k); }
        result(r,c) = sum;
    }}
}

//-----------------------------------------------------------------------------
// Matrix-vector multiplication: result = A * v.
//-----------------------------------------------------------------------------
void MatVecMult(VectorView & result, const Matrix & A, const VectorView & v)
{
    const int nrows = A.NRows();
    const int ncols = A.NCols();
    assert_true((result.Size() == nrows) && (v.Size() == ncols));
    for (int r = 0; r < nrows; ++r) {
        double sum = 0.0;
        for (int c = 0; c < ncols; ++c) { sum += A(r,c) * v(c); }
        result(r) = sum;
    }
}

//-----------------------------------------------------------------------------
// Add vectors: result = a + b.
//-----------------------------------------------------------------------------
void AddVectors(VectorView & result,
                const VectorView & a, const VectorView & b)
{
    assert_true(result.SameSize(a) && result.SameSize(b));
    std::transform(a.begin(), a.end(), b.begin(), result.begin(),
                   std::plus<double>());
}

//-----------------------------------------------------------------------------
// Subtract vectors: result = a - b.
//-----------------------------------------------------------------------------
void SubtractVectors(VectorView & result,
                     const VectorView & a, const VectorView & b)
{
    assert_true(result.SameSize(a) && result.SameSize(b));
    std::transform(a.begin(), a.end(), b.begin(), result.begin(),
                   std::minus<double>());
}

//-----------------------------------------------------------------------------
// Add matrices: result = A + B.
//-----------------------------------------------------------------------------
void AddMatrices(Matrix & result, const Matrix & A, const Matrix & B)
{
    assert_true(result.SameSize(A) && result.SameSize(B));
    std::transform(A.begin(), A.end(), B.begin(), result.begin(),
                   std::plus<double>());
}

//-----------------------------------------------------------------------------
// Subtract matrices: result = A - B.
//-----------------------------------------------------------------------------
void SubtractMatrices(Matrix & result, const Matrix & A, const Matrix & B)
{
    assert_true(result.SameSize(A) && result.SameSize(B));
    std::transform(A.begin(), A.end(), B.begin(), result.begin(),
                   std::minus<double>());
}

//-----------------------------------------------------------------------------
// Function initializes the object by the value specified (default is zero).
//-----------------------------------------------------------------------------
void Fill(VectorView & v, double vfill)
{
    std::fill(v.begin(), v.end(), vfill);
}

//-----------------------------------------------------------------------------
// Function initializes the identity matrix.
//-----------------------------------------------------------------------------
void MakeIdentityMatrix(Matrix & A)
{
    std::fill(A.begin(), A.end(), 0.0);
    for (int i = 0; i < A.NRows(); ++i) { A(i,i) = 1.0; }
}

//-----------------------------------------------------------------------------
// Function computes transposed matrix.
//-----------------------------------------------------------------------------
void GetTransposed(Matrix & At, const Matrix & A)
{
    const int nrows = A.NRows();
    const int ncols = A.NCols();
	assert_true(A.SameSizeTr(At));
    for (int r = 0; r < nrows; ++r) {
    for (int c = 0; c < ncols; ++c) { At(c,r) = A(r,c); }}
}

//-----------------------------------------------------------------------------
// Due to round-off errors a matrix supposed to be symmetric can loose this
// property. The function brings the matrix back to symmetry.
//-----------------------------------------------------------------------------
void Symmetrize(Matrix & A)
{
    const int nrows = A.NRows();
    assert_true(A.IsSquare());
    for (int i = 0;     i < nrows; ++i) {
    for (int j = i + 1; j < nrows; ++j) {
        A(j,i) = A(i,j) = 0.5 * (A(j,i) + A(i,j));
    }}
}

//-----------------------------------------------------------------------------
// Multiplying object by a scalar: v = v * mult.
//-----------------------------------------------------------------------------
void ScalarMult(VectorView & v, const double mult)
{
    std::transform(v.begin(), v.end(), v.begin(),
                    [mult](double x) { return x*mult; });
}

//-----------------------------------------------------------------------------
// L2 norm of an object |a|; Frobenius norm for matrices.
//-----------------------------------------------------------------------------
double Norm(const VectorView & v)
{
    double sum = 0.0;
    std::for_each(v.begin(), v.end(), [&](double x) { sum += x*x; });
    return std::sqrt(std::fabs(sum));
}

//-----------------------------------------------------------------------------
// L2 norm of objects difference: |a - b|; Frobenius norm for matrices.
//-----------------------------------------------------------------------------
double NormDiff(const VectorView & a, const VectorView & b)
{
    assert_true(a.SameSize(b));     // weak verification in case of matrices
    double sum = 0.0;
    auto ia = a.begin(), ib = b.begin();
    for (; ia != a.end(); ++ia, ++ib) { sum += std::pow(*ia - *ib, 2); }
    assert_true(ib == b.end());
    return std::sqrt(std::fabs(sum));
}

//-----------------------------------------------------------------------------
// Function returns the trace of a square matrix.
//-----------------------------------------------------------------------------
double Trace(const Matrix & A)
{
    const int nrows = A.NRows();
    assert_true(A.IsSquare());
    double sum = 0.0;
    for (int i = 0; i < nrows; ++i) { sum += A(i,i); }
    return sum;
}

//-----------------------------------------------------------------------------
// Change vector/matrix sign in-place: v = -v.
//-----------------------------------------------------------------------------
void Negate(VectorView & v)
{
    const int size = v.Size();
    for (int i = 0; i < size; ++i) { v(i) = - v(i); }
}

//-----------------------------------------------------------------------------
// Function generates a random object with either normal (mu=0, sigma=1) or
// uniform (0..1) distribution of entry values.
//-----------------------------------------------------------------------------
void MakeRandom(VectorView & v, const char type)
{
    const int size = v.Size();
    std::mt19937_64 gen(RandomSeed());
    if (type == 'n') {                              // normal, mu=0, sigma=1
        std::normal_distribution<double> distrib;
        for (int i = 0; i < size; ++i) { v(i) = distrib(gen); }
    } else if (type == 'u') {                       // uniform, [0..1]
        std::uniform_real_distribution<double> distrib;
        for (int i = 0; i < size; ++i) { v(i) = distrib(gen); }
    } else {
        assert_true(0) << "unknown distribution" << std::endl;
    }
}

//-----------------------------------------------------------------------------
// Function checks there is no NAN values among vector/matrix entries.
//-----------------------------------------------------------------------------
bool CheckNoNan(const VectorView & v)
{
    for (auto i = v.begin(); i != v.end(); ++i) {
        if (std::isnan(*i))
            return false;
    }
    return true;
}

} // namespace amdados

