#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------
#include "allscale/api/user/data/grid.h"
#include "amdados/app/static_grid.h"

using namespace allscale::api::user;

template<size_t LENGTH>
using VECTOR = allscale::utils::grid<double, LENGTH>;

template<size_t NROWS, size_t NCOLS>
using MATRIX = allscale::utils::grid<double, NROWS, NCOLS>;

//-------------------------------------------------------------------------------------------------
// Function checks that specified objects are two distinct intances of some class or type.
//-------------------------------------------------------------------------------------------------
template<typename A, typename B>
inline bool CheckDistinctObjects(const A & a, const B & b)
{
    return (static_cast<const void*>(&a) != static_cast<const void*>(&b));
}


//-------------------------------------------------------------------------------------------------
// Matrix multiplication: result = A * B.
// \param  result  out: NROWS-x-NCOLS matrix.
// \param  A       NROWS-x-MSIZE matrix.
// \param  B       MSIZE-x-NCOLS matrix.
//-------------------------------------------------------------------------------------------------

template<size_t NROWS, size_t MSIZE, size_t NCOLS>
void MatMult(allscale::utils::grid<double, NROWS, NCOLS> & result,
       const allscale::utils::grid<double, NROWS, MSIZE> & A,
       const allscale::utils::grid<double, MSIZE, NCOLS> & B)
{
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            double sum = 0.0;
            for (int k = 0; k < static_cast<int>(MSIZE); ++k) {
                sum += A[{r,k}] * B[{k,c}];
            }
            result[{r,c}] = sum;
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Matrix multiplication: result = A * B^t. Matrix B is implicitly transposed.
// \param  result  out: NROWS-x-NCOLS matrix.
// \param  A       NROWS-x-MSIZE matrix.
// \param  B       NCOLS-x-MSIZE matrix.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t MSIZE, size_t NCOLS>
void MatMultTransposed(allscale::utils::grid<double, NROWS, NCOLS> & result,
                 const allscale::utils::grid<double, NROWS, MSIZE> & A,
                 const allscale::utils::grid<double, NCOLS, MSIZE> & B)
{
    assert(static_cast<const void*>(&result) != static_cast<const void*>(&A));
    assert(static_cast<const void*>(&result) != static_cast<const void*>(&B));
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            double sum = 0.0;
            for (int k = 0; k < static_cast<int>(MSIZE); ++k) {
                sum += A[{r,k}] * B[{c,k}];     // note: B is transposed
            }
            result[{r,c}] = sum;
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Matrix-vector multiplication: result = A * v.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void MatVecMult(allscale::utils::grid<double, NROWS>        & result,
          const allscale::utils::grid<double, NROWS, NCOLS> & A,
          const allscale::utils::grid<double, NCOLS>        & v)
{
    assert(static_cast<const void*>(&result) != static_cast<const void*>(&v));
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        double sum = 0.0;
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            sum += A[{r,c}] * v[{c}];
        }
        result[{r}] = sum;
    }
}

//-------------------------------------------------------------------------------------------------
// Add vectors: result = a + b.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
void AddVectors(allscale::utils::grid<double, LENGTH> & result,
          const allscale::utils::grid<double, LENGTH> & a,
          const allscale::utils::grid<double, LENGTH> & b)
{
    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        result[{i}] = a[{i}] + b[{i}];
    }
}

//-------------------------------------------------------------------------------------------------
// Subtract vectors: result = a - b.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
void SubtractVectors(allscale::utils::grid<double, LENGTH> & result,
               const allscale::utils::grid<double, LENGTH> & a,
               const allscale::utils::grid<double, LENGTH> & b)
{
    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        result[{i}] = a[{i}] - b[{i}];
    }
}

//-------------------------------------------------------------------------------------------------
// Add matrices: result = A + B.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void AddMatrices(allscale::utils::grid<double, NROWS, NCOLS> & result,
           const allscale::utils::grid<double, NROWS, NCOLS> & A,
           const allscale::utils::grid<double, NROWS, NCOLS> & B)
{
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            result[{r,c}] = A[{r,c}] + B[{r,c}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Subtract matrices: result = A - B.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void SubtractMatrices(allscale::utils::grid<double, NROWS, NCOLS> & result,
                const allscale::utils::grid<double, NROWS, NCOLS> & A,
                const allscale::utils::grid<double, NROWS, NCOLS> & B)
{
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            result[{r,c}] = A[{r,c}] - B[{r,c}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function initializes the vector by the value specified (default is zero).
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
void FillVector(allscale::utils::grid<double, LENGTH> & v, double vfill = 0.0)
{
    for (int k = 0; k < static_cast<int>(LENGTH); ++k) {
        v[{k}] = vfill;
    }
}

//-------------------------------------------------------------------------------------------------
// Function initializes the matrix by the value specified (default is zero).
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void FillMatrix(allscale::utils::grid<double, NROWS, NCOLS> & A, double vfill = 0.0)
{
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            A[{r,c}] = vfill;
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function initializes the identity matrix.
//-------------------------------------------------------------------------------------------------
template<size_t MSIZE>
void MakeIdentityMatrix(allscale::utils::grid<double, MSIZE, MSIZE> & A)
{
    for (int r = 0; r < static_cast<int>(MSIZE); ++r) {
        for (int c = 0; c < static_cast<int>(MSIZE); ++c) {
            A[{r,c}] = 0.0;
        }
        A[{r,r}] = 1.0;
    }
}

//-------------------------------------------------------------------------------------------------
// Function computes transposed matrix.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void GetTransposed(allscale::utils::grid<double, NCOLS, NROWS> & At,
             const allscale::utils::grid<double, NROWS, NCOLS> & A)
{
    assert(static_cast<const void*>(&At) != static_cast<const void*>(&A));
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            At[{c,r}] = A[{r,c}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Due to round-off errors a matrix supposed to be symmetric can loose this property.
// The function brings the matrix back to symmetry.
//-------------------------------------------------------------------------------------------------
template<size_t MSIZE>
void Symmetrize(allscale::utils::grid<double, MSIZE, MSIZE> & A)
{
    for (int i = 0; i < static_cast<int>(MSIZE); ++i) {
        for (int j = i + 1; j < static_cast<int>(MSIZE); ++j) {
            double v = 0.5 * (A[{j,i}] + A[{i,j}]);
            A[{j,i}] = v;
            A[{i,j}] = v;
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Multiplying by a scalar: A = A * mult.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void MatScalarMult(allscale::utils::grid<double, NROWS, NCOLS> & A, const double mult)
{
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            A[{r,c}] *= mult;
        }
    }
}

//-------------------------------------------------------------------------------------------------
// L2 norm of a vector |a|.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
double NormVec(const allscale::utils::grid<double, LENGTH> & a)
{
    double sqSum = 0.0;
    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        double v = a[{i}];
        sqSum += v * v;
    }
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// L2 norm of vector difference: |a - b|.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
double NormVecDiff(const allscale::utils::grid<double, LENGTH> & a,
                   const allscale::utils::grid<double, LENGTH> & b)
{
    double sqSum = 0.0;
    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        double v = a[{i}] - b[{i}];
        sqSum += v * v;
    }
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// Frobenius norm of a matrix: |A|_fro.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double NormMat(const allscale::utils::grid<double, NROWS, NCOLS> & A)
{
    double sqSum = 0.0;
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            double v = A[{r,c}];
            sqSum += v * v;
        }
    }
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// Frobenius norm of matrix difference: |A - B|_fro.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double NormMatDiff(const allscale::utils::grid<double, NROWS, NCOLS> & A,
                   const allscale::utils::grid<double, NROWS, NCOLS> & B)
{
    double sqSum = 0.0;
    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            double v = A[{r,c}] - B[{r,c}];
            sqSum += v * v;
        }
    }
    return std::sqrt(std::fabs(sqSum));
}


////-------------------------------------------------------------------------------------------------
//// Change vector sign in-place: v = -v.
////-------------------------------------------------------------------------------------------------
//template<size_t LENGTH>
//void NegateVector(allscale::utils::grid<double, LENGTH> & v)
//{
//    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
//        v[{i}] = - v[{i}];
//    }
//}
