#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

namespace amdados {
namespace app {
namespace utils {

template<size_t LENGTH>
    using Vector = allscale::utils::grid<double, LENGTH>;

template<size_t NROWS, size_t NCOLS>
    using Matrix = allscale::utils::grid<double, NROWS, NCOLS>;

template<size_t LENGTH>
    using VecPtr = std::unique_ptr< allscale::utils::grid<double, LENGTH> >;

template<size_t NROWS, size_t NCOLS>
    using MatPtr = std::unique_ptr< allscale::utils::grid<double, NROWS, NCOLS> >;

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
void MatMult(      Matrix<NROWS, NCOLS> & result,
             const Matrix<NROWS, MSIZE> & A,
             const Matrix<MSIZE, NCOLS> & B)
{
    assert(CheckDistinctObjects(result, A));
    assert(CheckDistinctObjects(result, B));

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
void MatMultTransposed(      Matrix<NROWS, NCOLS> & result,
                       const Matrix<NROWS, MSIZE> & A,
                       const Matrix<NCOLS, MSIZE> & B)
{
    assert(CheckDistinctObjects(result, A));
    assert(CheckDistinctObjects(result, B));

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
void MatVecMult(      Vector<NROWS>        & result,
                const Matrix<NROWS, NCOLS> & A,
                const Vector<NCOLS>        & v)
{
    assert(CheckDistinctObjects(result, v));

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
void AddVectors(      Vector<LENGTH> & result,
                const Vector<LENGTH> & a,
                const Vector<LENGTH> & b)
{
    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        result[{i}] = a[{i}] + b[{i}];
    }
}

//-------------------------------------------------------------------------------------------------
// Subtract vectors: result = a - b.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
void SubtractVectors(      Vector<LENGTH> & result,
                     const Vector<LENGTH> & a,
                     const Vector<LENGTH> & b)
{
    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        result[{i}] = a[{i}] - b[{i}];
    }
}

//-------------------------------------------------------------------------------------------------
// Add matrices: result = A + B.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void AddMatrices(      Matrix<NROWS, NCOLS> & result,
                 const Matrix<NROWS, NCOLS> & A,
                 const Matrix<NROWS, NCOLS> & B)
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
void SubtractMatrices(      Matrix<NROWS, NCOLS> & result,
                      const Matrix<NROWS, NCOLS> & A,
                      const Matrix<NROWS, NCOLS> & B)
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
void FillVector(Vector<LENGTH> & v, double vfill = 0.0)
{
    for (int k = 0; k < static_cast<int>(LENGTH); ++k) {
        v[{k}] = vfill;
    }
}

//-------------------------------------------------------------------------------------------------
// Function initializes the matrix by the value specified (default is zero).
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void FillMatrix(Matrix<NROWS, NCOLS> & A, double vfill = 0.0)
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
void MakeIdentityMatrix(Matrix<MSIZE, MSIZE> & A)
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
void GetTransposed(      Matrix<NCOLS, NROWS> & At,
                   const Matrix<NROWS, NCOLS> & A)
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
void Symmetrize(Matrix<MSIZE, MSIZE> & A)
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
void MatScalarMult(Matrix<NROWS, NCOLS> & A, const double mult)
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
double NormVec(const Vector<LENGTH> & a)
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
double NormVecDiff(const Vector<LENGTH> & a, const Vector<LENGTH> & b)
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
double NormMat(const Matrix<NROWS, NCOLS> & A)
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
double NormMatDiff(const Matrix<NROWS, NCOLS> & A, const Matrix<NROWS, NCOLS> & B)
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

//-------------------------------------------------------------------------------------------------
// Function returns the trace of a square matrix.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
double Trace(const Matrix<NROWS, NCOLS> & A)
{
    static_assert(NROWS == NCOLS, "square matrix is expected");
    double sum = 0.0;
    for (int i = 0; i < static_cast<int>(NROWS); ++i) {
        sum += A[{i,i}];
    }
    return sum;
}

//-------------------------------------------------------------------------------------------------
// Function generates a random vector.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
void MakeRandomVector(Vector<LENGTH> & v)
{
    std::mt19937                           gen(std::time(nullptr));
    std::uniform_real_distribution<double> distrib(1.0, 2.0);

    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
        v[{i}] = distrib(gen);
    }
}

//-------------------------------------------------------------------------------------------------
// Function generates a random matrix.
//-------------------------------------------------------------------------------------------------
template<size_t NROWS, size_t NCOLS>
void MakeRandomMatrix(Matrix<NROWS,NCOLS> & A)
{
    std::mt19937                           gen(std::time(nullptr));
    std::uniform_real_distribution<double> distrib(1.0, 2.0);

    for (int r = 0; r < static_cast<int>(NROWS); ++r) {
        for (int c = 0; c < static_cast<int>(NCOLS); ++c) {
            A[{r,c}] = distrib(gen);
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function checks there is no NAN values on the grid.
//-------------------------------------------------------------------------------------------------
template<size_t LENGTH>
bool CheckNoNan(const Vector<LENGTH> & grid)
{
    // TODO: check size: must be LENGTH
    for (int i = 0; i < static_cast<int>(LENGTH); i++) {
        if (std::isnan(grid[{i}]))
            return false;
    }
    return true;
}

//-------------------------------------------------------------------------------------------------
// Function checks there is no NAN values on the grid.
//-------------------------------------------------------------------------------------------------
template<size_t NELEMS_X, size_t NELEMS_Y>
bool CheckNoNan(const Matrix<NELEMS_X,NELEMS_Y> & grid)
{
    // TODO: check grid sizes: must be NELEMS_X by NELEMS_Y
    for (int i = 0; i < static_cast<int>(NELEMS_X); i++) {
        for (int j = 0; j < static_cast<int>(NELEMS_Y); j++) {
            if (std::isnan(grid[{i,j}]))
                return false;
        }
    }
    return true;
}

} // end namespace utils
} // end namespace app
} // end namespace amdados

////-------------------------------------------------------------------------------------------------
//// Change vector sign in-place: v = -v.
////-------------------------------------------------------------------------------------------------
//template<size_t LENGTH>
//void NegateVector(Vector<LENGTH> & v)
//{
//    for (int i = 0; i < static_cast<int>(LENGTH); ++i) {
//        v[{i}] = - v[{i}];
//    }
//}
