#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

namespace amdados {
namespace app {
namespace utils {

#undef  __AMDADOS_BOUND_CHECKING__
#define __AMDADOS_BOUND_CHECKING__

#undef  __AMDADOS_ROW_MAJOR__
//#define __AMDADOS_ROW_MAJOR__

//=================================================================================================
// Base class for vector or matrix (of size 'SIZE') which cannot be instantiated on its
// own but can be used as a light-weight, vector-like view to derived class.
//=================================================================================================
template<int SIZE>
class VectorView
{
protected:
    double * data;                      // content

private:
    VectorView();                       // no default constructor
    VectorView(const VectorView &);     // no copy constructor

protected:
    // Constructor just sets up data pointer.
    inline explicit VectorView(double * p) : data(p) {}

public:
    // Copy operator.
    inline VectorView & operator=(const VectorView & vec) {
        std::copy(vec.data, vec.data + SIZE, data);
        return (*this);
    }

    // Indexing operator provides read only access to vector element.
    inline const double & operator()(int i) const {
#ifndef NDEBUG
#ifdef __AMDADOS_BOUND_CHECKING__
        if (!(static_cast<unsigned int>(i) < static_cast<unsigned int>(SIZE))) assert_true(0);
#endif
#endif
        return data[i];
    }

    // Indexing operator provides read/write access to vector element.
    inline double & operator()(int i) {
#ifndef NDEBUG
#ifdef __AMDADOS_BOUND_CHECKING__
        if (!(static_cast<unsigned int>(i) < static_cast<unsigned int>(SIZE))) assert_true(0);
#endif
#endif
        return data[i];
    }

    // Functions for iterating the vector.
    inline double * begin() { return data; }
    inline double * end()   { return data + SIZE; }

    inline const double * begin() const { return data; }
    inline const double * end()   const { return data + SIZE; }
};

//=================================================================================================
// Simple vector-like wrapper around a raw array of constant size.
// Rationale: (1) safe indexing operations; (2) no stack overflow for automatic objects.
//=================================================================================================
template<int SIZE>
class Vector : public VectorView<SIZE>
{
private: using VectorView<SIZE>::data;
public:
    // Default constructor creates vector filled by zeros.
    inline Vector() : VectorView<SIZE>(new double[SIZE]) {
        std::fill(data, data + SIZE, 0.0);
    }

    // Copy constructor (allocates memory).
    inline Vector(const Vector<SIZE> & vec) : VectorView<SIZE>(new double[SIZE]) {
        std::copy(vec.data, vec.data + SIZE, data);
    }

    // Destructor frees allocated memory.
    inline ~Vector() {
        delete [] data;
    }
};

//=================================================================================================
// Simple matrix-like wrapper around a raw array of constant size.
// Rationale: (1) safe indexing operations; (2) no stack overflow for automatic objects.
// Note, the matrix is row-major, i.e. the column index is faster than the row one.
// Note, implicit casting from matrix to vector view is allowed.
//=================================================================================================
template<int NROWS, int NCOLS>
class Matrix : public VectorView<NROWS * NCOLS>
{
public:  enum { SIZE = NROWS * NCOLS };
private: using VectorView<SIZE>::data;
private: VectorView<SIZE> & operator=(const VectorView<SIZE> &);    // vector copy is disabled
public:
    // Default constructor creates matrix filled by zeros.
    inline Matrix() : VectorView<SIZE>(new double[SIZE]) {
        std::fill(data, data + SIZE, 0.0);
    }

    // Copy constructor (allocates memory). Note, user can copy a matrix but not a vector.
    inline Matrix(const Matrix & mat) : VectorView<SIZE>(new double[SIZE]) {
        std::copy(mat.data, mat.data + SIZE, data);
    }

    // Destructor frees allocated memory.
    inline ~Matrix() {
        delete [] data;
    }

    // Copy operator. Note, user can copy a matrix but not a vector.
    inline Matrix & operator=(const Matrix & mat) {
        std::copy(mat.data, mat.data + SIZE, data);
        return (*this);
    }

    // Indexing operator provides read only access to matrix element.
    inline const double & operator()(int r, int c) const {
#ifndef NDEBUG
#ifdef __AMDADOS_BOUND_CHECKING__
        if (!((static_cast<unsigned int>(r) < static_cast<unsigned int>(NROWS)) &&
              (static_cast<unsigned int>(c) < static_cast<unsigned int>(NCOLS)))) assert_true(0);
#endif
#endif

#ifdef __AMDADOS_ROW_MAJOR__
        return data[r * NCOLS + c];
#else
        return data[r + NROWS * c];
#endif
    }

    // Indexing operator provides read/write access to matrix element.
    inline double & operator()(int r, int c) {
#ifndef NDEBUG
#ifdef __AMDADOS_BOUND_CHECKING__
        if (!((static_cast<unsigned int>(r) < static_cast<unsigned int>(NROWS)) &&
              (static_cast<unsigned int>(c) < static_cast<unsigned int>(NCOLS)))) assert_true(0);
#endif
#endif

#ifdef __AMDADOS_ROW_MAJOR__
        return data[r * NCOLS + c];
#else
        return data[r + NROWS * c];
#endif
    }

    // Function converts 2D index to a flat 1D one.
    static inline int sub2ind(int r, int c) {
#ifndef NDEBUG
#ifdef __AMDADOS_BOUND_CHECKING__
        if (!((static_cast<unsigned int>(r) < static_cast<unsigned int>(NROWS)) &&
              (static_cast<unsigned int>(c) < static_cast<unsigned int>(NCOLS)))) assert_true(0);
#endif
#endif

#ifdef __AMDADOS_ROW_MAJOR__
        return (r * NCOLS + c);
#else
        return (r + NROWS * c);
#endif
    }
};

//-------------------------------------------------------------------------------------------------
// Function checks that specified objects are two distinct instances of some class or type.
//-------------------------------------------------------------------------------------------------
template<typename A, typename B>
inline bool CheckDistinctObjects(const A & a, const B & b)
{
    return (static_cast<const void*>(&a) != static_cast<const void*>(&b));
}

//-------------------------------------------------------------------------------------------------
// Copy vector from Allscale vector.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void VectorFromAllscale(VectorView<SIZE> & v, const allscale::utils::grid<double,SIZE> & a)
{
    for (int i = 0; i < SIZE; ++i) { v(i) = a[{i}]; }
}

//-------------------------------------------------------------------------------------------------
// Copy Allscale vector from vector.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void AllscaleFromVector(allscale::utils::grid<double,SIZE> & a, const VectorView<SIZE> & v)
{
    for (int i = 0; i < SIZE; ++i) { a[{i}] = v(i); }
}

//-------------------------------------------------------------------------------------------------
// Copy matrix from Allscale matrix.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void MatrixFromAllscale(Matrix<NROWS,NCOLS> & m,
                        const allscale::utils::grid<double,NROWS,NCOLS> & a)
{
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) { m(r,c) = a[{r,c}]; }}
}

//-------------------------------------------------------------------------------------------------
// Copy Allscale matrix from matrix.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void AllscaleFromMatrix(allscale::utils::grid<double,NROWS,NCOLS> & a,
                        const Matrix<NROWS,NCOLS> & m)
{
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) { a[{r,c}] = m(r,c); }}
}

//-------------------------------------------------------------------------------------------------
// Matrix multiplication: result = A * B.
// \param  result  out: NROWS-x-NCOLS matrix.
// \param  A       NROWS-x-MSIZE matrix.
// \param  B       MSIZE-x-NCOLS matrix.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int MSIZE, int NCOLS>
void MatMult(      Matrix<NROWS,NCOLS> & result,
             const Matrix<NROWS,MSIZE> & A,
             const Matrix<MSIZE,NCOLS> & B)
{
    assert_true(CheckDistinctObjects(result, A));
    assert_true(CheckDistinctObjects(result, B));
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) {
        double sum = 0.0;
        for (int k = 0; k < MSIZE; ++k) { sum += A(r,k) * B(k,c); }
        result(r,c) = sum;
    }}
}

//-------------------------------------------------------------------------------------------------
// Matrix multiplication with transposition: result = A * B^t.
// Note, the matrix B is not explicitly transposed, rather it is traversed in "transposed" way.
// \param  result  out: NROWS-x-NCOLS matrix.
// \param  A       NROWS-x-MSIZE matrix.
// \param  B       NCOLS-x-MSIZE matrix.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int MSIZE, int NCOLS>
void MatMultTr(      Matrix<NROWS,NCOLS> & result,
               const Matrix<NROWS,MSIZE> & A,
               const Matrix<NCOLS,MSIZE> & B)
{
    assert_true(CheckDistinctObjects(result, A));
    assert_true(CheckDistinctObjects(result, B));
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) {
        double sum = 0.0;
        for (int k = 0; k < MSIZE; ++k) { sum += A(r,k) * B(c,k); } // get B as transposed
        result(r,c) = sum;
    }}
}

//-------------------------------------------------------------------------------------------------
// Matrix-vector multiplication: result = A * v.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void MatVecMult(      VectorView<NROWS>   & result,
                const Matrix<NROWS,NCOLS> & A,
                const VectorView<NCOLS>   & v)
{
    assert_true(CheckDistinctObjects(result, v));
    for (int r = 0; r < NROWS; ++r) {
        double sum = 0.0;
        for (int c = 0; c < NCOLS; ++c) { sum += A(r,c) * v(c); }
        result(r) = sum;
    }
}

//-------------------------------------------------------------------------------------------------
// Add vectors: result = a + b.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void AddVectors(      VectorView<SIZE> & result,
                const VectorView<SIZE> & a,
                const VectorView<SIZE> & b)
{
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<double>());
}

//-------------------------------------------------------------------------------------------------
// Subtract vectors: result = a - b.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void SubtractVectors(      VectorView<SIZE> & result,
                     const VectorView<SIZE> & a,
                     const VectorView<SIZE> & b)
{
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::minus<double>());
}

//-------------------------------------------------------------------------------------------------
// Add matrices: result = A + B.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void AddMatrices(      Matrix<NROWS,NCOLS> & result,
                 const Matrix<NROWS,NCOLS> & A,
                 const Matrix<NROWS,NCOLS> & B)
{
    std::transform(A.begin(), A.end(), B.begin(), result.begin(), std::plus<double>());
}

//-------------------------------------------------------------------------------------------------
// Subtract matrices: result = A - B.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void SubtractMatrices(      Matrix<NROWS,NCOLS> & result,
                      const Matrix<NROWS,NCOLS> & A,
                      const Matrix<NROWS,NCOLS> & B)
{
    std::transform(A.begin(), A.end(), B.begin(), result.begin(), std::minus<double>());
}

//-------------------------------------------------------------------------------------------------
// Function initializes the vector by the value specified (default is zero).
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void FillVector(VectorView<SIZE> & v, double vfill = 0.0)
{
    std::fill(v.begin(), v.end(), vfill);
}

//-------------------------------------------------------------------------------------------------
// Function initializes the matrix by the value specified (default is zero).
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void FillMatrix(Matrix<NROWS,NCOLS> & A, double vfill = 0.0)
{
    std::fill(A.begin(), A.end(), vfill);
}

//-------------------------------------------------------------------------------------------------
// Function initializes the identity matrix.
//-------------------------------------------------------------------------------------------------
template<int MSIZE>
void MakeIdentityMatrix(Matrix<MSIZE,MSIZE> & A)
{
    std::fill(A.begin(), A.end(), 0.0);
    for (int i = 0; i < MSIZE; ++i) { A(i,i) = 1.0; }
}

//-------------------------------------------------------------------------------------------------
// Function computes transposed matrix.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void GetTransposed(      Matrix<NCOLS,NROWS> & At,
                   const Matrix<NROWS,NCOLS> & A)
{
    assert_true(CheckDistinctObjects(At, A));
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) { At(c,r) = A(r,c); }}
}

//-------------------------------------------------------------------------------------------------
// Due to round-off errors a matrix supposed to be symmetric can loose this property.
// The function brings the matrix back to symmetry.
//-------------------------------------------------------------------------------------------------
template<int MSIZE>
void Symmetrize(Matrix<MSIZE,MSIZE> & A)
{
    for (int i = 0;     i < MSIZE; ++i) {
    for (int j = i + 1; j < MSIZE; ++j) { A(j,i) = A(i,j) = 0.5 * (A(j,i) + A(i,j)); }}
}

//-------------------------------------------------------------------------------------------------
// Multiplying vector by a scalar: v = v * mult.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void VecScalarMult(VectorView<SIZE> & v, const double mult)
{
    std::transform(v.begin(), v.end(), v.begin(), [mult](double x) { return x*mult; });
}

//-------------------------------------------------------------------------------------------------
// Multiplying matrix by a scalar: A = A * mult.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void MatScalarMult(Matrix<NROWS,NCOLS> & A, const double mult)
{
    std::transform(A.begin(), A.end(), A.begin(), [mult](double x) { return x*mult; });
}

//-------------------------------------------------------------------------------------------------
// L2 norm of a vector |a|.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
double NormVec(const VectorView<SIZE> & a)
{
    double sqSum = 0.0;
    std::for_each(a.begin(), a.end(), [&](double x) { sqSum += x*x; });
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// L2 norm of vector difference: |a - b|.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
double NormVecDiff(const VectorView<SIZE> & a, const VectorView<SIZE> & b)
{
    double sqSum = 0.0;
    auto ia = a.begin(), ib = b.begin();
    for (; ia != a.end(); ++ia, ++ib) { sqSum += std::pow(*ia - *ib, 2); }
    assert_true(ib == b.end());
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// Frobenius norm of a matrix: |A|_fro.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
double NormMat(const Matrix<NROWS, NCOLS> & A)
{
    double sqSum = 0.0;
    std::for_each(A.begin(), A.end(), [&](double x) { sqSum += x*x; });
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// Frobenius norm of matrix difference: |A - B|_fro.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
double NormMatDiff(const Matrix<NROWS,NCOLS> & A, const Matrix<NROWS,NCOLS> & B)
{
    double sqSum = 0.0;
    auto ia = A.begin(), ib = B.begin();
    for (; ia != A.end(); ++ia, ++ib) { sqSum += std::pow(*ia - *ib, 2); }
    assert_true(ib == B.end());
    return std::sqrt(std::fabs(sqSum));
}

//-------------------------------------------------------------------------------------------------
// Function returns the trace of a square matrix.
//-------------------------------------------------------------------------------------------------
template<int MSIZE>
double Trace(const Matrix<MSIZE,MSIZE> & A)
{
    double sum = 0.0;
    for (int i = 0; i < MSIZE; ++i) { sum += A(i,i); }
    return sum;
}

//-------------------------------------------------------------------------------------------------
// Function generates a random vector with either normal (mu=0, sigma=1) or uniform (0..1)
// distribution of entry values.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void MakeRandomVector(VectorView<SIZE> & v, const char type)
{
    std::mt19937 gen(std::time(nullptr));
    if (type == 'n') {                                              // normal, mu=0, sigma=1
        std::normal_distribution<double> distrib;
        for (int i = 0; i < SIZE; ++i) { v(i) = distrib(gen); }
    } else if (type == 'u') {                                       // uniform, [0..1]
        std::uniform_real_distribution<double> distrib;
        for (int i = 0; i < SIZE; ++i) { v(i) = distrib(gen); }
    } else {
        assert_true(0) << "unknown distribution" << endl;
    }
}

//-------------------------------------------------------------------------------------------------
// Function generates a random matrix uniform (0..1) distribution of entry values.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void MakeRandomMatrix(Matrix<NROWS,NCOLS> & A)
{
    std::mt19937                           gen(std::time(nullptr));
    std::uniform_real_distribution<double> distrib(1.0, 2.0);

    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) { A(r,c) = distrib(gen); }}
}

//-------------------------------------------------------------------------------------------------
// Function checks there is no NAN values among vector entries.
//-------------------------------------------------------------------------------------------------
template<int SIZE>
bool CheckNoNan(const VectorView<SIZE> & v)
{
    for (int i = 0; i < SIZE; i++) { if (std::isnan(v(i))) return false; }
    return true;
}

//-------------------------------------------------------------------------------------------------
// Function checks there is no NAN values among matrix entries.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
bool CheckNoNan(const Matrix<NROWS,NCOLS> & A)
{
    for (int i = 0; i < NROWS; i++) {
    for (int j = 0; j < NCOLS; j++) { if (std::isnan(A(i,j))) return false; }}
    return true;
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<int SIZE>
void WriteVector(const Vector<SIZE> & vec, const string & filename, int precision = 4)
{
    std::ofstream f(filename, std::ios::out | std::ios::trunc);
    assert_true(f.good());
    precision = std::min(std::max(precision, 2), 16);
    for (int i = 0; i < SIZE; i++) { f << std::setprecision(precision) << vec(i) << " "; }
    f << std::endl << std::flush;
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void WriteMatrix(const Matrix<NROWS,NCOLS> & A, const string & filename, int precision = 4)
{
    std::ofstream f(filename, std::ios::out | std::ios::trunc);
    assert_true(f.good());
    precision = std::min(std::max(precision, 2), 16);
    for (int r = 0; r < NROWS; r++) {
    for (int c = 0; c < NCOLS; c++) {
        f << std::setprecision(precision) << A(r,c) << " "; } f << std::endl; }
    f << std::flush;
}

} // end namespace utils
} // end namespace app
} // end namespace amdados

////-------------------------------------------------------------------------------------------------
//// Change vector sign in-place: v = -v.
////-------------------------------------------------------------------------------------------------
//template<int SIZE>
//void NegateVector(VectorView<SIZE> & v)
//{
//    for (int i = 0; i < SIZE; ++i) {
//        v(i) = - v(i);
//    }
//}
