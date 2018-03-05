//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

// Macros check vector/matrix sizes and produce error message with
// file/line information upon mismatch.
#ifndef NDEBUG
#if defined(VEC_CHECK_SIZE) || defined(MAT_CHECK_SIZE)
#error Macro redefinifion
#endif

#define VEC_CHECK_SIZE(vec, length) { \
    if (vec.Size() != length) assert_true(0) << "size mismatch" << std::endl; }

#define MAT_CHECK_SIZE(mat, nrows, ncols) { \
    if ((mat.NRows() != nrows) || (mat.NCols() != ncols)) \
        assert_true(0) << "size mismatch" << std::endl; }
#endif

//=============================================================================
// Base class for any vector or matrix which cannot be instantiated on its own
// but can be used as a light-weight, vector-like view to derived class.
//=============================================================================
class VectorView
{
protected:
    double * data;                      // content
    int      size;                      // vector length

protected:
    VectorView() : data(nullptr), size(0) {}
    VectorView(const VectorView &) = delete;

public:
    // Copy operator expects preallocated object as memory
    // (re)allocation is not allowed for viewer class.
    VectorView & operator=(const VectorView & vec)
    {
        if (IsDistinct(vec)) {
            assert_true(SameSize(vec));
            std::copy(vec.begin(), vec.end(), begin());
        }
        return (*this);
    }

    // Indexing operator provides read only access to vector element.
    const double & operator()(int i) const
    {
#ifndef NDEBUG
        if (!(static_cast<unsigned>(i) < static_cast<unsigned>(size)))
            assert_true(0);
#endif
        return data[i];
    }

    // Indexing operator provides read/write access to vector element.
    double & operator()(int i)
    {
#ifndef NDEBUG
        if (!(static_cast<unsigned>(i) < static_cast<unsigned>(size)))
            assert_true(0);
#endif
        return data[i];
    }

    // Functions for iterating the (constant) vector.
    double * begin() { return data; }
    double * end()   { return data + size; }

    const double * begin() const { return data; }
    const double * end()   const { return data + size; }

    // Function returns the size of this vector (view).
    int Size() const { return size; }

    // Function returns "true" if both vectors have the same length.
    bool SameSize(const VectorView & v) const { return (size == v.size); }

    // Function returns "true" if specified vector-view does not
    // refer to the same object.
    bool IsDistinct(const VectorView & v) const { return (data != v.data); }

    // Returns "true" if this vector-view refers to an empty object.
    bool Empty() const { return (data == nullptr); }
};

//=============================================================================
// Simple vector-like wrapper around a raw array of constant size.
//=============================================================================
class Vector : public VectorView
{
private:
    using VectorView::data;
    using VectorView::size;

public:
    // Default constructor.
    Vector() : VectorView()
    {
        data = nullptr;
        size = 0;
    }

    // Default constructor creates vector filled by zeros.
    Vector(int vec_size) : VectorView()
    {
        data = nullptr;
        size = 0;
        Resize(vec_size);
    }

    // Copy constructor (allocates memory).
    Vector(const Vector & vec) : VectorView()
    {
        data = nullptr;
        size = 0;
        Resize(vec.size, false);
        std::copy(vec.begin(), vec.end(), begin());
    }

    // Deallocates and clears this object.
    void Clear()
    {
        if (data != nullptr) { delete [] data; }
        data = nullptr;
        size = 0;
    }

    // Resizes this vector and optionally fills it up by zeros.
    void Resize(int new_size, bool fillzero = true)
    {
        if (size != new_size) {
            Clear();
            if (new_size > 0) {
                data = new double[new_size];
                size = new_size;
            }
        }
        if (fillzero) std::fill(begin(), end(), 0.0);
    }

    // Destructor frees allocated memory.
    ~Vector()
    {
        Clear();
    }

    // Copy operator reallocates vector if needed.
    VectorView & operator=(const VectorView & vec)
    {
        if (IsDistinct(vec)) {
            if (!SameSize(vec)) {
                Resize(vec.Size(), false);
            }
            std::copy(vec.begin(), vec.end(), begin());
        }
        return (*this);
    }
};

//=============================================================================
// Simple matrix-like wrapper around a raw array of constant size.
// Note, the matrix is row-major (the column index is faster than the row one).
// Note, implicit casting from matrix to vector-view is allowed.
//=============================================================================
class Matrix : public VectorView
{
private:
    using VectorView::data;
    using VectorView::size;
    int nrows;                  ///< number of rows
    int ncols;                  ///< number of columns

public:
    // Default constructor.
    Matrix() : VectorView(), nrows(0), ncols(0)
    {
        data = nullptr;
        size = 0;
    }

    // Default constructor creates matrix filled by zeros.
    Matrix(int numrows, int numcols) : VectorView(), nrows(0), ncols(0)
    {
        data = nullptr;
        size = 0;
        Resize(numrows, numcols);
    }

    // Copy constructor. Note, user can copy a matrix but not a vector.
    Matrix(const Matrix & mat) : VectorView(), nrows(0), ncols(0)
    {
        data = nullptr;
        size = 0;
        Resize(mat.nrows, mat.ncols, false);
        std::copy(mat.begin(), mat.end(), begin());
    }

    // Deallocates and clears this object.
    void Clear()
    {
        if (data != nullptr) { delete [] data; }
        data = nullptr;
        size = nrows = ncols = 0;
    }

    // Resizes this matrix and optionally fills it up by zeros.
    void Resize(int numrows, int numcols, bool fillzero = true)
    {
        if ((nrows != numrows) || (ncols != numcols)) {
            Clear();
            if ((numrows > 0) && (numcols > 0)) {
                size = numrows * numcols;
                nrows = numrows;
                ncols = numcols;
                assert_true(size > 0);
                data = new double[size];
            }
        }
        if (fillzero) std::fill(begin(), end(), 0.0);
    }

    // Destructor frees allocated memory.
    ~Matrix()
    {
        Clear();
    }

    // Copy operator reallocates matrix if needed.
    Matrix & operator=(const Matrix & mat)
    {
        if (IsDistinct(mat)) {
            if (!SameSize(mat)) {
                Resize(mat.nrows, mat.ncols, false);
            }
            std::copy(mat.begin(), mat.end(), begin());
        }
        return (*this);
    }

    // Note, user can copy a matrix but not a vector.
    VectorView & operator=(const VectorView &) = delete;

    // Indexing operator provides read only access to matrix element.
    const double & operator()(int r, int c) const
    {
#ifndef NDEBUG
        if (!((static_cast<unsigned>(r) < static_cast<unsigned>(nrows)) &&
              (static_cast<unsigned>(c) < static_cast<unsigned>(ncols))))
            assert_true(0);
#endif
        return data[r * ncols + c];
    }

    // Indexing operator provides read/write access to matrix element.
    double & operator()(int r, int c)
    {
#ifndef NDEBUG
        if (!((static_cast<unsigned>(r) < static_cast<unsigned>(nrows)) &&
              (static_cast<unsigned>(c) < static_cast<unsigned>(ncols))))
            assert_true(0);
#endif
        return data[r * ncols + c];
    }

    // Functions return the number of rows, the number of columns and
    // the total number of matrix elements respectively.
    int NRows() const { return nrows; }
    int NCols() const { return ncols; }
    int Size()  const { return size;  }

    // Returns "true" if both matrices have the same size.
    bool SameSize(const Matrix & mat) const
    {
        return ((nrows == mat.nrows) && (ncols == mat.ncols));
    }

    // Returns "true" if this matrix has the same size as transposed(mat).
    bool SameSizeTr(const Matrix & mat) const
    {
        return ((nrows == mat.ncols) && (ncols == mat.nrows));
    }

    // Returns "true" if matrix is square.
    bool IsSquare() const { return (nrows == ncols); }
};

////-----------------------------------------------------------------------------
//// Copy vector from Allscale vector.
////-----------------------------------------------------------------------------
//template<typename AllscaleGrid>
//void VectorFromAllscale(VectorView & v, const AllscaleGrid & a)
//{
//    const int size = a.size()[0];
//    assert_true(v.Size() == size);
//    for (int i = 0; i < size; ++i) { v(i) = a[{i}]; }
//}
//
////-----------------------------------------------------------------------------
//// Copy Allscale vector from vector.
////-----------------------------------------------------------------------------
//template<typename AllscaleGrid>
//void AllscaleFromVector(AllscaleGrid & a, const VectorView & v)
//{
//    const int size = a.size()[0];
//    assert_true(v.Size() == size);
//    for (int i = 0; i < size; ++i) { a[{i}] = v(i); }
//}

////-----------------------------------------------------------------------------
//// Copy matrix from Allscale matrix.
////-----------------------------------------------------------------------------
//template<typename AllscaleGrid>
//void MatrixFromAllscale(Matrix & m, const AllscaleGrid & a)
//{
//    const int nrows = a.size()[0];
//    const int ncols = a.size()[1];
//    assert_true((m.NRows() == nrows) && (m.NCols() == ncols));
//    for (int r = 0; r < nrows; ++r) {
//    for (int c = 0; c < ncols; ++c) { m(r,c) = a[{r,c}]; }}
//}
//
////-----------------------------------------------------------------------------
//// Copy Allscale matrix from matrix.
////-----------------------------------------------------------------------------
//template<typename AllscaleGrid>
//void AllscaleFromMatrix(AllscaleGrid & a, const Matrix & m)
//{
//    const int nrows = a.size()[0];
//    const int ncols = a.size()[1];
//    assert_true((m.NRows() == nrows) && (m.NCols() == ncols));
//    for (int r = 0; r < nrows; ++r) {
//    for (int c = 0; c < ncols; ++c) { a[{r,c}] = m(r,c); }}
//}

void MatMult(Matrix & result, const Matrix & A, const Matrix & B);

void MatMultTr(Matrix & result, const Matrix & A, const Matrix & B);

void MatVecMult(VectorView & result, const Matrix & A, const VectorView & v);

void AddVectors(VectorView & result,
                const VectorView & a, const VectorView & b);

void SubtractVectors(VectorView & result,
                     const VectorView & a, const VectorView & b);

void AddMatrices(Matrix & result,
                 const Matrix & A, const Matrix & B);

void SubtractMatrices(Matrix & result,
                             const Matrix & A, const Matrix & B);

void Fill(VectorView & v, double vfill = 0.0);

void MakeIdentityMatrix(Matrix & A);

void GetTransposed(Matrix & At, const Matrix & A);

void Symmetrize(Matrix & A);

void ScalarMult(VectorView & v, const double mult);

double Norm(const VectorView & v);

double NormDiff(const VectorView & a, const VectorView & b);

double Trace(const Matrix & A);

void Negate(VectorView & v);

void MakeRandom(VectorView & v, const char type);

bool CheckNoNan(const VectorView & v);

void WriteVector(const Vector & vec,
                        const std::string & filename, int precision = 4);

void WriteMatrix(const Matrix & A,
                        const std::string & filename, int precision = 4);

} // namespace amdados

