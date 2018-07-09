//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

#include <vector>

#include <allscale/utils/serializer.h>

namespace amdados {

	using index_t = int64_t;

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
// VectorView class, which is instantiated as a Vector and also serves as a
// base class for Matrix
//=============================================================================
class VectorView
{
protected:

	// content of this vector, memory management done by STL.
    std::vector<double> data;

public:

	// default constructor.
    VectorView() : data() {}

	// constructor, ensures data can fit size elements.
	VectorView(index_t size) : data(size) {}

	VectorView(const VectorView &) = default;

	VectorView(VectorView&& other) {
		data.swap(other.data);
	}

	~VectorView() = default;

    // Indexing operator provides read only access to vector elements.
	// Using a call operator here for unified element access with
	// multi-dimensional indices (e.g. (i,y)) in subclasses.
	const double& operator()(index_t i) const
    {
#ifndef NDEBUG
        if (!(static_cast<unsigned>(i) < static_cast<unsigned>(data.size())))
            assert_true(0);
#endif
        return data[i];
    }

    // Indexing operator provides read/write access to vector elements.
	// Using a call operator here for unified element access with
	// multi-dimensional indices (e.g. (i,y)) in subclasses.
	double& operator()(index_t i)
    {
#ifndef NDEBUG
        if (!(static_cast<unsigned>(i) < static_cast<unsigned>(data.size())))
            assert_true(0);
#endif
        return data[i];
    }

    // Functions for iterating the (constant) vector.
    double * begin() { return &*data.begin(); }
    double * end()   { return &*data.end(); }

    const double * begin() const { return &*data.cbegin(); }
    const double * end()   const { return &*data.cend(); }

    // Function returns the size of this vector (view).
	index_t Size() const { return data.size(); }

    // Function returns "true" if both vectors have the same length.
	// TODO: Quite expensive operation, maybe replace with something faster.
    bool SameSize(const VectorView & v) const { return (data.size() == v.data.size()); }

    // Returns "true" if this vector-view refers to an empty object.
    bool Empty() const { return (data.empty()); }

	// Deletes the content of this VectorView.
	void Clear() {
		data.clear();
	}

	// Resizes the content to fit the given size, optionally also clearing the content.
	void Resize(index_t new_size, bool fillzero = true) {
		if(fillzero) {
			data.clear();
		}
		data.resize(new_size);
	}

	VectorView & operator=(const VectorView & vec)
	{
		data = vec.data;
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& out, const VectorView& vv) {
		out << "VectorView [ ";
		out << "1x" << vv.Size();
		for(const auto& e : vv.data) {
			out << ", " << e;
		}
		out << " ]" << std::endl;
		return out;
	}

};

//=============================================================================
// Simple vector-like wrapper around a raw array of constant size.
//=============================================================================
using Vector = VectorView;

//=============================================================================
// Simple matrix-like wrapper around a raw array of constant size.
// Note, the matrix is row-major (the column index is faster than the row one).
// Note, implicit casting from matrix to vector-view is allowed.
//=============================================================================
class Matrix : public VectorView
{
private:
	index_t nrows; ///< number of rows
    index_t ncols; ///< number of columns

public:
    // Default constructor.
    Matrix() : VectorView(), nrows(0), ncols(0) { }

    // Default constructor creates matrix filled by zeros.
    Matrix(index_t numrows, index_t numcols) : VectorView(numrows*numcols), nrows(numrows), ncols(numcols) { }

    // Copy constructor. Note, user can copy a matrix but not a vector.
    Matrix(const Matrix & mat) : VectorView(mat), nrows(mat.nrows), ncols(mat.ncols) { }

	// Move constructor.
	Matrix(Matrix && mat) : VectorView(std::move(mat)), nrows(mat.nrows), ncols(mat.ncols)
	{
		mat.nrows = 0;
		mat.ncols = 0;
	}

	// Deallocates and clears this object.
    void Clear()
    {
		VectorView::Clear();
		nrows = 0;
		ncols = 0;
    }

    // Resizes this matrix and optionally clears the contents.
	void Resize(index_t numrows, index_t numcols, bool fillzero = true)
    {
		VectorView::Resize(numrows*numcols, fillzero);
        nrows = numrows;
        ncols = numcols;
    }

	~Matrix() = default;
    
    // Copy operator
    Matrix & operator=(const Matrix & mat)
    {
		VectorView::operator=(mat);
		nrows = mat.nrows;
		ncols = mat.ncols;
        return *this;
    }

    // Indexing operator provides read only access to matrix element.
	// Using a call operator here for unified element access with
	// single-dimensional indices (e.g. (i)) in parent class.
	const double& operator()(index_t r, index_t c) const
    {
#ifndef NDEBUG
        if (!((static_cast<unsigned>(r) < static_cast<unsigned>(nrows)) &&
              (static_cast<unsigned>(c) < static_cast<unsigned>(ncols))))
            assert_true(0);
#endif
        return VectorView::operator()(r * ncols + c);
    }

    // Indexing operator provides read/write access to matrix element.
	// Using a call operator here for unified element access with
	// single-dimensional indices (e.g. (i)) in parent class.
	double& operator()(index_t r, index_t c)
    {
#ifndef NDEBUG
        if (!((static_cast<unsigned>(r) < static_cast<unsigned>(nrows)) &&
              (static_cast<unsigned>(c) < static_cast<unsigned>(ncols))))
            assert_true(0);
#endif
        return VectorView::operator()(r * ncols + c);
    }

    // Returns the number of rows of this matrix.
	index_t NRows() const { return nrows; }

	// Returns the number of columns of this matrix.
	index_t NCols() const { return ncols; }

	// Returns the total number of elements of this matrix.
	index_t Size() const { return nrows * ncols; }

    // Returns "true" if both matrices have the same number of rows and columns, respectively.
    bool SameSize(const Matrix & mat) const
    {
        return ((nrows == mat.nrows) && (ncols == mat.ncols));
    }

    // Returns "true" if this matrix has the same number of rows and columns as transposed(mat).
    bool SameSizeTr(const Matrix & mat) const
    {
        return ((nrows == mat.ncols) && (ncols == mat.nrows));
    }

    // Returns "true" if matrix is square.
    bool IsSquare() const { return (nrows == ncols); }

	friend std::ostream& operator<<(std::ostream& out, const Matrix& m) {
		out << "Matrix [ ";
		out << m.NRows() << "x";
		out << m.NCols() << "; ";
		out << m.data;
		out << " ]" << std::endl;
		return out;
	}

	// Serialization: Load a given Matrix
	static Matrix load(allscale::utils::ArchiveReader& reader) {
		Matrix res;
		res.nrows = reader.read<int>();
		res.ncols = reader.read<int>();
		res.data = reader.read<std::vector<double>>();
		return res;
	}

	// Serialization: Store this Matrix
	void store(allscale::utils::ArchiveWriter& writer) const {
		writer.write(nrows);
		writer.write(ncols);
		writer.write(data);
	}

};

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

