//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

#ifndef AMDADOS_PLAIN_MPI
namespace allscale {
namespace utils { class ArchiveReader; class ArchiveWriter; }}
#endif

namespace amdados {

using index_t = int64_t;

//=============================================================================
// Vector of double-precision floating-point values.
//=============================================================================
class Vector
{
protected:
	// Content of this vector, memory management done by STL.
    std::vector<double> m_data;

public:
	// Default constructor.
    Vector() : m_data() {}

	// Constructor, ensures data can fit size elements.
	Vector(index_t size) : m_data(static_cast<size_t>(size)) {}

	// Constructor.
	Vector(const Vector &) = default;

	// Move constructor.
	Vector(Vector && other) {
		m_data.swap(other.m_data);
	}

	// Destructor.
	~Vector() = default;

    // Indexing operator provides read only access to vector elements.
	// Using a call operator here for unified element access with
	// multi-dimensional indices (e.g. (i,y)) in subclasses.
	const double & operator()(index_t i) const {
#ifndef NDEBUG
        if (!(static_cast<unsigned>(i) < static_cast<unsigned>(m_data.size())))
            assert_true(0);
#endif
        return m_data[static_cast<size_t>(i)];
    }

    // Indexing operator provides read/write access to vector elements.
	// Using a call operator here for unified element access with
	// multi-dimensional indices (e.g. (i,y)) in subclasses.
	double & operator()(index_t i) {
#ifndef NDEBUG
        if (!(static_cast<unsigned>(i) < static_cast<unsigned>(m_data.size())))
            assert_true(0);
#endif
        return m_data[static_cast<size_t>(i)];
    }

    // Functions for iterating the (constant) vector.
    double * begin() { return &*m_data.begin(); }
    double * end()   { return &*m_data.end();   }

    const double * begin() const { return &*m_data.cbegin(); }
    const double * end()   const { return &*m_data.cend();   }

    // Function returns the size of this vector (view).
	index_t Size() const { return static_cast<index_t>(m_data.size()); }

    // Function returns "true" if both vectors have the same length.
	// TODO: Quite expensive operation, maybe replace with something faster.
    bool SameSize(const Vector & v) const {
        return (m_data.size() == v.m_data.size());
    }

    // Returns "true" if this vector-view refers to an empty object.
    bool Empty() const { return (m_data.empty()); }

    // Returns "true" if "v" is not the same object as "*this".
    bool IsDistinct(const Vector & v) const {
        return (&m_data != &(v.m_data));
    }

	// Deletes the content of this Vector.
	void Clear() {
		m_data.clear();
	}

	// Resizes the content to fit the given size, optionally also clearing the content.
	void Resize(index_t new_size, bool fillzero = true) {
		if (fillzero) { m_data.clear(); }
		m_data.resize(static_cast<size_t>(new_size));
	}

	// Copy operator.
	Vector & operator=(const Vector & vec) {
		m_data = vec.m_data;
		return *this;
	}

	// Prints this vector. TODO: "beautify" the output.
	friend std::ostream & operator<<(std::ostream & out, const Vector & v) {
		out << "Vector [ ";
		out << "1x" << v.Size();
		for(const auto & e : v.m_data) { out << ", " << e; }
		out << " ]" << std::endl;
		return out;
	}
};

//=============================================================================
// Matrix of double-precision floating-point values.
// Note, the matrix is row-major (the column index is faster than the row one).
// Note, implicit casting from matrix to vector is allowed.
//=============================================================================
class Matrix : public Vector
{
private:
	index_t m_nrows;    ///< number of rows
    index_t m_ncols;    ///< number of columns

public:
    // Default constructor.
    Matrix() : Vector(), m_nrows(0), m_ncols(0) { }

    // Default constructor creates matrix filled by zeros.
    Matrix(index_t numrows, index_t numcols)
    : Vector(numrows*numcols), m_nrows(numrows), m_ncols(numcols) { }

    // Copy constructor. Note, user can copy a matrix but not a vector.
    Matrix(const Matrix & mat)
    : Vector(mat), m_nrows(mat.m_nrows), m_ncols(mat.m_ncols) { }

	// Move constructor.
	Matrix(Matrix && mat)
    : Vector(std::move(mat)), m_nrows(mat.m_nrows), m_ncols(mat.m_ncols)
	{
		mat.m_nrows = 0;
		mat.m_ncols = 0;
	}

	// Destructor.
	~Matrix() = default;

	// Deallocates and clears this object.
    void Clear()
    {
		Vector::Clear();
		m_nrows = 0;
		m_ncols = 0;
    }

    // Resizes this matrix and optionally clears the contents.
	void Resize(index_t numrows, index_t numcols, bool fillzero = true)
    {
		Vector::Resize(numrows * numcols, fillzero);
        m_nrows = numrows;
        m_ncols = numcols;
    }

    // Copy operator.
    Matrix & operator=(const Matrix & mat)
    {
		Vector::operator=(mat);
		m_nrows = mat.m_nrows;
		m_ncols = mat.m_ncols;
        return *this;
    }

    // Indexing operator provides read only access to matrix element.
	// Using a call operator here for unified element access with
	// single-dimensional indices (e.g. (i)) in parent class.
	const double & operator()(index_t r, index_t c) const
    {
#ifndef NDEBUG
        if (!((static_cast<unsigned>(r) < static_cast<unsigned>(m_nrows)) &&
              (static_cast<unsigned>(c) < static_cast<unsigned>(m_ncols))))
            assert_true(0);
#endif
        return Vector::operator()(r * m_ncols + c);
    }

    // Indexing operator provides read/write access to matrix element.
	// Using a call operator here for unified element access with
	// single-dimensional indices (e.g. (i)) in parent class.
	double & operator()(index_t r, index_t c)
    {
#ifndef NDEBUG
        if (!((static_cast<unsigned>(r) < static_cast<unsigned>(m_nrows)) &&
              (static_cast<unsigned>(c) < static_cast<unsigned>(m_ncols))))
            assert_true(0);
#endif
        return Vector::operator()(r * m_ncols + c);
    }

    // Returns the number of rows of this matrix.
	index_t NRows() const { return m_nrows; }

	// Returns the number of columns of this matrix.
	index_t NCols() const { return m_ncols; }

	// Returns the total number of elements of this matrix.
	index_t Size() const { return m_nrows * m_ncols; }

    // Swaps content of two matrices.
    void swap(Matrix & x) {
        this->m_data.swap(x.m_data);
        std::swap(this->m_nrows, x.m_nrows);
        std::swap(this->m_ncols, x.m_ncols);
    }

    // Returns "true" for matrices with the same dimensions.
    bool SameSize(const Matrix & mat) const {
        return ((m_nrows == mat.m_nrows) && (m_ncols == mat.m_ncols));
    }

    // Returns "true" if this matrix has the same number of rows and columns
    // as transposed(mat).
    bool SameSizeTr(const Matrix & mat) const {
        return ((m_nrows == mat.m_ncols) && (m_ncols == mat.m_nrows));
    }

    // Returns "true" if matrix is square.
    bool IsSquare() const { return (m_nrows == m_ncols); }

    // Prints this matrix. TODO: "beautify" the output.
	friend std::ostream & operator<<(std::ostream & out, const Matrix & m) {
		out << "Matrix [ ";
		out << m.NRows() << "x";
		out << m.NCols() << "; " << std::endl;
        for (index_t r = 0; r < m.NRows(); ++r) {
            for (index_t c = 0; c < m.NCols(); ++c) { out << m(r, c) << "  "; }
            out << std::endl;
        }
		out << " ]" << std::endl;
		return out;
	}

#ifndef AMDADOS_PLAIN_MPI
	// Serialization: load this matrix.
	static Matrix load(::allscale::utils::ArchiveReader & reader);
	// Serialization: store this matrix.
	void store(::allscale::utils::ArchiveWriter & writer) const;
#endif
};

void MatMult(Matrix & result, const Matrix & A, const Matrix & B);

void MatMultTr(Matrix & result, const Matrix & A, const Matrix & B);

void MatVecMult(Vector & result, const Matrix & A, const Vector & v);

void AddVectors(Vector & result, const Vector & a, const Vector & b);

void SubtractVectors(Vector & result, const Vector & a, const Vector & b);

void AddMatrices(Matrix & result, const Matrix & A, const Matrix & B);

void SubtractMatrices(Matrix & result, const Matrix & A, const Matrix & B);

void Fill(Vector & v, double vfill = 0.0);

void MakeIdentityMatrix(Matrix & A);

void GetTransposed(Matrix & At, const Matrix & A);

void Symmetrize(Matrix & A);

void ScalarMult(Vector & v, const double mult);

double Norm(const Vector & v);

double NormDiff(const Vector & a, const Vector & b);

double Trace(const Matrix & A);

void Negate(Vector & v);

void MakeRandom(Vector & v, const char type);

bool CheckNoNan(const Vector & v);

void WriteVector(const Vector & vec,
                        const std::string & filename, int precision = 4);

void WriteMatrix(const Matrix & A,
                        const std::string & filename, int precision = 4);

} // namespace amdados

