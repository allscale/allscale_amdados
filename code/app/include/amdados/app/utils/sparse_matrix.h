//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2014
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {
namespace utils {

//=================================================================================================
// Triplet is used to initialise entries of a sparse matrix.
//=================================================================================================
class Triplet
{
public:
    typedef double value_t;
    typedef long   index_t;

    // Constructor.
    inline Triplet() : mRow(0), mCol(0), mVal(0) {}

    // Constructor.
    inline Triplet(index_t r, index_t c, value_t v) : mRow(r), mCol(c), mVal(v) {}

    // Setter.
    inline void Set(index_t r, index_t c, value_t v) {
        mRow = r;
        mCol = c;
        mVal = v;
    }

    // Operator is designed to support triplet ordering (sorting).
    inline bool operator<(const Triplet & x) const {
        return ((mRow < x.mRow) || ((mRow == x.mRow) && (mCol < x.mCol)));
    }

    // Function returns the row index of the element.
    inline index_t row() const { return mRow; }

    // Function returns the column index of the element.
    inline index_t col() const { return mCol; }

    // Function returns the value of the element.
    inline value_t value() const { return mVal; }

private:
    index_t mRow;       ///< row index
    index_t mCol;       ///< column index
    value_t mVal;       ///< value
};

//=================================================================================================
/** Sparse matrix with very minimum functionality. */
//=================================================================================================
template<int NROWS, int NCOLS>
class SpMatrix
{
public:
    using value_t = Triplet::value_t;
    using index_t = Triplet::index_t;
    using triplet_t = Triplet;

    SpMatrix();
    void SetFromTriplets(std::vector<triplet_t> & triplets, MemoryPolicy mempol);

    long NRows() const { return static_cast<long>(NROWS); }     ///< number of rows
    long NCols() const { return static_cast<long>(NCOLS); }     ///< number of columns

    template<int NROWS_other, int NCOLS_other>
        friend void SparseMulVector(VectorView< NROWS_other >                  & result,
                                    const SpMatrix< NROWS_other, NCOLS_other > & A,
                                    const VectorView< NCOLS_other >            & vec);

    template<int MSIZE, int NROWS_other, int NCOLS_other>
        friend void SparseMulDense(Matrix< NROWS_other, NCOLS_other >   & result,
                                   const SpMatrix< NROWS_other, MSIZE > & A,
                                   const Matrix< MSIZE, NCOLS_other >   & B);

    template<int MSIZE, int NROWS_other, int NCOLS_other>
        friend void DenseMulSparse(Matrix< NROWS_other, NCOLS_other >   & result,
                                   const Matrix< NROWS_other, MSIZE >   & A,
                                   const SpMatrix< MSIZE, NCOLS_other > & B);

    template<int MSIZE, int NROWS_other, int NCOLS_other>
        friend void DenseMulSparseTr(Matrix< NROWS_other, NCOLS_other >   & result,
                                     const Matrix< NROWS_other, MSIZE >   & A,
                                     const SpMatrix< NCOLS_other, MSIZE > & B);

    template<int NROWS_other, int NCOLS_other>
        friend void MatVecMult(VectorView< NROWS_other >                  & result,
                               const SpMatrix< NROWS_other, NCOLS_other > & A,
                               const VectorView< NCOLS_other >            & vec);

    template<int MSIZE, int NROWS_other, int NCOLS_other>
        friend void MatMult(Matrix< NROWS_other, NCOLS_other >   & result,
                            const SpMatrix< NROWS_other, MSIZE > & A,
                            const Matrix< MSIZE, NCOLS_other >   & B);

    template<int MSIZE, int NROWS_other, int NCOLS_other>
        friend void MatMult(Matrix< NROWS_other, NCOLS_other >   & result,
                            const Matrix< NROWS_other, MSIZE >   & A,
                            const SpMatrix< MSIZE, NCOLS_other > & B);

    template<int MSIZE, int NROWS_other, int NCOLS_other>
        friend void MatMultTr(Matrix< NROWS_other, NCOLS_other >   & result,
                              const Matrix< NROWS_other, MSIZE >   & A,
                              const SpMatrix< NCOLS_other, MSIZE > & B);

private:
    // Column index and value of matrix element.
    struct ColumnAndValue {
        value_t val;                            ///< value of matrix element
        index_t col;                            ///< column index
        ColumnAndValue() : val(0), col(-1) {}   ///< default constructor
    };

    std::vector<ColumnAndValue>  mColVals;      ///< column indexes and entry values
    std::vector<ColumnAndValue*> mRows;         ///< pointers to matrix rows

private:
    SpMatrix(const SpMatrix &);                 // non-copyable
    SpMatrix & operator=(const SpMatrix &);     // non-copyable
};

//-------------------------------------------------------------------------------------------------
// Constructor.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
SpMatrix<NROWS,NCOLS>::SpMatrix() : mColVals(), mRows()
{
    static_assert(sizeof(value_t) == sizeof(index_t), "type sizes mismatch");
}

//-------------------------------------------------------------------------------------------------
// Function initializes this matrix from an array of triplets.
// A T T E N T I O N: the input array "triplets" will be sorted upon completion.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void SpMatrix<NROWS,NCOLS>::SetFromTriplets(std::vector<triplet_t> & triplets, MemoryPolicy mempol)
{
    assert((0 < NROWS) && (0 < NCOLS) && !triplets.empty());

    // Sort triplets in row-major order.
    std::sort(triplets.begin(), triplets.end());

    // Count the number of elements with unique row-column index pairs and check the bounds.
    index_t N_unique = 1;
    {
        auto curr = triplets.begin();
        assert(IsInsideRange<index_t>(curr->row(), 0, static_cast<index_t>(NROWS)));
        assert(IsInsideRange<index_t>(curr->col(), 0, static_cast<index_t>(NCOLS)));
        for (auto prev = curr++; curr != triplets.end(); prev = curr++) {
            if ((prev->row() != curr->row()) ||
                (prev->col() != curr->col())) {
                ++N_unique;
            }
            assert(IsInsideRange<index_t>(curr->row(), 0, static_cast<index_t>(NROWS)));
            assert(IsInsideRange<index_t>(curr->col(), 0, static_cast<index_t>(NCOLS)));
        }
    }

    // Allocate array for holding the sparse matrix infrastructure.
    mColVals.clear();
    mRows.clear();
    if (mempol == MemoryPolicy::RELEASE_MEMORY) {       // deallocate completely
        mColVals.shrink_to_fit();
        mRows.shrink_to_fit();
    }
    mColVals.resize(N_unique);
    mRows.resize(NROWS + 1);

    // Build the sparse matrix summing up the duplicate entries.
    auto v = mColVals.begin();
    index_t rowNo = 0;
    auto it = triplets.begin();
    while (rowNo < static_cast<index_t>(NROWS)) {
        mRows[rowNo] = &(*v);
        while ((it != triplets.end()) && (rowNo == it->row())) {
            index_t colNo = it->col();
            double sum = it->value();
            while ((++it != triplets.end()) && (rowNo == it->row()) && (colNo == it->col())) {
                sum += it->value();
            }
            v->val = static_cast<value_t>(sum);
            v->col = colNo;
            ++v;
        }
        ++rowNo;
    }
    assert(rowNo == static_cast<index_t>(NROWS));
    assert((v == mColVals.end()) && (it == triplets.end()));
    mRows[NROWS] = &(*v);
}

//-------------------------------------------------------------------------------------------------
// Sparse times dense vector multiplication: dense_result = sparse_A * dense_vector.
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
void SparseMulVector(VectorView<NROWS> & result,
                     const SpMatrix<NROWS,NCOLS> & A, const VectorView<NCOLS> & vec)
{
    assert_true(CheckDistinctObjects(result, vec));
    for (int r = 0; r < NROWS; ++r) {
        double sum = 0.0;
        for (auto it = A.mRows[r]; it != A.mRows[r + 1]; ++it) {
            sum += it->val * vec(static_cast<int>(it->col));
        }
        result(r) = sum;
    }
}
template<int NROWS, int NCOLS>
void MatVecMult(VectorView<NROWS> & result,
        const SpMatrix<NROWS,NCOLS> & A, const VectorView<NCOLS> & vec)
{
    SparseMulVector(result, A, vec);
}

//-------------------------------------------------------------------------------------------------
// Sparse times dense matrix multiplication: dense_result = sparse_A * dense_B.
//-------------------------------------------------------------------------------------------------
template<int MSIZE, int NROWS, int NCOLS>
void SparseMulDense(Matrix<NROWS,NCOLS> & result,
                    const SpMatrix<NROWS,MSIZE> & A, const Matrix<MSIZE,NCOLS> & B)
{
    assert_true(CheckDistinctObjects(result, B));
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) {
        double sum = 0.0;
        for (auto it = A.mRows[r]; it != A.mRows[r + 1]; ++it) {
            sum += it->val * B(static_cast<int>(it->col),c);
        }
        result(r,c) = sum;
    }}
}
template<int MSIZE, int NROWS, int NCOLS>
void MatMult(Matrix<NROWS,NCOLS> & result,
        const SpMatrix<NROWS,MSIZE> & A, const Matrix<MSIZE,NCOLS> & B)
{
    SparseMulDense(result, A, B);
}

//-------------------------------------------------------------------------------------------------
// Dense times sparse matrix multiplication: dense_result = dense_A * sparse_B.
//-------------------------------------------------------------------------------------------------
template<int MSIZE, int NROWS, int NCOLS>
void DenseMulSparse(Matrix<NROWS,NCOLS> & result,
                    const Matrix<NROWS,MSIZE> & A, const SpMatrix<MSIZE,NCOLS> & B)
{
    assert_true(CheckDistinctObjects(result, A));
    FillMatrix(result, 0.0);
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < MSIZE; ++c) {
        for (auto it = B.mRows[c]; it != B.mRows[c + 1]; ++it) {
            result(r,static_cast<int>(it->col)) += A(r,c) * it->val;
        }
    }}
}
template<int MSIZE, int NROWS, int NCOLS>
void MatMult(Matrix<NROWS,NCOLS> & result,
        const Matrix<NROWS,MSIZE> & A, const SpMatrix<MSIZE,NCOLS> & B)
{
    DenseMulSparse(result, A, B);
}

//-------------------------------------------------------------------------------------------------
// Dense times sparse^T matrix multiplication: dense_result = dense_A * transposed(sparse_B).
// Note, the matrix B is not explicitly transposed, rather traversed accordingly.
//-------------------------------------------------------------------------------------------------
template<int MSIZE, int NROWS, int NCOLS>
void DenseMulSparseTr(Matrix<NROWS,NCOLS> & result,
                      const Matrix<NROWS,MSIZE> & A, const SpMatrix<NCOLS,MSIZE> & B)
{
    assert_true(CheckDistinctObjects(result, A));
    for (int r = 0; r < NROWS; ++r) {
    for (int c = 0; c < NCOLS; ++c) {
        double sum = 0.0;
        for (auto it = B.mRows[c]; it != B.mRows[c + 1]; ++it) {
            sum += A(r,static_cast<int>(it->col)) * it->val;
        }
        result(r,c) = sum;
    }}
}
template<int MSIZE, int NROWS, int NCOLS>
void MatMultTr(Matrix<NROWS,NCOLS> & result,
               const Matrix<NROWS,MSIZE> & A, const SpMatrix<NCOLS,MSIZE> & B)
{
    DenseMulSparseTr(result, A, B);
}

//-------------------------------------------------------------------------------------------------
// Function generates a random sparse matrix given density.
// Triplets are returned for testing (to compare against Armadillo, for example).
//-------------------------------------------------------------------------------------------------
template<int NROWS, int NCOLS>
std::vector<Triplet> MakeRandomSpMatrix(SpMatrix<NROWS,NCOLS> & A, double density)
{
    using index_t = Triplet::index_t;
    using value_t = Triplet::value_t;

    const index_t NR = static_cast<index_t>(NROWS);
    const index_t NC = static_cast<index_t>(NCOLS);

    density = Bound(density, 0.01, 0.3);    // validate
    const index_t N = static_cast<index_t>(std::ceil(NR * NC * density));

    std::vector<Triplet> triplets;
    triplets.reserve(N);

    std::mt19937                            gen(std::time(nullptr));
    std::uniform_int_distribution<index_t>  row_distrib(0, NR - 1);
    std::uniform_int_distribution<index_t>  col_distrib(0, NC - 1);
    std::uniform_real_distribution<value_t> val_distrib(1.0, 2.0);
    for (index_t k = 0; k < N; ++k) {
        triplets.push_back(Triplet(row_distrib(gen), col_distrib(gen), val_distrib(gen)));
    }
    A.SetFromTriplets(triplets, MemoryPolicy::RELEASE_MEMORY);
    return triplets;
}

} // end namespace utils
} // end namespace app
} // end namespace amdados



/*// Function returns true if both triplets address the same sparse matrix entry.*/
/*static bool IsSameEntry(const Triplet & a, const Triplet & b) {*/
/*return ((a.mRow == b.mRow) && (a.mCol == b.mCol));*/
/*}*/

/*//-------------------------------------------------------------------------------------------------*/
/**//** Function returns the diagonal of this matrix as a vector. All diagonal elements must exist. */
/*//-------------------------------------------------------------------------------------------------*/
/*template<int NROWS, int NCOLS, typename VALUE, typename INDEX>*/
/*template<typename VECTOR>*/
/*void SpMatrix<VALUE,INDEX>::GetDiagonal(VECTOR & diag) const*/
/*{*/
/*assert((0 < NROWS) && (NROWS == NCOLS));*/
/*if (index_t(diag.size()) != NROWS) diag.resize(NROWS);*/
/*index_t count = 0;*/
/*for (index_t rowNo = 0; rowNo < NROWS; ++rowNo) {*/
/*const ColumnAndValue * r = mRows[rowNo];*/
/*const ColumnAndValue * e = mRows[rowNo+1];*/
/*for (; r < e; ++r) {*/
/*if (r->col == rowNo) {*/
/*diag[rowNo] = r->val;*/
/*++count;*/
/*break;*/
/*}*/
/*}*/
/*}*/
/*assert(count == NROWS);*/
/*}*/

