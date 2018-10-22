//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
//             Fearghal O'Donncha, feardonn@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

class SubDomain;

//=============================================================================
// Global grid of subdomains, where individual subdomain can be instantiated
// only in a single process. If a subdomain does not belong to this process,
// it is represented by NULL pointer in the corresponding grid cell.
//=============================================================================
class MpiGrid
{
private:
    typedef ::std::vector<SubDomain*> subdom_array_t;

    subdom_array_t m_subdoms;   ///< pointers to all the subdomains
    int_array_t    m_sd_ranks;  ///< subdomain ranks
    int            m_nprocs;    ///< number of concurrent processes
    int            m_rank;      ///< rank of this process
    int64_t        m_first;     ///< index of first subdomain of this process
    int64_t        m_last;      ///< index of last subdomain of this process
    size2d_t       m_size;      ///< grid size in number of subdomains

public:
//-----------------------------------------------------------------------------
// Constructor initializes grid structure and a sub-set of subdomains
// belonging to this process.
//-----------------------------------------------------------------------------
MpiGrid(const Configuration & conf)
    : m_subdoms()
    , m_nprocs(0)
    , m_rank(0)
    , m_first(0)
    , m_last(0)
    , m_size(0,0)
{
    bool ok1 = (MPI_Comm_size(MPI_COMM_WORLD, &m_nprocs) == MPI_SUCCESS);
    bool ok2 = (MPI_Comm_rank(MPI_COMM_WORLD, &m_rank) == MPI_SUCCESS);
    assert_true(ok1 && ok2) << "failed to query MPI engine";
    m_size = size2d_t(conf.asInt("num_subdomains_x"),
                      conf.asInt("num_subdomains_y"));
    int64_t N = m_size.x * m_size.y;

    // Although we use long integers, in many places it is assumed that
    // "int" should be sufficient to accommodate sizes.
    assert_true(N < std::numeric_limits<int>::max()) <<
        "too many subdomains: 'int' range was exceeded";

    // Get MPI_TAG_UB attribute which limits the maximum value stored in tag.
    long tag_ub = 0;
    {
        int flag = 0;
        void * pval = nullptr;
        MPI_CHECK(MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &pval, &flag));
        if (!flag) {
            assert_true(0) << "could not get MPI_TAG_UB";
        } else {
            tag_ub = *(reinterpret_cast<int*>(pval));
            assert_true(!(tag_ub < 32767))
                << "too-small value (%d) for MPI_TAG_UB: " << tag_ub;
            MY_LOG(INFO) << "MPI_TAG_UB: " << tag_ub;
        }
    }

    // We use tag for identification of subdomain and its boundaries.
    // This value could be potentially overflowed on very large problems.
    assert_true(N * (int)NSides < tag_ub) <<
        "too many subdomains: the tag value in MPI message can be overflowed"
        "\nMPI_TAG_UB = " << tag_ub << "\n\n";

    // Lambda function converts rank to flatten grid position.
    auto Rank2Pos = [=](int r) -> int64_t { return ((r * N) / m_nprocs); };

    // Flat-index interval occupied by subdomains attached to this process.
    m_first = Rank2Pos(m_rank);
    m_last  = Rank2Pos(m_rank + 1);

    // Create array of all the subdomains. Those attached to other processes
    // will be NULL. Those attached to this process will be set up later.
    m_subdoms.resize((size_t)N);
    std::fill(m_subdoms.begin(), m_subdoms.end(), nullptr);

    // Compute ranks of all the subdomains including those on other processes.
    m_sd_ranks.resize((size_t)N);
    for (int r = 0; r < m_nprocs; ++r) {
        std::fill(m_sd_ranks.begin() + Rank2Pos(r),
                  m_sd_ranks.begin() + Rank2Pos(r + 1), r);
    }
}

//-----------------------------------------------------------------------------
// Destructor releases the sub-set of subdomains of this process.
//-----------------------------------------------------------------------------
virtual ~MpiGrid()
{
}

//-----------------------------------------------------------------------------
// Returns pointer to a subdomain at the cell (x,y); can be NULL.
//-----------------------------------------------------------------------------
SubDomain * getSubdomain(const point2d_t & pos) const
{
    return m_subdoms[(size_t)sub2ind(pos)];
}

//-----------------------------------------------------------------------------
// Returns pointer to a subdomain at flat index position; can be NULL.
//-----------------------------------------------------------------------------
SubDomain * getSubdomain(long flat_index) const
{
    assert_true(static_cast<size_t>(flat_index) < m_subdoms.size());
    return m_subdoms[static_cast<size_t>(flat_index)];
}

//-----------------------------------------------------------------------------
// Returns pointer to a subdomain at the cell (x,y); can be NULL.
//-----------------------------------------------------------------------------
void setSubdomain(const point2d_t & pos, SubDomain * sd)
{
    size_t i = (size_t)sub2ind(pos);
    assert_true(m_subdoms[i] == nullptr) << "subdomain must be set just once";
    m_subdoms[i] = sd;
}

//-----------------------------------------------------------------------------
// Returns process rank of a subdomain at the cell (x,y).
//-----------------------------------------------------------------------------
int getRank(const point2d_t pos) const
{
    return m_sd_ranks[(size_t)sub2ind(pos)];
}

//-----------------------------------------------------------------------------
// Returns rank of this process.
//-----------------------------------------------------------------------------
int myRank() const
{
    return m_rank;
}

//-----------------------------------------------------------------------------
// Returns grid size in number of subdomains.
//-----------------------------------------------------------------------------
size2d_t getGridSize() const
{
    return m_size;
}

//-----------------------------------------------------------------------------
// Returns flat index of a cell with coordinates (x,y) on the grid.
//-----------------------------------------------------------------------------
long sub2ind(long x, long y) const
{
    return base_sub2ind(x, y, m_size.x, m_size.y);
}

//-----------------------------------------------------------------------------
// Returns flat index of a cell with coordinates (p.x,p.y) on the grid.
//-----------------------------------------------------------------------------
long sub2ind(const point2d_t & p) const
{
    return base_sub2ind(p.x, p.y, m_size.x, m_size.y);
}

//-----------------------------------------------------------------------------
// Returns cell coordinates (x,y) on the grid given flat index.
//-----------------------------------------------------------------------------
point2d_t ind2sub(long idx) const
{
    return base_ind2sub(idx, m_size.x, m_size.y);
}

//-----------------------------------------------------------------------------
// Returns "true" for sub-domain located inside the whole domain.
//-----------------------------------------------------------------------------
bool is_inside(const point2d_t & subdom_pos) const
{
    return ((0 <= subdom_pos.x) && (subdom_pos.x < m_size.x) &&
            (0 <= subdom_pos.y) && (subdom_pos.y < m_size.y));
}

//-----------------------------------------------------------------------------
// Function loops over all grid cells (subdomains) and invokes lambda function
// of the following signature: void op(point2d_t subdom_position,
//                                     int       subdom_rank);
// Note, both local and remote subdomains will be listed.
//-----------------------------------------------------------------------------
template<typename OPERATION>
void forAll(const OPERATION & op)
{
    point2d_t pos;
    for (pos.x = 0; pos.x < m_size.x; pos.x++) {
    for (pos.y = 0; pos.y < m_size.y; pos.y++) {
        op(pos, getRank(pos));
    }}
}

//-----------------------------------------------------------------------------
// Function loops over all local (!) grid cells (subdomains) and invokes lambda
// function with the following signature: void op(SubDomain * sd);
//-----------------------------------------------------------------------------
template<typename OPERATION>
void forAllLocal(const OPERATION & op)
{
    for (int64_t k = m_first; k < m_last; ++k) {
        op(m_subdoms[(size_t)k]);
    }
}

//-----------------------------------------------------------------------------
// Function loops over all local (!) grid cells (subdomains) and invokes lambda
// function of the following signature: void op(const SubDomain * sd);
// Mind "const" modifier of this function.
//-----------------------------------------------------------------------------
template<typename OPERATION>
void forAllLocal(const OPERATION & op) const
{
    for (int64_t k = m_first; k < m_last; ++k) {
        op(m_subdoms[(size_t)k]);
    }
}

};  // class MpiGrid

}   // namespace amdados







////-----------------------------------------------------------------------------
//// Returns "iterator" to the first subdomain of this process.
////-----------------------------------------------------------------------------
//SubDomain ** begin()
//{
//    return (m_subdoms.data() + m_first);
//}
//
////-----------------------------------------------------------------------------
//// Returns "iterator" to the next after last subdomain of this process.
////-----------------------------------------------------------------------------
//SubDomain ** end()
//{
//    return (m_subdoms.data() + m_last);
//}
