//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {
namespace utils {

//-------------------------------------------------------------------------------------------------
// Flags to decide whether to keep already allocated memory or release it back to the system.
//-------------------------------------------------------------------------------------------------
enum class MemoryPolicy
{
    RETAIN_MEMORY,
    RELEASE_MEMORY
};

//-------------------------------------------------------------------------------------------------
// Function returns true is a value hits closed interval: v in [vmin .. vmax].
//-------------------------------------------------------------------------------------------------
template<typename T>
bool IsBounded(const T & v, const T & vmin, const T & vmax)
{
    return ((vmin <= v) && (v <= vmax));
}

//-------------------------------------------------------------------------------------------------
// Function returns true is an integer (!) value hits semi-open interval: v in [vmin .. vmax).
//-------------------------------------------------------------------------------------------------
template<typename T>
bool IsInsideRange(const T & v, const T & vmin, const T & vmax)
{
    static_assert(numeric_limits<T>::is_integer, "integer type is expected");
    return ((vmin <= v) && (v < vmax));
}

//-------------------------------------------------------------------------------------------------
// Function clamps the value to the interval [vmin .. vmax] and returns the new value.
//-------------------------------------------------------------------------------------------------
template<typename T>
T Bound(const T & v, const T & vmin, const T & vmax)
{
    assert(vmin <= vmax);
    return std::min(std::max(v, vmin), vmax);
}

//-------------------------------------------------------------------------------------------------
// Function returns a squared value.
//-------------------------------------------------------------------------------------------------
template<typename T>
T Square(const T & v)
{
    return (v * v);
}

//-------------------------------------------------------------------------------------------------
// Functor converts 2D index (x,y) into a plain one.
//-------------------------------------------------------------------------------------------------
template<int SizeX, int SizeY>
struct Sub2Ind {
    Sub2Ind() {
        static_assert(SizeX * SizeY < static_cast<int>(numeric_limits<int>::max()), "overflow");
    }

    int operator()(int ix, int iy) const {
        assert_true((0 <= ix) && (ix < SizeX));
        assert_true((0 <= iy) && (iy < SizeY));
        //return (x + static_cast<int>(SizeX) * y);
        return (ix * SizeY + iy);   // TODO: we already have observations based
    }                               // on this indexing: y changes faster, swap?
};

//-------------------------------------------------------------------------------------------------
// Functor converts a plane index into 2D one (x,y).
//-------------------------------------------------------------------------------------------------
template<int SizeX, int SizeY>
struct Ind2Sub {
    Ind2Sub() {
        static_assert(SizeX * SizeY < static_cast<int>(numeric_limits<int>::max()), "overflow");
    }

    void operator()(const int idx, int & x, int & y) const {
        assert_true((0 <= idx) && (idx < SizeX * SizeY));
        std::div_t divresult = std::div(idx, SizeY);    // if y is the fastest
        x = divresult.quot;
        y = divresult.rem;
    }
};

//-------------------------------------------------------------------------------------------------
// Function reshapes a vector into 2D grid structure,
// in Matlab notation: grid = reshape(vec, [SizeX, SizeY]).
//-------------------------------------------------------------------------------------------------
template<int SizeX, int SizeY, typename GRID>
void Reshape1Dto2D(GRID & grid, const allscale::utils::grid<double, SizeX * SizeY> & vec)
{
    // TODO: check grid sizes: must be SizeX by SizeY
    Sub2Ind<SizeX, SizeY> sub2ind;
    for (int i = 0; i < static_cast<int>(SizeX); i++) {
        for (int j = 0; j < static_cast<int>(SizeY); j++) {
            grid[{i,j}] = vec[{sub2ind(i,j)}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function unrolls 2D grid structure into a vector, in Matlab notation: vec = grid(:).
//-------------------------------------------------------------------------------------------------
template<int SizeX, int SizeY, typename GRID>
void Reshape2Dto1D(allscale::utils::grid<double, SizeX * SizeY> & vec, const GRID & grid)
{
    // TODO: check grid sizes: must be SizeX by SizeY
    Sub2Ind<SizeX, SizeY> sub2ind;
    for (int i = 0; i < static_cast<int>(SizeX); i++) {
        for (int j = 0; j < static_cast<int>(SizeY); j++) {
            vec[{sub2ind(i,j)}] = grid[{i,j}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function creates a new directory or does nothing if it already exists.
//-------------------------------------------------------------------------------------------------
inline void MakeDirectory(const char * dir)
{
    assert_true(dir != nullptr);
    std::string cmd("mkdir -p ");   // TODO: not portable, use STL "experimental" instead
    cmd += dir;
    int retval = std::system(cmd.c_str());  // TODO: mutex
    retval = std::system("sync");
    (void)retval;
}

} // end namespace utils
} // end namespace app
} // end namespace allscale

/*#ifndef NDEBUG*/
    /*for (int i = 0; i < MSIZE; ++i) {*/
    /*for (int j = 0; j < MSIZE; ++j) {*/
    /*std::cout << std::setw(7) << std::setprecision(4) << M[{i,j}] << " ";*/
    /*}*/
    /*std::cout << endl;*/
    /*}*/
    /*std::cout << endl;*/


//std::cout << idx[0] << "  " << idx[1] << endl << flush;
//if (idx[0]==0 && idx[1]==0)
//{
    //std::fstream f("/tmp/mat.txt", std::ios::trunc | std::ios::out);
    //for (int i = 0; i < SUB_PROBLEM_SIZE; ++i) {
    //for (int j = 0; j < SUB_PROBLEM_SIZE; ++j) {
        //f << B[{i,j}] << " ";
    //}
    //f << endl << flush;
    //}
    //std::cout << "matrix was written" << endl << flush;
//}
//std::system("sleep 100");
//assert_true(0);

/*#endif*/


