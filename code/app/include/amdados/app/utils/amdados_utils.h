#pragma once
#include "amdados/app/utils/matrix.h"

using namespace allscale::api::user;

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
    static_assert(std::numeric_limits<T>::is_integer, "integer type is expected");
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

double mu1(double timestep){
    return  -0.6 * sin(timestep/10 - M_PI) * 0.2;
}

double mu2(double timestep){
    return  -1.2 * sin(timestep/5 - M_PI) * 0.2;
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

//-------------------------------------------------------------------------------------------------
// Functor converts 2D index (x,y) into a plain one.
//-------------------------------------------------------------------------------------------------
template<size_t SizeX, size_t SizeY>
struct Sub2Ind {
    Sub2Ind() {
        static_assert(SizeX * SizeY < static_cast<size_t>(std::numeric_limits<int>::max()),
                      "overflow");
    }

    int operator()(int x, int y) const {
        assert(static_cast<unsigned int>(x) < SizeX);
        assert(static_cast<unsigned int>(y) < SizeY);
        //return (x + static_cast<int>(SizeX) * y);
        return (x * static_cast<int>(SizeY) + y);   // TODO: we already have observations based
    }                                               // on this indexing: y changes faster, swap?
};

//-------------------------------------------------------------------------------------------------
// Function reshapes a vector into 2D grid structure,
// in Matlab notation: grid = reshape(vec, [NELEMS_X, NELEMS_Y]).
//-------------------------------------------------------------------------------------------------
template<size_t NELEMS_X, size_t NELEMS_Y, typename GRID>
void Reshape1Dto2D(GRID & grid, const Vector<NELEMS_X * NELEMS_Y> & vec)
{
    // TODO: check grid sizes: must be NELEMS_X by NELEMS_Y
    Sub2Ind<NELEMS_X, NELEMS_Y> sub2ind;
    for (int i = 0; i < static_cast<int>(NELEMS_X); i++) {
        for (int j = 0; j < static_cast<int>(NELEMS_Y); j++) {
            grid[{i,j}] = vec[{sub2ind(i,j)}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// Function unrolls 2D grid structure into a vector, in Matlab notation: vec = grid(:).
//-------------------------------------------------------------------------------------------------
template<size_t NELEMS_X, size_t NELEMS_Y, typename GRID>
void Reshape2Dto1D(Vector<NELEMS_X * NELEMS_Y> & vec, const GRID & grid)
{
    // TODO: check grid sizes: must be NELEMS_X by NELEMS_Y
    Sub2Ind<NELEMS_X, NELEMS_Y> sub2ind;
    for (int i = 0; i < static_cast<int>(NELEMS_X); i++) {
        for (int j = 0; j < static_cast<int>(NELEMS_Y); j++) {
            vec[{sub2ind(i,j)}] = grid[{i,j}];
        }
    }
}

//-------------------------------------------------------------------------------------------------
// As temporary measure make explicit loop that sets all values to zero.
//-------------------------------------------------------------------------------------------------
template<size_t SizeX, size_t SizeY>
void SetAllValsToZero(allscale::utils::grid<double,SizeX,SizeY>& init)
{
    for(int i=0; i<(int)SizeX; ++i) {      //SizeX and SizeY are number of elements in X and Y
        for(int j=0; j<(int)SizeY; ++j) {      //SizeX and SizeY are number of elements in X and
            init[{i,j}] = 1000.;   // Set all vals to zero
        }
    }

}

//-------------------------------------------------------------------------------------------------
// Function creates a new directory or does nothing if it already exists.
//-------------------------------------------------------------------------------------------------
void MakeDirectory(const char * dir)
{
    assert_true(dir != nullptr);
    std::string cmd("mkdir -p ");   // TODO: not portable
    cmd += dir;
    int retval = std::system(cmd.c_str());  // TODO: mutex
    retval = std::system("sync");
    (void)retval;
}

} // end namespace utils
} // end namespace app
} // end namespace allscale
