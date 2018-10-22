//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

void PrintErrorAndExit(const char * err_msg);

namespace amdados {

enum Directions { Up = 0, Down = 1, Left = 2, Right = 3, NSides = 4 };

//=============================================================================
// 2D point structure.
//=============================================================================
struct Point2D
{
    long x, y;

    Point2D() : x(0), y(0) {}

    Point2D(long _x, long _y) : x(_x), y(_y) {}

    bool operator==(const Point2D & p) const {
        return ((p.x == x) && (p.y == y));
    }

	bool operator<(const Point2D & p) const {
		return (x < p.x) || (x == p.x && y < p.y);
	}

    Point2D neighbour(Directions dir) const {
        switch(dir) {
            case Left:  return Point2D(x - 1, y);
            case Right: return Point2D(x + 1, y);
            case Down:  return Point2D(x, y - 1);
            case Up:    return Point2D(x, y + 1);
            default:    return Point2D(x, y);
        }
        return Point2D(x,y);
    }

    Point2D left()  const { return Point2D(x - 1, y); }
    Point2D right() const { return Point2D(x + 1, y); }
    Point2D down()  const { return Point2D(x, y - 1); }
    Point2D up()    const { return Point2D(x, y + 1); }

	friend std::ostream & operator<<(std::ostream & out, const Point2D& p) {
		return out << "P(" << p.x << "," << p.y << ")";
	}
};

typedef Point2D point2d_t;
typedef Point2D size2d_t;

typedef ::std::vector<Point2D> point_array_t;
typedef ::std::vector<int>     int_array_t;
typedef ::std::vector<float>   float_array_t;
typedef ::std::vector<double>  double_array_t;

// Flow components (flow_x, flow_y).
typedef std::pair<double,double> flow_t;

//-----------------------------------------------------------------------------
// Function returns flat index of a point with coordinates (x,y) on the grid.
// N O T E: here Y is the fastest coordinate.
//-----------------------------------------------------------------------------
inline long base_sub2ind(long x, long y, long nx, long ny)
{
    (void)nx;(void)ny;
#ifndef NDEBUG
    assert_true((0 <= x) && (x < nx) && (0 <= y) && (y < ny));
#endif
    return (x * ny + y);
}

//-----------------------------------------------------------------------------
// Function maps a subdomain local coordinates to global ones.
// \param  query_point     point in question local to subdomain.
// \param  subdomain_pos   position of the subdomain inside a grid structure.
// \param  subdomain_size  size of the subdomain.
//-----------------------------------------------------------------------------
inline point2d_t Sub2Glo(const point2d_t & query_point,
                         const point2d_t & subdomain_pos,
                         const size2d_t  & subdomain_size)
{
    const auto x = query_point.x;
    const auto y = query_point.y;
#ifndef NDEBUG
    if (!((static_cast<size_t>(x) < static_cast<size_t>(subdomain_size.x)) &&
          (static_cast<size_t>(y) < static_cast<size_t>(subdomain_size.y))))
        assert_true(0);
#endif
    return point2d_t(x + subdomain_pos.x * subdomain_size.x,
                     y + subdomain_pos.y * subdomain_size.y);
}

//-----------------------------------------------------------------------------
// Function returns point's coordinates (x,y) on the grid given flat index.
// N O T E: here Y is the fastest coordinate.
//-----------------------------------------------------------------------------
inline point2d_t base_ind2sub(long idx, long nx, long ny)
{
    (void)nx;(void)ny;
    long x = idx / ny;
    long y = idx % ny;
#ifndef NDEBUG
    assert_true(static_cast<size_t>(idx) < static_cast<size_t>(nx * ny));
    assert_true((0 <= x) && (x < nx) && (0 <= y) && (y < ny));
#endif
    return point2d_t(x,y);
}

//-----------------------------------------------------------------------------
// Function returns rank of this process.
//-----------------------------------------------------------------------------
inline int GetRank()
{
    static int rank = -1;
    if (rank >= 0) return rank;
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS) return rank;
    std::cerr << "ERROR: MPI_Comm_rank() failed" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    return -1;
}

//-----------------------------------------------------------------------------
// Function checks the return code of MPI function and aborts on error.
//-----------------------------------------------------------------------------
inline void MpiCheck(int retval, const char * file, int line)
{
    if (retval == MPI_SUCCESS) return;
    MY_LOG(ERROR) << "MPI error at " << file << ":" << line
                  << "\nreturn code: " << retval;
    PrintErrorAndExit(nullptr);
}

#ifdef MPI_CHECK
#error macro MPI_CHECK redefinition
#endif
#define MPI_CHECK(func) MpiCheck(func, __FILE__, __LINE__)

}   // namespace amdados

