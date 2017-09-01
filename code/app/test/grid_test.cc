//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "allscale/api/user/operator/pfor.h"
#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"

namespace amdados {
namespace app {

using namespace ::allscale::utils;
using ::allscale::api::user::pfor;

// Number of subdomains in each dimension.
const int NUM_DOMAINS_X = 12;
const int NUM_DOMAINS_Y = 8;

// Number of elements (or nodal points) in each subdomain in each dimension.
const int NELEMS_X = 9;
const int NELEMS_Y = 11;

// Set up the configuration of a grid cell (static).
// With this type we can define a multi-resolution grid.
using sub_domain_config_t = CellConfig<
    layers<                             //  1000m x 1000m each subdomain covers
        layer<NELEMS_X,NELEMS_Y>,       //  10 x 10  100m nodes each consisting of
        layer<5,5>,                     //   5 x  5   20m nodes each consisting of
        layer<5,5>                      //   5 x  5    4m nodes
    >
>;

// Assign each layer a level corresponding to coarsest to finest.
enum {
    L_100m = 2,
    L_20m = 1,
    L_4m = 0,
};

// Total number of elements (or nodal points) in the entire domain in each dimension.
const int GLOBAL_NELEMS_X = NELEMS_X * NUM_DOMAINS_X;
const int GLOBAL_NELEMS_Y = NELEMS_Y * NUM_DOMAINS_Y;

// Position, index or size in 2D.
using point2d_t = ::allscale::api::user::data::GridPoint<2>;
using size2d_t = point2d_t;

// This is more elaborated grid of all the subdomain structures.
// These special subdomains can handle multi-resolution case.
using domain_t = ::allscale::api::user::data::Grid< Cell<double,sub_domain_config_t>, 2 >;

// Subdomain is a layer in a grid cell.
using subdomain_t = ::allscale::utils::grid<double,NELEMS_X,NELEMS_Y>;

// Origin and global grid size. The latter grid is the grid of subdomains,
// where the logical coordinates give a subdomain indices in each dimension.
const point2d_t Origin = {0, 0};
const size2d_t  SubDomGridSize = {NUM_DOMAINS_X, NUM_DOMAINS_Y};

const int _X_ = 0;  // index of abscissa
const int _Y_ = 1;  // index of ordinate

// Unique ID from subdomain location.
int GetId(point2d_t idx)
{
    bool ok1 = ((0 <= idx[_X_]) && (idx[_X_] < NUM_DOMAINS_X));
    bool ok2 = ((0 <= idx[_Y_]) && (idx[_Y_] < NUM_DOMAINS_Y));
    assert_true(ok1 && ok2);
    return (idx[_X_] * NUM_DOMAINS_Y + idx[_Y_]);
}

// Unique ID from subdomain location and its boundary side index.
int GetId(point2d_t idx, Direction dir)
{
    return (GetId(idx) + (static_cast<int>(dir) + 1) * NUM_DOMAINS_X * NUM_DOMAINS_Y);
}

} // namespace app
} // namespace amdados

//-------------------------------------------------------------------------------------------------
// Test the grid of subdomains.
//-------------------------------------------------------------------------------------------------
void TestGrid()
{
    using namespace amdados::app;

    // Origin and the number of subdomains in both directions.
    const int Ox = Origin[_X_];  const int Nx = SubDomGridSize[_X_];
    const int Oy = Origin[_Y_];  const int Ny = SubDomGridSize[_Y_];

    ::allscale::api::user::data::Grid<int,2> grid(SubDomGridSize);

#if 0
    // Segmentation fault as expected.
    for (int x = -10*SubDomGridSize[_X_]; x < 10*SubDomGridSize[_X_]; ++x) {
    for (int y = -10*SubDomGridSize[_Y_]; y < 10*SubDomGridSize[_Y_]; ++y) {
        grid[{x,y}] = x + y;
    }}
#endif

#if 1
    // Passes - {x,y} is a correct ordering.
    // Loop x, then y.
    int count = 0;
    for (int x = Ox; x < Nx; ++x) {
    for (int y = Oy; y < Ny; ++y) {
        grid[{x,y}] = count++;
    }}
    count = 0;
    for (int x = Ox; x < Nx; ++x) {
    for (int y = Oy; y < Ny; ++y) {
        bool ok = (grid[{x,y}] == count++);
        EXPECT_TRUE(ok);
    }}

    // Loop y, then x.
    count = 0;
    for (int y = Oy; y < Ny; ++y) {
    for (int x = Ox; x < Nx; ++x) {
        grid[{x,y}] = count++;
    }}
    count = 0;
    for (int y = Oy; y < Ny; ++y) {
    for (int x = Ox; x < Nx; ++x) {
        bool ok = (grid[{x,y}] == count++);
        EXPECT_TRUE(ok);
    }}
#else
    // Fails - {y,x} is a wrong ordering; this test can overrun buffer, beware.
    int count = 0;
    for (int x = Ox; x < Nx; ++x) {
    for (int y = Oy; y < Ny; ++y) { grid[{y,x}] = count++; }}
    count = 0;
    for (int x = Ox; x < Nx; ++x) {
    for (int y = Oy; y < Ny; ++y) {
        bool ok = (grid[{y,x}] == count++);
        EXPECT_TRUE(ok);
    }}
#endif
}

//-------------------------------------------------------------------------------------------------
// Testing the subdomain indexing.
//-------------------------------------------------------------------------------------------------
void TestSubdomainIndexing()
{
    using namespace amdados::app;

    domain_t domain(SubDomGridSize);
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        domain[idx].setActiveLayer(L_100m);
        subdomain_t & subdom = domain[idx].getLayer<L_100m>();

        // Test rationale: out-of-bounds index checking. This test must cause
        // segmentation fault but it does not because modular arithmetic prevents
        // error checking, impossible to track down the out-of-bounds errors.
        {
            for (int x = -10*NELEMS_X; x < 10*NELEMS_X; ++x) {
            for (int y = -10*NELEMS_Y; y < 10*NELEMS_Y; ++y) {
                subdom[{x,y}] = 1.0;
            }}
        }

        // Test rationale: check the sequence of numbers is written and read correctly.
        // This passes: X first, Y second - {x,y} ordering is correct.
        {
            // ----- Loop x, then y.
            int count = 0;
            for (int x = 0; x < NELEMS_X; ++x) {
            for (int y = 0; y < NELEMS_Y; ++y) {
                subdom[{x,y}] = count++;
            }}
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);
            count = 0;
            for (int x = 0; x < NELEMS_X; ++x) {
            for (int y = 0; y < NELEMS_Y; ++y) {
                bool ok = (subdom[{x,y}] == count++);
                EXPECT_TRUE(ok);
            }}
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);

            // ----- Loop y, then x.
            count = 0;
            for (int y = 0; y < NELEMS_Y; ++y) {
            for (int x = 0; x < NELEMS_X; ++x) {
                subdom[{x,y}] = count++;
            }}
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);
            count = 0;
            for (int y = 0; y < NELEMS_Y; ++y) {
            for (int x = 0; x < NELEMS_X; ++x) {
                bool ok = (subdom[{x,y}] == count++);
                EXPECT_TRUE(ok);
            }}
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);
        }

        // Test rationale: check the sequence of numbers is written and read correctly.
        // This fails: Y first, X second - {y,x} ordering is wrong.
        {
            // ----- Loop x, then y.
            int count = 0;
            for (int x = 0; x < NELEMS_X; ++x) {
            for (int y = 0; y < NELEMS_Y; ++y) {
                subdom[{y,x}] = count++;
            }}
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);
            count = 0;
            bool ok = true;
            for (int x = 0; x < NELEMS_X; ++x) {
            for (int y = 0; y < NELEMS_Y; ++y) {
                if (subdom[{y,x}] != count) ok = false;
                ++count;
            }}
            EXPECT_FALSE(ok);
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);

            // ----- Loop y, then x.
            count = 0;
            for (int y = 0; y < NELEMS_Y; ++y) {
            for (int x = 0; x < NELEMS_X; ++x) {
                subdom[{y,x}] = count++;
            }}
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);
            count = 0;
            ok = true;
            for (int y = 0; y < NELEMS_Y; ++y) {
            for (int x = 0; x < NELEMS_X; ++x) {
                if (subdom[{y,x}] != count) ok = false;
                ++count;
            }}
            EXPECT_FALSE(ok);
            EXPECT_TRUE(count == NELEMS_X * NELEMS_Y);
        }
    });
}

//-------------------------------------------------------------------------------------------------
// Testing the subdomain by reading/writing between subdomain boundaries.
//-------------------------------------------------------------------------------------------------
void TestInterSubdomain()
{
    using namespace amdados::app;

    // Origin and the number of subdomains in both directions.
    const int Ox = Origin[_X_];  const int Nx = SubDomGridSize[_X_];
    const int Oy = Origin[_Y_];  const int Ny = SubDomGridSize[_Y_];

    // Index increments and corresponding directions of a neighbour (remote) subdomain.
    Direction remote_dir[4];
    point2d_t remote_idx[4];

    remote_dir[Left ] = Right;  remote_idx[Left ] = point2d_t{-1,0};
    remote_dir[Right] = Left;   remote_idx[Right] = point2d_t{+1,0};
    remote_dir[Down ] = Up;     remote_idx[Down ] = point2d_t{0,-1};
    remote_dir[Up   ] = Down;   remote_idx[Up   ] = point2d_t{0,+1};

    EXPECT_TRUE(static_cast<int>(std::min(std::min(Left,Right), std::min(Down,Up))) == 0);
    EXPECT_TRUE(static_cast<int>(std::max(std::max(Left,Right), std::max(Down,Up))) == 3);

    domain_t domain(SubDomGridSize);
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        domain[idx].setActiveLayer(L_100m);
        subdomain_t & subdom = domain[idx].getLayer<L_100m>();

        // Fill up the subdomain by unique id.
        for (int x = 0; x < NELEMS_X; ++x) {
        for (int y = 0; y < NELEMS_Y; ++y) {
            subdom[{x,y}] = GetId(idx);
        }}

        // Fill up each subdomain border by unique id.
        for (int x = 1; x < NELEMS_X-1; ++x) {
            subdom[{x,0}]          = GetId(idx, Down);
            subdom[{x,NELEMS_Y-1}] = GetId(idx, Up);
        }
        for (int y = 1; y < NELEMS_Y-1; ++y) {
            subdom[{0,y}]          = GetId(idx, Left);
            subdom[{NELEMS_X-1,y}] = GetId(idx, Right);
        }
    });

    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        // For all the borders ...
        for (Direction dir : { Up, Down, Left, Right }) {
            // Skip global border.
            if ((dir == Left)  && (idx[0] == Ox))   continue;
            if ((dir == Right) && (idx[0] == Nx-1)) continue;
            if ((dir == Down)  && (idx[1] == Oy))   continue;
            if ((dir == Up)    && (idx[1] == Ny-1)) continue;

            // Check remote border has expected values.
            int len = ((dir == Up) || (dir == Down)) ? NELEMS_X : NELEMS_Y;
            point2d_t neighbour_idx = idx + remote_idx[dir];
            Direction neighbour_dir = remote_dir[dir];
            std::vector<double> border = domain[neighbour_idx].getBoundary(neighbour_dir);
            EXPECT_TRUE(border.size() == static_cast<size_t>(len));
            EXPECT_TRUE(border[0] == GetId(neighbour_idx));
            for (int k = 1; k < len-1; ++k) {
                EXPECT_TRUE(border[k] == GetId(neighbour_idx, neighbour_dir));
            }
            EXPECT_TRUE(border[len-1] == GetId(neighbour_idx));
        }
    });
}

//-------------------------------------------------------------------------------------------------
// Testing a subdomain by reading/writing to its own borders.
//-------------------------------------------------------------------------------------------------
void TestSubdomainBorders()
{
    using namespace amdados::app;

    domain_t domain(SubDomGridSize);
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        domain[idx].setActiveLayer(L_100m);
        subdomain_t & subdom = domain[idx].getLayer<L_100m>();

        // Fill up the subdomain by zeros.
        for (int x = 0; x < NELEMS_X; ++x) {
        for (int y = 0; y < NELEMS_Y; ++y) {
            subdom[{x,y}] = 0;
        }}

        // For all the borders ...
        for (Direction dir : { Up, Down, Left, Right }) {
            // Initialize and set border values.
            int len = ((dir == Up) || (dir == Down)) ? NELEMS_X : NELEMS_Y;
            std::vector<double> border(len);
            for (int k = 0; k < len; ++k) {
                border[k] = k + static_cast<int>(31*(dir + 1));
            }
            domain[idx].setBoundary(dir, border);

            // Read back and compare for equality.
            std::vector<double> replica = domain[idx].getBoundary(dir);
            EXPECT_TRUE(border.size() == replica.size());
            bool ok = std::equal(border.begin(), border.end(), replica.begin());
            EXPECT_TRUE(ok);

            // Ckeck that border values were properly written by reading the subdomain directly.
            switch (dir) {
                case Down: {
                    for (int x = 0; x < NELEMS_X; ++x) {
                        bool ok = (subdom[{x,0}] == border[x]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
                case Up: {
                    for (int x = 0; x < NELEMS_X; ++x) {
                        bool ok = (subdom[{x,NELEMS_Y-1}] == border[x]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
                case Left: {
                    for (int y = 0; y < NELEMS_Y; ++y) {
                        bool ok = (subdom[{0,y}] == border[y]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
                case Right: {
                    for (int y = 0; y < NELEMS_Y; ++y) {
                        bool ok = (subdom[{NELEMS_X-1,y}] == border[y]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
            }
        }
    });
}

//-------------------------------------------------------------------------------------------------
// Function tests grid indexing and subdomain layout.
//-------------------------------------------------------------------------------------------------
TEST(GridTest, Basic)
{
    TestGrid();
    TestSubdomainIndexing();
    TestSubdomainBorders();
    TestInterSubdomain();
}

