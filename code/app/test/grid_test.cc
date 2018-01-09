//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "allscale/api/user/data/adaptive_grid.h"

namespace amdados {

using namespace ::allscale::utils;
using ::allscale::api::user::algorithm::pfor;
using ::allscale::api::user::data::Direction;
using ::allscale::api::user::data::Direction::Up;
using ::allscale::api::user::data::Direction::Down;
using ::allscale::api::user::data::Direction::Left;
using ::allscale::api::user::data::Direction::Right;

// Number of subdomains in each dimension.
const int NUM_DOMAINS_X = 13;
const int NUM_DOMAINS_Y = 7;

// Number of elements (or nodal points) in each subdomain in each dimension.
const int SUBDOMAIN_X = 111;
const int SUBDOMAIN_Y = 17;

// Set up the configuration of a grid cell (static).
// With this type we can define a multi-resolution grid.
using sub_domain_config_t = ::allscale::api::user::data::CellConfig<2,
    ::allscale::api::user::data::layers<            //  1000m x 1000m each subdomain covers
        ::allscale::api::user::data::layer<SUBDOMAIN_X,SUBDOMAIN_Y>,// 10x10 100m nodes each consisting of
        ::allscale::api::user::data::layer<5,5>,              //  5x5   20m nodes each consisting of
        ::allscale::api::user::data::layer<5,5>               //  5x5    4m nodes
    >
>;

// Assign each layer a level corresponding to coarsest to finest.
enum {
    L_100m = 2,
    L_20m = 1,
    L_4m = 0,
};

// Position, index or size in 2D.
using point2d_t = ::allscale::api::user::data::GridPoint<2>;
using size2d_t = point2d_t;

// This is more elaborated grid of all the subdomain structures.
// These special subdomains can handle multi-resolution case.
using domain_t = ::allscale::api::user::data::Grid<
                    ::allscale::api::user::data::AdaptiveGridCell<double,sub_domain_config_t>,2>;

// Subdomain is a layer in a grid cell.
using subdomain_t = ::allscale::utils::StaticGrid<double,size_t(SUBDOMAIN_X),size_t(SUBDOMAIN_Y)>;

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

} // namespace amdados

//-------------------------------------------------------------------------------------------------
// Test the grid of subdomains.
//-------------------------------------------------------------------------------------------------
void TestGrid()
{
    std::cout << "TestGrid() ..." << std::endl;
    using namespace amdados;
    // Set to "true" for the hard test, but this can corrupt memory.
    const bool TestOutOfBounds = false;
    const bool TestYXOrder = false;

    // Origin and the number of subdomains in both directions.
    const int Ox = Origin[_X_];  const int Nx = SubDomGridSize[_X_];
    const int Oy = Origin[_Y_];  const int Ny = SubDomGridSize[_Y_];

    ::allscale::api::user::data::Grid<int,2> grid(SubDomGridSize);

    // Segmentation fault as expected.
    if (TestOutOfBounds)
    {
        for (int x = -10*SubDomGridSize[_X_]; x < 10*SubDomGridSize[_X_]; ++x) {
        for (int y = -10*SubDomGridSize[_Y_]; y < 10*SubDomGridSize[_Y_]; ++y) {
            grid[{x,y}] = x + y;
        }}
    }

    if (TestYXOrder)    // Fails - {y,x} is a wrong ordering; BEWARE: this test can overrun buffer.
    {
        int count = 0;
        for (int x = Ox; x < Nx; ++x) {
        for (int y = Oy; y < Ny; ++y) { grid[{y,x}] = count++; }}
        count = 0;
        for (int x = Ox; x < Nx; ++x) {
        for (int y = Oy; y < Ny; ++y) {
            bool ok = (grid[{y,x}] == count++);
            EXPECT_TRUE(ok);
        }}
    }
    else                // Passes - {x,y} is a correct ordering.
    {
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
    }
}

//-------------------------------------------------------------------------------------------------
// Testing the subdomain indexing.
//-------------------------------------------------------------------------------------------------
void TestSubdomainIndexing()
{
    std::cout << "TestSubdomainIndexing() ..." << std::endl;
    using namespace amdados;
    // Set to "true" for the hard test, but this can corrupt memory.
    const bool TestOutOfBounds = false;
    const bool TestYXOrder = false;

    domain_t domain1(SubDomGridSize);
    domain_t domain2(SubDomGridSize);
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        domain1[idx].setActiveLayer(L_100m);
        domain2[idx].setActiveLayer(L_100m);
        subdomain_t & subdom1 = domain1[idx].getLayer<L_100m>();
        subdomain_t & subdom2 = domain2[idx].getLayer<L_100m>();

        // Test rationale: out-of-bounds index checking. This test must cause
        // segmentation fault and in the new API version it does (but not in the old one).
        if (TestOutOfBounds)
        {
            for (int x = -100*SUBDOMAIN_X; x <= +100*SUBDOMAIN_X; ++x) {
            for (int y = -100*SUBDOMAIN_Y; y <= +100*SUBDOMAIN_Y; ++y) {
                subdom1[{x,y}] = 1.0;
            }}
        }

        // Test rationale: check the sequence of numbers is written and read correctly.
        // This passes: X first, Y second - {x,y} ordering is correct.
        {
            // ----- Loop x, then y.
            {
                domain1[idx].forAllActiveNodes([](double & value) { value = 0.0; });
                domain2[idx].forAllActiveNodes([](double & value) { value = 0.0; });

                int count = 0;
                for (int x = 0; x < SUBDOMAIN_X; ++x) {
                for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                    subdom1[{x,y}] = count;
 if (TestYXOrder) { subdom2[{y,x}] = count; }
                    ++count;
                }}
                EXPECT_TRUE(count == SUBDOMAIN_X * SUBDOMAIN_Y);

                bool ok1 = true, ok2 = true;
                count = 0;
                for (int x = 0; x < SUBDOMAIN_X; ++x) {
                for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                    if (!(subdom1[{x,y}] == count)) ok1 = false;
 if (TestYXOrder) { if (!(subdom2[{y,x}] == count)) ok2 = false; } else ok2 = false;
                    ++count;
                }}
                EXPECT_TRUE(ok1 && !ok2);
                EXPECT_TRUE(count == SUBDOMAIN_X * SUBDOMAIN_Y);
            }

            // ----- Loop y, then x.
            {
                domain1[idx].forAllActiveNodes([](double & value) { value = 0.0; });
                domain2[idx].forAllActiveNodes([](double & value) { value = 0.0; });

                int count = 0;
                for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                for (int x = 0; x < SUBDOMAIN_X; ++x) {
                    subdom1[{x,y}] = count;
 if (TestYXOrder) { subdom2[{y,x}] = count; }
                    ++count;
                }}
                EXPECT_TRUE(count == SUBDOMAIN_X * SUBDOMAIN_Y);

                bool ok1 = true, ok2 = true;
                count = 0;
                for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                for (int x = 0; x < SUBDOMAIN_X; ++x) {
                    if (!(subdom1[{x,y}] == count)) ok1 = false;
 if (TestYXOrder) { if (!(subdom2[{y,x}] == count)) ok2 = false; } else ok2 = false;
                    ++count;
                }}
                EXPECT_TRUE(ok1 && !ok2);
                EXPECT_TRUE(count == SUBDOMAIN_X * SUBDOMAIN_Y);
            }
        }
    });

    if (TestYXOrder) {
        std::cout << "Testing the wrong Y-X indexing order can corrupt memory, so we stop here"
                  << std::endl;
        std::exit(0);
    }
}

//-------------------------------------------------------------------------------------------------
// Testing the subdomain by reading/writing between subdomain boundaries.
//-------------------------------------------------------------------------------------------------
void TestInterSubdomain()
{
    std::cout << "TestInterSubdomain() ..." << std::endl;
    using namespace amdados;

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
        //for (int x = 0; x < SUBDOMAIN_X; ++x) {
        //for (int y = 0; y < SUBDOMAIN_Y; ++y) {
            //subdom[{x,y}] = GetId(idx);
        //}}
        domain[idx].forAllActiveNodes([&](double & value) { value = GetId(idx); });

        // Fill up each subdomain border by unique id, excluding (!) the corner points.
        for (int x = 1; x < SUBDOMAIN_X-1; ++x) {
            subdom[{x,0}]          = GetId(idx, Down);
            subdom[{x,SUBDOMAIN_Y-1}] = GetId(idx, Up);
        }
        for (int y = 1; y < SUBDOMAIN_Y-1; ++y) {
            subdom[{0,y}]          = GetId(idx, Left);
            subdom[{SUBDOMAIN_X-1,y}] = GetId(idx, Right);
        }
    });

    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        // For all the borders ...
        for (Direction dir : { Up, Down, Left, Right }) {
            // Skip global border.
            if ((dir == Left)  && (idx[_X_] == Ox))   continue;
            if ((dir == Right) && (idx[_X_] == Nx-1)) continue;
            if ((dir == Down)  && (idx[_Y_] == Oy))   continue;
            if ((dir == Up)    && (idx[_Y_] == Ny-1)) continue;

            // Check remote border has expected values.
            int len = ((dir == Up) || (dir == Down)) ? SUBDOMAIN_X : SUBDOMAIN_Y;
            point2d_t neighbour_idx = idx + remote_idx[dir];
            Direction neighbour_dir = remote_dir[dir];
            std::vector<double> border = domain[neighbour_idx].getBoundary(neighbour_dir);
            EXPECT_TRUE(border.size() == static_cast<size_t>(len));
            EXPECT_TRUE(border[0] == GetId(neighbour_idx));             // corner point
            for (int k = 1; k < len-1; ++k) {
                EXPECT_TRUE(border[k] == GetId(neighbour_idx, neighbour_dir));
            }
            EXPECT_TRUE(border[len-1] == GetId(neighbour_idx));         // corner point
        }
    });
}

//-------------------------------------------------------------------------------------------------
// Testing a subdomain by reading/writing to its own borders.
//-------------------------------------------------------------------------------------------------
void TestSubdomainBorders()
{
    std::cout << "TestSubdomainBorders() ..." << std::endl;
    using namespace amdados;

    domain_t domain(SubDomGridSize);
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        domain[idx].setActiveLayer(L_100m);
        subdomain_t & subdom = domain[idx].getLayer<L_100m>();

        // Fill up the subdomain by zeros.
        for (int x = 0; x < SUBDOMAIN_X; ++x) {
        for (int y = 0; y < SUBDOMAIN_Y; ++y) {
            subdom[{x,y}] = 0;
        }}

        // For all the borders ...
        for (Direction dir : { Up, Down, Left, Right }) {
            // Initialize and set border values.
            int len = ((dir == Up) || (dir == Down)) ? SUBDOMAIN_X : SUBDOMAIN_Y;
            std::vector<double> border(len);
            for (int k = 0; k < len; ++k) {
                border[k] = k + GetId(idx, dir);
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
                    for (int x = 0; x < SUBDOMAIN_X; ++x) {
                        bool ok = (subdom[{x,0}] == border[x]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
                case Up: {
                    for (int x = 0; x < SUBDOMAIN_X; ++x) {
                        bool ok = (subdom[{x,SUBDOMAIN_Y-1}] == border[x]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
                case Left: {
                    for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                        bool ok = (subdom[{0,y}] == border[y]);
                        EXPECT_TRUE(ok);
                    }
                }
                break;
                case Right: {
                    for (int y = 0; y < SUBDOMAIN_Y; ++y) {
                        bool ok = (subdom[{SUBDOMAIN_X-1,y}] == border[y]);
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
    using namespace amdados;
    std::cout << "Grid of " << NUM_DOMAINS_X << "x"
                            << NUM_DOMAINS_Y << " subdomains" << std::endl;
    std::cout << "Subdomain size " << SUBDOMAIN_X << "x" << SUBDOMAIN_Y << std::endl;

    TestGrid();
    TestSubdomainIndexing();
    TestSubdomainBorders();
    TestInterSubdomain();
}

