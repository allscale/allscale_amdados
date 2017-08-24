//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include "allscale/api/user/operator/pfor.h"
#include "allscale/utils/assert.h"
#include "amdados/app/amdados_grid.h"

namespace amdados {
namespace app {

using namespace ::allscale::utils;
using namespace ::amdados::app;
using ::allscale::api::user::pfor;
using ::allscale::api::user::data::Grid;
using ::allscale::api::user::data::GridPoint;

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

// Origin and global grid size. The latter grid is the grid of subdomains,
// where the logical coordinates give a subdomain indices in each dimension.
const point2d_t Origin = {0, 0};
const size2d_t  SubDomGridSize = {NUM_DOMAINS_X, NUM_DOMAINS_Y};

const int _X_ = 0;  // index of abscissa
const int _Y_ = 1;  // index of ordinate

} // namespace app
} // namespace amdados

int FailTest1()
{
    using namespace ::amdados::app;

    {
        ::allscale::api::user::data::Grid<int,2> grid(SubDomGridSize);

#if 0
        // Segmentation fault as expected.
        for (int x = -10*SubDomGridSize[_X_]; x < 10*SubDomGridSize[_X_]; ++x) {
        for (int y = -10*SubDomGridSize[_Y_]; y < 10*SubDomGridSize[_Y_]; ++y) {
            grid[{x,y}] = x + y;
        }}
#endif

#if 1
        // Passes - {x,y} is correct.
        int count = 0;
        for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
        for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { grid[{x,y}] = count++; }}
        count = 0;
        for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
        for (int y = 0; y < SubDomGridSize[_Y_]; ++y) {
            bool ok = (grid[{x,y}] == count++);
            assert_true(ok);
        }}
#else
        // Fails - {y,x} is wrong.
        int count = 0;
        for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
        for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { grid[{y,x}] = count++; }}
        count = 0;
        for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
        for (int y = 0; y < SubDomGridSize[_Y_]; ++y) {
            bool ok = (grid[{y,x}] == count++);
            assert_true(ok);
        }}
#endif
    }

    domain_t glo_tmp(SubDomGridSize);
    pfor(Origin, SubDomGridSize, [&](const point2d_t & idx) {
        glo_tmp[idx].setActiveLayer(L_100m);
        for (Direction dir : { Up, Down, Left, Right }) {
            auto & tmp = glo_tmp[idx];
            ::allscale::utils::grid<double,NELEMS_X,NELEMS_Y> & subdom = tmp.getLayer<L_100m>();

            {
                // Test rationale: out-of-bounds index checking.
                // This fails: modular arithmetic prevents error checking, impossible to
                // track down the out-of-bounds errors.
                for (int x = -10*NELEMS_X; x < 10*NELEMS_X; ++x) {
                for (int y = -10*NELEMS_Y; y < 10*NELEMS_Y; ++y) {
                    subdom[{x,y}] = double(dir + 10);
                }}
            }

            {
                // Test rationale: check the sequence of numbers is written and read correctly.
                // This passes: X first, Y second - {x,y} ordering is correct.
                int count = 0;
                for (int x = 0; x < NELEMS_X; ++x) {
                for (int y = 0; y < NELEMS_Y; ++y) {
                    subdom[{x,y}] = double(count);
                    ++count;
                }}
                assert_true(count == NELEMS_X * NELEMS_Y);
                count = 0;
                for (int x = 0; x < NELEMS_X; ++x) {
                for (int y = 0; y < NELEMS_Y; ++y) {
                    bool ok = (subdom[{x,y}] == double(count));
                    assert_true(ok);
                    ++count;
                }}
                assert_true(count == NELEMS_X * NELEMS_Y);
            }

            {
                // Test rationale: check the sequence of numbers is written and read correctly.
                // This fails: Y first, X second - {y,x} ordering is wrong.
                int count = 0;
                for (int x = 0; x < NELEMS_X; ++x) {
                for (int y = 0; y < NELEMS_Y; ++y) {
                    subdom[{y,x}] = double(count);
                    ++count;
                }}
                assert_true(count == NELEMS_X * NELEMS_Y);
                count = 0;
                for (int x = 0; x < NELEMS_X; ++x) {
                for (int y = 0; y < NELEMS_Y; ++y) {
                    bool ok = (subdom[{y,x}] == double(count));
                    assert_true(ok);
                    ++count;
                }}
                assert_true(count == NELEMS_X * NELEMS_Y);
            }

/*
            auto local_boundary = tmp.getBoundary(dir);
            std::fill(local_boundary.begin(), local_boundary.end(), double(dir));
            tmp.setBoundary(dir, local_boundary);
            assert_true(std::equal(local_boundary.begin(), local_boundary.end(),
                                   tmp.getBoundary(dir).begin()));
*/
        }
    });

    return 1;
 }

