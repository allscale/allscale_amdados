//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {

//#################################################################################################
// P R I M A R Y  constants and types.
//#################################################################################################

// Number of subdomains in each dimension.
const int NUM_DOMAINS_X = 12;
const int NUM_DOMAINS_Y = 8;

// Number of elements (or nodal points) in each subdomain in each dimension.
const int NELEMS_X = 9;
const int NELEMS_Y = 11;

// Number of observation points (number of sensors) in X and Y dimensions in the whole domain.
// In C++ code we actually use the product OBSERVATION_NX*OBSERVATION_NY, see the constanst below
// nase NUM_SUBDOMAIN_OBSERVATIONS. However, these two constant are retained for compatibility
// with Matlab/Python code.
const int OBSERVATION_NX = 30;
const int OBSERVATION_NY = 30;

// Set up the configuration of a grid cell (static).
// With this type we can define a multi-resolution grid.
using sub_domain_config_t = ::allscale::api::user::data::CellConfig<2,
    ::allscale::api::user::data::layers<            //  1000m x 1000m each subdomain covers
        ::allscale::api::user::data::layer<NELEMS_X,NELEMS_Y>,// 10x10 100m nodes each consisting of
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

//#################################################################################################
// D E P E N D A N T  constants and types.
//#################################################################################################

// Total number of elements (or nodal points) in the entire domain in each dimension.
const int GLOBAL_NELEMS_X = NELEMS_X * NUM_DOMAINS_X;
const int GLOBAL_NELEMS_Y = NELEMS_Y * NUM_DOMAINS_Y;

// Number of elements (or nodal points) in a subdomain.
const int SUB_PROBLEM_SIZE = NELEMS_X * NELEMS_Y;

// Number of available observations in each subdomain: ceil((global #observations / #subdomains).
const int NUM_SUBDOMAIN_OBSERVATIONS =
    (OBSERVATION_NX * OBSERVATION_NY + (NUM_DOMAINS_X * NUM_DOMAINS_Y) - 1) /
                                       (NUM_DOMAINS_X * NUM_DOMAINS_Y);

// Position, index or size in 2D.
using point2d_t = ::allscale::api::user::data::GridPoint<2>;
using size2d_t = point2d_t;

// Position, index or size in 3D.
using point3d_t = ::allscale::api::user::data::GridPoint<3>;
using size3d_t = point3d_t;

// This is more elaborated grid of all the subdomain structures.
// These special subdomains can handle multi-resolution case.
using domain_t = ::allscale::api::user::data::Grid<
                    ::allscale::api::user::data::AdaptiveGridCell<double,sub_domain_config_t>,2>;

// Origin and global grid size. The latter grid is the grid of subdomains,
// where the logical coordinates give a subdomain indices in each dimension.
const point2d_t Origin = {0, 0};
const size2d_t  SubDomGridSize = {NUM_DOMAINS_X, NUM_DOMAINS_Y};

const int _X_ = 0;  // index of abscissa
const int _Y_ = 1;  // index of ordinate

//#if 1
//#define SUB2IND(x, y, SizeX, SizeY) ((x) * (SizeY) + (y))       // row major
//#else
//#define SUB2IND(x, y, SizeX, SizeY) ((x) + (SizeX) * (y))       // column major
//#endif

//-----------------------------------------------------------------------------
// Function converts 2D index to a flat 1D one.
// A T T E N T I O N: ordinate changes faster than abscissa.
//                    This is in agreement with row-major Matrix class.
//-----------------------------------------------------------------------------
inline int Sub2Ind(int x, int y)
{
#ifndef NDEBUG
    if (!((static_cast<unsigned>(x) < static_cast<unsigned>(NELEMS_X)) &&
          (static_cast<unsigned>(y) < static_cast<unsigned>(NELEMS_Y))))
        assert_true(0);
#endif
    return (x * NELEMS_Y + y);
}

//-----------------------------------------------------------------------------
// Function maps the subdomain local abscissa to global one.
//-----------------------------------------------------------------------------
inline int Sub2GloX(const point2d_t & subdomain, const int x)
{
#ifndef NDEBUG
    if (!(static_cast<unsigned>(x) < static_cast<unsigned>(NELEMS_X)))
        assert_true(0);
#endif
    return (x + subdomain[_X_] * NELEMS_X);
}
//-----------------------------------------------------------------------------
// Function maps the subdomain local ordinate to global one.
//-----------------------------------------------------------------------------
inline int Sub2GloY(const point2d_t & subdomain, const int y)
{
#ifndef NDEBUG
    if (!(static_cast<unsigned>(y) < static_cast<unsigned>(NELEMS_Y)))
        assert_true(0);
#endif
    return (y + subdomain[_Y_] * NELEMS_Y);
}

} // namespace app
} // namespace amdados

