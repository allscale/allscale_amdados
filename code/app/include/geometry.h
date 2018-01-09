//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

//#############################################################################
// P R I M A R Y  constants and types.
//#############################################################################

// Number of nodal points in a subdomain in each dimension.
const int SUBDOMAIN_X = 11;
const int SUBDOMAIN_Y = 11;

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

//#############################################################################
// D E P E N D A N T  constants and types.
//#############################################################################

// Number of elements (or nodal points) in a subdomain.
const int SUB_PROBLEM_SIZE = SUBDOMAIN_X * SUBDOMAIN_Y;

// Position, index or size in 2D.
typedef ::allscale::api::user::data::GridPoint<2> point2d_t;
typedef ::allscale::api::user::data::GridPoint<2> size2d_t;
typedef ::std::vector<point2d_t>                  point_array_t;

// This is more elaborated grid of all the subdomain structures.
// These special subdomains can handle multi-resolution case.
using domain_t = ::allscale::api::user::data::Grid<
                    ::allscale::api::user::data::AdaptiveGridCell<
                                            double,sub_domain_config_t>, 2>;

const int _X_ = 0;      // index of abscissa
const int _Y_ = 1;      // index of ordinate

/**
 * Function converts 2D index to a flat 1D one.
 * A T T E N T I O N: ordinate changes faster than abscissa.
 *                    This is in agreement with row-major Matrix class.
 */
inline int Sub2Ind(int x, int y)
{
#ifndef NDEBUG
    if (!((static_cast<unsigned>(x) < static_cast<unsigned>(SUBDOMAIN_X)) &&
          (static_cast<unsigned>(y) < static_cast<unsigned>(SUBDOMAIN_Y))))
        assert_true(0);
#endif
    return (x * SUBDOMAIN_Y + y);
}

/**
 * Function return the global sizes of a field.
 */
template<typename T>
point2d_t GlobalSize(const ::allscale::api::user::data::Grid<T,2> & field)
{
    return point2d_t(field.size()[_X_] * SUBDOMAIN_X,
                     field.size()[_Y_] * SUBDOMAIN_Y);
}

/**
 * Function maps the subdomain local abscissa to global one.
 */
inline int Sub2GloX(const point2d_t & subdomain, const int x)
{
#ifndef NDEBUG
    if (!(static_cast<unsigned>(x) < static_cast<unsigned>(SUBDOMAIN_X)))
        assert_true(0);
#endif
    return (x + subdomain[_X_] * SUBDOMAIN_X);
}

/**
 * Function maps the subdomain local ordinate to global one.
 */
inline int Sub2GloY(const point2d_t & subdomain, const int y)
{
#ifndef NDEBUG
    if (!(static_cast<unsigned>(y) < static_cast<unsigned>(SUBDOMAIN_Y)))
        assert_true(0);
#endif
    return (y + subdomain[_Y_] * SUBDOMAIN_Y);
}

/**
 * Function converts global coordinates to the subdomain index on the grid.
 */
inline point2d_t Glo2CellIndex(const point2d_t & p)
{
    return point2d_t(p.x / SUBDOMAIN_X, p.y / SUBDOMAIN_Y);
}

/**
 * Function converts global coordinates to the local ones inside a subdomain.
 */
inline point2d_t Glo2Sub(const point2d_t & p)
{
    return point2d_t(p.x % SUBDOMAIN_X, p.y % SUBDOMAIN_Y);
}

} // namespace amdados

