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

// Maximum number of nodal points in a subdomain in each dimension.
//const int SUBDOMAIN_X = 16;
//const int SUBDOMAIN_Y = 16;

// Set up the configuration of a grid cell (static).
// With this type we can define a multi-resolution grid.
#if 0
// layer: 0, size: [16,16]
// layer: 1, size: [ 8, 8]
// layer: 2, size: [ 4, 4]
// layer: 3, size: [ 2, 2]
// layer: 4, size: [ 1, 1]
typedef ::allscale::api::user::data::CellConfig<2,
            ::allscale::api::user::data::layers<
                    ::allscale::api::user::data::layer<2,2>,
                    ::allscale::api::user::data::layer<2,2>,
                    ::allscale::api::user::data::layer<2,2>,
                    ::allscale::api::user::data::layer<2,2>
            >
> subdomain_config_t;

// Assign each layer a level corresponding to coarsest to finest.
enum {
    LayerFine = 0,
    LayerMid  = 1,
    LayerLow  = 2
};
#else
// layer: 0, size: [16,16]
// layer: 1, size: [ 8, 8]
// layer: 2, size: [ 1, 1]
typedef ::allscale::api::user::data::CellConfig<2,
            ::allscale::api::user::data::layers<
                    ::allscale::api::user::data::layer<1,1>,
                    ::allscale::api::user::data::layer<8,8>,
                    ::allscale::api::user::data::layer<2,2>
            >
> subdomain_config_t;

// Assign each layer a level corresponding to coarsest to finest.
enum {
    LayerFine = 0,
    LayerLow  = 1
};
#endif

//const int _X_ = 0;      // index of abscissa
//const int _Y_ = 1;      // index of ordinate

//#############################################################################
// D E R I V E D  constants and types.
//#############################################################################

// Number of elements (or nodal points) in a subdomain.
//const int SUB_PROBLEM_SIZE = SUBDOMAIN_X * SUBDOMAIN_Y;

// 2D point, also "index" in parallel for (pfor) loops.
typedef ::allscale::api::user::data::GridPoint<2> point2d_t;
typedef ::std::vector<point2d_t>                  point_array_t;

// 2D size (same as 2D point).
typedef ::allscale::api::user::data::GridPoint<2> size2d_t;

// A grid cell that constitutes a sub-domain.
typedef ::allscale::api::user::data::AdaptiveGridCell<double,
                                            subdomain_config_t> subdomain_t;

// Collection of sub-domains constitutes the whole domain.
typedef ::allscale::api::user::data::Grid<subdomain_t, 2> domain_t;

/**
 * Function maps a subdomain local coordinates to global ones.
 * \param  query_point     point in question local to subdomain.
 * \param  subdomain_pos   position of the subdomain inside a grid structure.
 * \param  subdomain_size  size of the subdomain.
 */
inline point2d_t Sub2Glo(const point2d_t & query_point,
                         const point2d_t & subdomain_pos,
                         const size2d_t  & subdomain_size)
{
    const long x = query_point.x;
    const long y = query_point.y;
#ifndef NDEBUG
    if (!((static_cast<size_t>(x) < static_cast<size_t>(subdomain_size.x)) &&
          (static_cast<size_t>(y) < static_cast<size_t>(subdomain_size.y))))
        assert_true(0);
#endif
    return point2d_t(x + subdomain_pos.x * subdomain_size.x,
                     y + subdomain_pos.y * subdomain_size.y);
}

} // namespace amdados
