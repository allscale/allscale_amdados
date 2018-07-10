//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {

// I M P O R T A N T: cell layout must agree with the parameters specified
// in configuration file. Namely, "subdomain_x" and "subdomain_y" must be
// equal to the number of points in the finest layer (of index 0).
//
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

// Levels of resolution of subdomain layers.
enum {
    LayerFine = 0,
    LayerLow  = 1
};

// 2D point, also "index" in parallel for (pfor) loops.
typedef ::allscale::api::user::data::GridPoint<2> point2d_t;
typedef ::std::vector<point2d_t>                  point_array_t;

// 2D size (same as 2D point).
typedef ::allscale::api::user::data::GridPoint<2> size2d_t;

// A grid cell constitutes a sub-domain.
typedef ::allscale::api::user::data::AdaptiveGridCell<double,
                                            subdomain_config_t> subdomain_t;

// Grid of sub-domains constitutes the whole domain.
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

//-----------------------------------------------------------------------------
// Function returns the grid size as a number of subdomains in both dimensions.
//-----------------------------------------------------------------------------
inline point2d_t GetGridSize(const Configuration & conf)
{
    return point2d_t(conf.asInt("num_subdomains_x"),
                     conf.asInt("num_subdomains_y"));
}

} // namespace amdados
