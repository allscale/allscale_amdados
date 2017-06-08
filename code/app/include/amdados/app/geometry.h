//-----------------------------------------------------------------------------
// Author    : Fearghal O'Donncha, feardonn@ie.ibm.com
//             Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {

// Number of sub-domains in each dimension.
const int NUM_DOMAINS_X = 10;
const int NUM_DOMAINS_Y = 10;

// Number of elements (or nodal points) in each sub-domain in each dimension.
const int NELEMS_X = 10;
const int NELEMS_Y = 10;

// Set up the configuration of a grid cell (static).
using sub_domain_config_t = CellConfig<
    layers<                             //  1000m x 1000m each sub-domain covers
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

} // namespace app
} // namespace amdados

