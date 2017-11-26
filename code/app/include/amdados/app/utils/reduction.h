//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once

namespace amdados {
namespace app {
namespace utils {

//@{
//-------------------------------------------------------------------------------------------------
// Functions for global reduction across all the subdomains.
// A T T E N T I O N: these functions must be used ONLY for testing, debugging or visualization.
//-------------------------------------------------------------------------------------------------
inline double ReduceMean(const ::allscale::api::user::data::Grid<double,2> & grid)
{
    double sum = 0.0;
    for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
    for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { sum += grid[{x,y}]; }}
    return (sum / static_cast<double>(SubDomGridSize[_X_] * SubDomGridSize[_Y_]));
}
inline double ReduceAbsMin(const ::allscale::api::user::data::Grid<double,2> & grid)
{
    double v = std::fabs(grid[{0,0}]);
    for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
    for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { v = std::min(v, std::fabs(grid[{x,y}])); }}
    return v;
}
inline double ReduceAbsMax(const ::allscale::api::user::data::Grid<double,2> & grid)
{
    double v = std::fabs(grid[{0,0}]);
    for (int x = 0; x < SubDomGridSize[_X_]; ++x) {
    for (int y = 0; y < SubDomGridSize[_Y_]; ++y) { v = std::max(v, std::fabs(grid[{x,y}])); }}
    return v;
}
//@}

} // end namespace utils
} // end namespace app
} // end namespace allscale
