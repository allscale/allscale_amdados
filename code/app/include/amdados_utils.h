//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

#include "geometry.h"

namespace amdados {

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
const char PathSep = '\\';
#else
const char PathSep = '/';
#endif

typedef ::std::vector<float>         float_array_t;
typedef ::std::vector<double>        double_array_t;
typedef ::std::vector<unsigned char> ubyte_array_t;

const double TINY = ::std::numeric_limits<double>::min() /
         ::std::pow(::std::numeric_limits<double>::epsilon(),3);

class Configuration;

/**
 * Function returns true is a value hits closed interval: v in [vmin .. vmax].
 */
template<typename T>
bool IsBounded(const T & v, const T & vmin, const T & vmax)
{
    return ((vmin <= v) && (v <= vmax));
}

/**
 * Function returns true is an integer (!) value hits semi-open interval:
 * v in [vmin .. vmax).
 */
template<typename T>
bool IsInsideRange(const T & v, const T & vmin, const T & vmax)
{
    static_assert(std::numeric_limits<T>::is_integer,
                  "integer type is expected");
    return ((vmin <= v) && (v < vmax));
}

/**
 * Function clamps the value to the interval [vmin .. vmax]
 * and returns the new value.
 */
template<typename T>
T Bound(const T & v, const T & vmin, const T & vmax)
{
    assert(vmin <= vmax);
    return std::min(std::max(v, vmin), vmax);
}

/*
 * Function rounds the value to the nearest integer.
 */
inline int Round(double val)
{
    assert_true(std::fabs(val) < std::numeric_limits<int>::max());
    return static_cast<int>(std::floor(val + 0.5));
}

void CheckFileExists(const Configuration & conf, const std::string & filename);

uint64_t RandomSeed();

std::string MakeFileName(const Configuration & conf, const std::string & what);

point2d_t GetGridSize(const Configuration & conf);

} // namespace amdados

