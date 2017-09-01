//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <random>
#include <memory>
#include <cmath>
#include <cstdlib>
#include <functional>

#ifndef IBM_NOINLINE
#ifdef _MSC_VER                                         // assume MS Visual C++
    #define IBM_NOINLINE __declspec(noinline)
#else                                                   // assume GCC/CLang
    #define IBM_NOINLINE __attribute__ ((noinline))
#endif
#endif // IBM_NOINLINE

namespace amdados {
namespace app {

using ::std::endl;
using ::std::flush;
using ::std::numeric_limits;
/*using ::std::unique_ptr;*/
using ::std::string;

} // end namespace app
} // end namespace amdados

// The macro can enforce inlining even in debugging mode.
#ifdef IBM_ALWAYS_INLINE
#error IBM_ALWAYS_INLINE is already defined
#endif
#if defined(__GNUC__) || defined(__GNUG__)
#define IBM_ALWAYS_INLINE __attribute__((always_inline))
#else
#define IBM_ALWAYS_INLINE
#endif

inline IBM_ALWAYS_INLINE void CheckRange1D(int x, int Nx)
{
    if (!(static_cast<unsigned>(x) < static_cast<unsigned>(Nx)))
        assert_true(0);
}

inline IBM_ALWAYS_INLINE void CheckRange2D(int x, int y, int Nx, int Ny)
{
    if (!((static_cast<unsigned>(x) < static_cast<unsigned>(Nx)) &&
          (static_cast<unsigned>(y) < static_cast<unsigned>(Ny))))
        assert_true(0);
}

