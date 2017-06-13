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
using ::std::unique_ptr;
using ::std::string;

} // end namespace app
} // end namespace amdados
