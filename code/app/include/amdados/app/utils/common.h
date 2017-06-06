#pragma once
//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <limits>
#include <random>
#include <memory>
#include <cassert>
#include <cmath>

#ifndef IBM_NOINLINE
#ifdef _MSC_VER                                         // assume MS Visual C++
    #define IBM_NOINLINE __declspec(noinline)
#else                                                   // assume GCC/CLang
    #define IBM_NOINLINE __attribute__ ((noinline))
#endif
#endif // IBM_NOINLINE

