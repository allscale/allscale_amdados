//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

/////////////////////////////////////////////////////////////////////
// Error handling, messaging and logging facilities for debugging. //
/////////////////////////////////////////////////////////////////////

#if defined(AMDADOS_DEBUGGING)
#error AMDADOS_DEBUGGING macro redefinition
#endif

// A T T E N T I O N: Comment out the line below in the final release.
//#define AMDADOS_DEBUGGING



#ifdef AMDADOS_DEBUGGING

#include <cstdlib>
#include <cstdio>
#include <unistd.h>            // for access(), mkstemp()
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>
#include <mutex>

#if defined(MY_INFO) || defined(MY_TIME_IT) || \
    defined(MY_ERR) || defined(MY_TRY) || defined(MY_CATCH)
#error macro redefinition
#endif

#warning !!! AMDADOS_DEBUGGING is enabled !!!


#define MY_INFO(fmt, ...) { fprintf(stdout, fmt, __VA_ARGS__); \
                            fprintf(stdout, "\n"); \
                            fflush(stdout); }

#define MY_ERR(fmt, ...) { fprintf(stdout, "ERROR: "); \
                           fprintf(stdout, fmt, __VA_ARGS__); \
                           fprintf(stdout, "\n"); \
                           fflush(stdout); }

// Simple timer.
class M_y_T_i_m_e_r {
private:
    std::chrono::high_resolution_clock::time_point m_start_time;
public:
    explicit M_y_T_i_m_e_r(const char * text) {
        std::cout << text << std::endl << std::flush;
        m_start_time = std::chrono::high_resolution_clock::now();
    }
    ~M_y_T_i_m_e_r() {
        std::cout << "execution time: "
            << (1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - m_start_time).count())
            << " seconds" << std::endl << std::endl << std::flush;
    }
};
#define MY_TIME_IT(text) M_y_T_i_m_e_r  _t_i_m_e_r_(text); (void)_t_i_m_e_r_;

#define MY_TRY   try
#define MY_CATCH catch (const std::domain_error & e) { \
    std::cout << std::endl << "domain error: " << e.what() << std::endl; } \
                 catch (const std::runtime_error & e) { \
    std::cout << std::endl << "runtime error: " << e.what() << std::endl; } \
                 catch (const std::exception & e) { \
    std::cout << std::endl << "exception: " << e.what() << std::endl; } \
                 catch (...) { \
    std::cout << std::endl << "Unsupported exception" << std::endl; }

#else   // !AMDADOS_DEBUGGING

#define MY_INFO(x, ...)     { (void)(x); }
#define MY_TIME_IT(x)       { (void)(x); }
#define MY_TRY
#define MY_CATCH

#endif  // AMDADOS_DEBUGGING

