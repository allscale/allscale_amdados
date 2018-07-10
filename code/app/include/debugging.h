//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
//-----------------------------------------------------------------------------

#pragma once

#if defined(AMDADOS_DEBUGGING)
#error AMDADOS_DEBUGGING macro redefinition
#endif
// A T T E N T I O N: Comment out the line below in the final release.
#define AMDADOS_DEBUGGING

#if defined(MY_LOG) || defined(MY_TIME_IT)
#error macro redefinition
#endif

///////////////////////////////////////////////////////////////////////////////
// Logger.
///////////////////////////////////////////////////////////////////////////////

#if defined(AMDADOS_DEBUGGING) || defined(AMDADOS_PLAIN_MPI)

#include <iomanip>

//=============================================================================
// Simple logging implementation. It is NOT thread safe, but could easily be.
//=============================================================================
class MyLogger
{
public:
    enum LoggingType { MY_LOGGING_INFO, MY_LOGGING_WARNING, MY_LOGGING_ERROR };

    explicit MyLogger(LoggingType t) : m_type(t), m_ss() {}

    ~MyLogger() {
        extern std::fstream gLogFile;       // global log-file
        time_t t = std::time(nullptr);
        struct tm * tm_info = std::localtime(&t);   // not thread-safe
        const char * type =
            ((m_type == MY_LOGGING_INFO)    ? " [INFO] "    :
            ((m_type == MY_LOGGING_WARNING) ? " [WARNING] " : " [ERROR] "));
        gLogFile << std::put_time(tm_info, "%d-%m-%Y %H-%M-%S")
                 << type << m_ss.str() << std::endl;
#ifdef AMDADOS_PLAIN_MPI
        int rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::cout << std::put_time(tm_info, "%d-%m-%Y %H-%M-%S")
                      << type << m_ss.str() << std::endl;
        }
#endif
    }

    MyLogger & operator<<(const char * s) {
        if (s != nullptr) m_ss << s;
        return (*this);
    }

    template<typename VALUE_TYPE>
    MyLogger & operator<<(const VALUE_TYPE & v) {
        m_ss << v;
        return (*this);
    }

private:
    LoggingType       m_type;   ///< type of logging information
    std::stringstream m_ss;     ///< temporary buffer for the output
};

#else

//=============================================================================
// Stub class does nothing.
//=============================================================================
class MyLogger
{
public:
    enum LoggingType { MY_LOGGING_INFO, MY_LOGGING_WARNING, MY_LOGGING_ERROR };
    explicit MyLogger(LoggingType) {}
    ~MyLogger() {}
    MyLogger & operator<<(const char *) { return (*this); }
    template<typename VALUE_TYPE>
    MyLogger & operator<<(const VALUE_TYPE &) { return (*this); }
};

#endif  // AMDADOS_DEBUGGING || AMDADOS_PLAIN_MPI

#define MY_LOG(x) MyLogger(MyLogger::MY_LOGGING_##x)

///////////////////////////////////////////////////////////////////////////////
// Error handler.
///////////////////////////////////////////////////////////////////////////////

#if defined(AMDADOS_PLAIN_MPI)

#if defined(assert_true)
#error macro redefinition
#endif

//=============================================================================
// Simple MPI-friendly substitution to Allscale API macro "assert_true".
//=============================================================================
class MyAssertTrue
{
public:
    explicit MyAssertTrue(bool cond, const char * file, int line)
    : m_ss(), m_cond(cond) {
        if (!m_cond) {
            int rank = -1;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MY_LOG(ERROR) << "at: " << file << ":" << line
                          << ", rank: " << rank;
        }
    }

    ~MyAssertTrue() {
        void PrintErrorAndExit(const char *);               // declaration
        if (!m_cond) PrintErrorAndExit(nullptr);
    }

    MyAssertTrue & operator<<(const char * s) {
        if (!m_cond && (s != nullptr)) { m_ss << s; }
        return (*this);
    }

    template<typename VALUE_TYPE>
    MyAssertTrue & operator<<(const VALUE_TYPE & v) {
        if (!m_cond) { m_ss << v; }
        return (*this);
    }

private:
    std::stringstream m_ss;     // temporary buffer for the message
    bool              m_cond;   // status of condition
};
#define assert_true(cond) MyAssertTrue(cond, __FILE__, __LINE__)

#endif  // AMDADOS_PLAIN_MPI

///////////////////////////////////////////////////////////////////////////////
// Timer.
///////////////////////////////////////////////////////////////////////////////

#if defined(AMDADOS_PLAIN_MPI)

//=============================================================================
// Simple timer.
//=============================================================================
class MyTimer
{
private:
    double m_start_time;
public:
    explicit MyTimer(const char * text) {
        MY_LOG(INFO) << ((text == nullptr) ? "" : text);
        m_start_time = MPI_Wtime();
    }
    ~MyTimer() {
        MY_LOG(INFO) << "execution time: "
                     << (MPI_Wtime() - m_start_time) << " seconds";
    }
};
#define MY_TIME_IT(text) MyTimer  _t_i_m_e_r_(text); (void)_t_i_m_e_r_;

#elif defined(AMDADOS_DEBUGGING)

//#include <cstdlib>
//#include <cstdio>
//#include <unistd.h>            // for access(), mkstemp()
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <iomanip>
//#include <stdexcept>
#include <chrono>
//#include <mutex>

#warning !!! AMDADOS_DEBUGGING is enabled !!!

//=============================================================================
// Simple timer.
//=============================================================================
class MyTimer
{
private:
    std::chrono::high_resolution_clock::time_point m_start_time;
public:
    explicit MyTimer(const char * text) {
        MY_LOG(INFO) << text;
        m_start_time = std::chrono::high_resolution_clock::now();
    }
    ~MyTimer() {
        MY_LOG(INFO) << "execution time: "
            << (1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - m_start_time).count())
            << " seconds";
    }
};
#define MY_TIME_IT(text) MyTimer  _t_i_m_e_r_(text); (void)_t_i_m_e_r_;

#else   // !AMDADOS_DEBUGGING and !AMDADOS_PLAIN_MPI

#define MY_TIME_IT(x) { (void)(x); }

#endif  // AMDADOS_DEBUGGING || AMDADOS_PLAIN_MPI
