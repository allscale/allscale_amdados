//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <map>
#include "allscale/utils/assert.h"
#include "allscale/api/user/data/adaptive_grid.h"
#include "../include/debugging.h"
#include "../include/geometry.h"
#include "../include/amdados_utils.h"
#include "../include/configuration.h"

namespace amdados {

/**
 * Function returns "true" if the string "str" ends with "ending".
 */
bool EndsWith(const std::string & str, const std::string & ending)
{
    const size_t ns = str.size();
    const size_t ne = ending.size();
    if (ns < ne) return false;
    return (str.compare(ns - ne, ne, ending) == 0);
}

/**
 * Function checks if the file exists and prints error message if it does not.
 */
void CheckFileExists(const Configuration & conf, const std::string & filename)
{
    (void) conf; (void) filename;
#ifdef AMDADOS_DEBUGGING
    if (!std::fstream(filename, std::ios::in).good()) {
        MY_ERR("failed to open file (not existing?): %s", filename.c_str());
        std::exit(1);
    }
#endif
}

/**
 * Function queries the high-resolution clock until it has changed. The last
 * obtained value is then used as a good seed guaranteed to be unique.
 */
uint64_t RandomSeed()
{
    using hrclock_t = std::chrono::high_resolution_clock;
    uint64_t prev = hrclock_t::now().time_since_epoch().count();
    uint64_t seed = 0;
    while ((seed = hrclock_t::now().time_since_epoch().count()) == prev) {}
    return seed;
}

/**
 * Function returns the file name for sensor locations, analytic solution,
 * simulation solution or field given configuration settings. Important,
 * grid resolution and the number of time steps (except for sensor locations)
 * are encrypted into the file name. The latter makes the results obtained
 * with different settings to be distinguishable.
 */
std::string MakeFileName(const Configuration & conf,
                         const std::string & entity, const char * suffix)
{
    if (!((entity == "sensors")  ||
          (entity == "analytic") ||
          (entity == "solution") ||
          (entity == "field"))) {
        MY_INFO("%s", "allowed entity strings: "
                "{'sensors', 'analytic', 'simulation', 'field'}");
        assert_true(0) << "wrong entity to make a file name for";
    }

    int Nx = conf.asInt("num_subdomains_x") * conf.asInt("subdomain_x");
    int Ny = conf.asInt("num_subdomains_y") * conf.asInt("subdomain_y");

    std::stringstream filename;
    filename << conf.asString("output_dir") << PathSep << entity
             << "_Nx" << Nx << "_Ny" << Ny;
    if (entity != "sensors") {
        filename << "_Nt" << conf.asInt("Nt");
    }

    if (suffix == nullptr) {
//        assert_true(entity != "field")
//                << "file suffix is always expected for the field entity";
        filename << ".txt";
    } else {
        if (EndsWith(suffix, ".avi")) {
            filename << suffix;
        } else if (EndsWith(suffix, ".png") || EndsWith(suffix, ".txt")) {
            filename << "_" << suffix;
        } else {
            assert_true(0) << "file suffix does not have expected extension";
        }
    }

    const char * fmt = (entity == "field") ? "\n%s%s" : "%s%s";
    MY_INFO(fmt, "File name: ", filename.str().c_str()); (void)fmt;
    return filename.str();
}

/**
 * Function returns the grid size as a number of subdomains
 * in both dimensions.
 */
point2d_t GetGridSize(const Configuration & conf)
{
    return point2d_t(conf.asInt("num_subdomains_x"),
                     conf.asInt("num_subdomains_y"));
}

} // namespace amdados
