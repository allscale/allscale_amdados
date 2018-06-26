//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017-2018
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
 * Function checks if the file exists and prints error message if it does not.
 */
#ifdef AMDADOS_DEBUGGING
void CheckFileExists(const Configuration & conf, const std::string & filename)
{
    if (!std::fstream(filename, std::ios::in).good()) {
        MY_ERR("failed to open file (not existing?): %s", filename.c_str());
        std::exit(1);
    }
}
#else
void CheckFileExists(const Configuration &, const std::string &) {}
#endif

/**
 * Function queries the high-resolution clock until it has changed. The last
 * obtained value is then used as a good seed guaranteed to be unique.
 */
uint64_t RandomSeed()
{
    //using hrclock_t = std::chrono::high_resolution_clock;
    //uint64_t prev = hrclock_t::now().time_since_epoch().count();
    //uint64_t seed = 0;
    //while ((seed = hrclock_t::now().time_since_epoch().count()) == prev) {}
    //return seed;
	return 2063;
}

/**
 * Function returns the file name for: (1) sensor locations ("sensors"), (2)
 * analytic solution ("analytic") or (3) state field ("field") given
 * configuration settings. Important, grid resolution and the number of time
 * steps (except for sensor locations) are encrypted into the file name. This
 * helps to distinguish simulations with different settings.
 */
std::string MakeFileName(const Configuration & conf, const std::string & what)
{
    int Nx = conf.asInt("num_subdomains_x") * conf.asInt("subdomain_x");
    int Ny = conf.asInt("num_subdomains_y") * conf.asInt("subdomain_y");

    std::stringstream filename;
    filename << conf.asString("output_dir") << PathSep << what
             << "_Nx" << Nx << "_Ny" << Ny;

    if (what == "sensors") {
        filename << ".txt";
    } else if (what == "analytic") {
        filename << "_Nt" << conf.asInt("Nt") << ".txt";
    } else if (what == "field") {
        filename << "_Nt" << conf.asInt("Nt") << ".bin";
    } else {
        assert_true(0) << "unknown entity to make a file name from";
    }

    MY_INFO(((what == "field") ? "\n%s%s" : "%s%s"),
            "File name: ", filename.str().c_str());

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
