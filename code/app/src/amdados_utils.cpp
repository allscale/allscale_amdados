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
 * Function clears the output directory removing all the files generated
 * by the previous call to this application.
 * N O T E: this function is designed for debugging.
 *          In production mode it does nothing.
 */
//void CleanOutputDir(const std::string & dir, const char * ext)
//{
//    (void) dir; (void) ext;
//#ifdef AMDADOS_DEBUGGING
//    int retval = 0;
//    assert_true(!dir.empty());
//    std::string cmd = std::string("/bin/rm -fr ") + dir + PathSep;
//    retval = std::system("sync");
//    if (ext == nullptr) {
//        retval = std::system((cmd + "field*.png").c_str());
//        retval = std::system((cmd + "field*.pgm").c_str());
//        retval = std::system((cmd + "field*.jpg").c_str());
//        retval = std::system((cmd + "field*.avi").c_str());
//        retval = std::system((cmd + "*.log").c_str());
//        retval = std::system((cmd + "*.out").c_str());
//        retval = std::system((cmd + "*profile*.txt").c_str());
//    } else {
//        assert_true((std::strlen(ext) == 4) && (ext[0] == '.'))
//                << "wrong extension; expected format: '.abc'";
//        retval = std::system(((cmd + "*") + ext).c_str());
//    }
//    retval = std::system("sync");
//    (void) retval;
//#endif
//}

/**
 * Function creates a video file from a sequence of field states written into
 * image files using 'ffmpeg' utility installed system-wide. By the end of
 * video creation, all *.pgm files are removed in order to save space.
 */
void MakeVideo(const Configuration & conf, const std::string & file_title)
{
    (void) conf; (void) file_title;
#ifdef AMDADOS_DEBUGGING
    MY_INFO("%s", "\n\n")
    const int framerate = 24;
    int retval = 0;
    std::stringstream cmd, wildcards;
    wildcards << conf.asString("output_dir") << PathSep
              << file_title << "*.pgm";
    // Generate AVI file, if 'ffmpeg' was installed.
    cmd << "ffmpeg -y -f image2 -framerate " << framerate
        << " -pattern_type glob -i '" << wildcards.str() << "' "
        << MakeFileName(conf, file_title, ".avi");
    retval = std::system(cmd.str().c_str());
    retval = std::system("sync");
    // Remove the individual field images as they were encoded into AVI file.
    cmd.str(std::string());
    cmd << "/bin/rm -fr " << wildcards.str();
    retval = std::system(cmd.str().c_str());
    retval = std::system("sync");
    (void) retval;
    MY_INFO("%s", "\n\n")
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
        assert_true(entity != "field")
                << "file suffix is always expected for the field entity";
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
