//-----------------------------------------------------------------------------
// Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
// Copyright : IBM Research Ireland, 2017
//-----------------------------------------------------------------------------

#include <cstdlib>
#include <string>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include "allscale/utils/assert.h"

namespace amdados {
namespace app {
namespace utils {

//-----------------------------------------------------------------------------
// Function creates a new directory or does nothing if it already exists.
// TODO: this will not work on Windows, use STL "experimental" instead.
//-----------------------------------------------------------------------------
void MakeDirectory(const char * dir)
{
    assert_true(dir != nullptr);
    std::string cmd("mkdir -p ");
    cmd += dir;
    int retval = std::system(cmd.c_str());
    retval = std::system("sync");
    (void)retval;
}

//-----------------------------------------------------------------------------
// Function creates an output directory inside the current one, which is
// supposed to be the project root folder. If the directory is already exist
// all its content will be deleted.
// TODO: this will not work on Windows, use STL "experimental" instead.
//-----------------------------------------------------------------------------
void CreateAndCleanOutputDir(const std::string & dir)
{
    assert_true(!dir.empty());
    std::string cmd("mkdir -p ");
    cmd += dir;
    int retval = std::system(cmd.c_str());
    retval = std::system("sync");
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*.png").c_str());
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*.pgm").c_str());
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*.jpg").c_str());
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*.avi").c_str());
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*.log").c_str());
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*.out").c_str());
    retval = std::system(
            (std::string("/bin/rm -fr ") + dir + "/*profile*.txt").c_str());
    retval = std::system("sync");
    (void) retval;
}

//-----------------------------------------------------------------------------
// Function opens a text file in a temporary directory; useful for unit tests.
//-----------------------------------------------------------------------------
void OpenTextFileForUnitTest(std::fstream & file, const char * basename)
{
    file.flush();
    file.close();
    assert_true(basename != nullptr);
    std::string filename = basename;
#if defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#warning "Still need to decide what is the best practice for Windows"
#else
    int retval = std::system("(mkdir -p /tmp/amdados) && sync");
    (void) retval;
    filename.insert(0, "/tmp/amdados/");
#endif
    file.open(filename, std::ios::out | std::ios::trunc);
    assert_true(file.good()) << "failed to open: " << filename << std::endl;
}

//-----------------------------------------------------------------------------
// Function queries the high-resolution clock until it has changed.
// The last obtained value is then a good seed guaranteed to be unique.
//-----------------------------------------------------------------------------
uint64_t RandomSeed()
{
    using hrclock_t = std::chrono::high_resolution_clock;
    uint64_t prev = hrclock_t::now().time_since_epoch().count();
    uint64_t seed = 0;
    while ((seed = hrclock_t::now().time_since_epoch().count()) == prev) {}
    return seed;
}

} // end namespace utils
} // end namespace app
} // end namespace allscale
