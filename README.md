# amdados

Description goes here...

## Quickstart

Ensure you have GCC 6 (or 7) installed and set as your default C/C++ compiler.
Furthermore CMake 3.6 (or later) is required for the build and testing process.
Simply execute the following commands to build the project and run all tests.

    $ mkdir build
    $ cd build
    $ cmake ../code
    $ make -j8
    $ ctest -j8

## Advanced Options

### Configuration

Following options can be supplied to CMake

| Option                  | Values          |
| ----------------------- | --------------- |
| -DCMAKE_BUILD_TYPE      | Release / Debug |
| -DBUILD_SHARED_LIBS     | ON / OFF        |
| -DBUILD_TESTS           | ON / OFF        |
| -DBUILD_DOCS            | ON / OFF        |
| -DBUILD_COVERAGE        | ON / OFF        |
| -DUSE_ASSERT            | ON / OFF        |
| -DALLSCALE_CHECK_BOUNDS | ON / OFF        |
| -DUSE_VALGRIND          | ON / OFF        |
| -DUSE_ALLSCALECC        | ON / OFF        |
| -DENABLE_PROFILING      | ON / OFF        |
| -DTHIRD_PARTY_DIR       | \<path\>        |

The files `cmake/build_settings.cmake` and `code/CMakeLists.txt` state their
default value.

## Using the AllScale Compiler

To use the AllScale compiler, you first have to setup the required dependencies
in order to build it. Some packages should be installed with your system's
package manager, others must be installed via the provided dependency
installer.  See `scripts/dependencies/README.md` for more details. The
following commands assume you are using Ubuntu 16.04 LTS -- adjust the lines if
you are using a different distribution.

    $ sudo apt update
    $ sudo apt install bison build-essential cmake flex git groff-base \
      libcpufreq-dev libgmp10 libhwloc-dev libjemalloc-dev libluajit-5.1-dev \
      libpapi-dev m4 pkg-config python ruby time wget zlib1g
    $ scripts/dependencies/installer llvm boost cudd gmp zlib ghc cabal

To build this project using the AllScale compiler, simply set the corresponding
CMake option. You may want to use a separate build directory to easily switch
between GCC and AllScaleCC.

    $ mkdir build_with_allscalecc
    $ cd build_with_allscalecc
    $ ../scripts/dependencies/third_party_linker
    $ cmake -DUSE_ALLSCALECC=ON ../code
    $ make -j8
    $ ctest -j8

## Development

### Adding new Modules

The setup script can be run again to add new modules, just provide the same
project name.

    $ scripts/setup/run amdados frontend backend utils

### Adding new Parts to Modules

There is a utility script to add new *parts* to an existing module. The project
name and module name must be provided followed by a list of *parts* to
generate. Folders will be created along the way.

    $ scripts/setup/add_part amdados frontend sema extensions/malloc_extension

This will add the files `sema.h`, `sema.cpp` and `sema_test.cc` to the
*frontend* module. Furthermore new subfolders `extensions` will be created
containing `malloc_extension.h`, `malloc_extension.cpp` and
`malloc_extension_test.cc` in their respective subdirectories.

### Using a present AllScale Compiler / API project

The following two CMake options allow you to manually set the AllScale Compiler
and AllScale API to use. Note that `USE_ALLSCALECC` must *not* be enabled for
this to work.

| Option                  | Values                         |
| ----------------------- | ------------------------------ |
| -DOVERRIDE_ALLSCALECC   | *path to `allscalecc` binary*  |
| -DOVERRIDE_ALLSCALE_API | *path to AllScale API project* |

### Executable Bit

When working on Windows via SMB share, consider setting following Git setting.

    $ git config core.filemode false

### Licensor

A script, together with a Git hook, is provided to automatically add a license
header to each source file upon commit. See `scripts/license`.

### Eclipse Project

    $ cmake -G "Eclipse CDT4 - Unix Makefiles" /path/to/project

### Visual Studio Solution

    $ cmake -G "Visual Studio 15 2017 Win64" -DBUILD_SHARED_LIBS=OFF Z:\path\to\project

Add path for third-party libraries when needed.

## Troubleshooting

### Getting GCC 5 / CMake 3.5 / Valgrind (for Testing)

The dependency installer can setup these required tools for you. Its README
(`scripts/dependencies/README.md`) holds the details.

It is preferred to use the operating system's package manager, if applicable.

### No Source Folder in Eclipse Project

Make sure your build folder is located outside the source folder. Eclipse is
not capable of dealing with such a setup correctly.

### Coverage

Building the coverage us currently only supported on Linux, as Perl and Bash
are required. To build and view the coverage set the corresponding CMake flag
to `ON` and run:

    $ make
    $ make coverage
    $ xdg-open coverage/index.html
    
### Armadillo Library

The Armadillo matrix library is used only (!) for testing of matrix operations.
It is not included in the repository (see .gitignore).
Armadillo is installed along with other dependencies, see below.
As of 29.03.2018 Armadillo is included in the repository, but this is not
necessary. It can be always installed from the official web-site, see
./scripts/armadillo.sh

### Dependecies

In order to compile the project on machine without Internet connection
(e.g. blade behind firewall), some prerequisites should be downloaded
beforehand. From the project root folder, one should execute:
    bash ./scripts/download.sh
This will install all necessary dependencies (including the Armadillo library)
into the folder "api".

### Amdados application

For complete description see "./doc/amdados_report.pdf".
As a quick start, one can run the application with default settings,
section 2.4.1 in the report:

    bash standard.build.sh
    /bin/rm -fr output/*                    # optional cleanup
    python3 python/ObservationsGenerator.py --config amdados.conf
    ./build/app/amdados --config amdados.conf
    python3 python/Visualize.py --field output/field_Nx*_Ny*_Nt*.bin

where "field_Nx*_Ny_Nt*.bin" should be replaced with the actual file name with
prefix "field".

B E W A R E:
simulation might be very long (~ 1 day) on the machine with few CPU cores.

### Prerequisites

(1) Python of version 3.5+ is needed to generate observations, to find
    ground-truth solution and for testing & visualization,
    see "./doc/amdados_report.pdf".

    We recommend multiplatform package manager Anaconda ver.3:
    https://www.anaconda.com/download/

    Anaconda includes numpy, scipy, matplotlib, ffmpeg and their dependencies.

(2) The latest Allscale API and Armadillo library (used in unit tests)
    can be downloaded and installed as mentioned above:
        bash ./scripts/download.sh

(3) a) On Mac OS one probably wants to install "XQuartz" (via brew) for
    visualization with Python.
    b) Another useful tool is "mplayer" - a command-line player that
    shows the solution rolling out in time.

