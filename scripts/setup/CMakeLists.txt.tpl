cmake_minimum_required(VERSION 3.5)
project(%PROJECT% VERSION 0.0.0 LANGUAGES C CXX)

# -- Module Path
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# -- Prefix Path
file(GLOB prefix_paths ${PROJECT_SOURCE_DIR}/third_party/*)
list(APPEND CMAKE_PREFIX_PATH ${prefix_paths})

# -- Project Settings
include(build_settings)
include(boost_settings)
include(doxygen)

# -- Dependencies
#set(GMP_VERSION 6.0.0 CACHE STRING "GMP Version")

# -- CMake Modules
include(file_globs)
include(add_unittest)
include(msvc_file_completion)
include(allscale)

# -- Project Modules
