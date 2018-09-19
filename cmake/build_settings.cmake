option(BUILD_SHARED_LIBS "Link libraries dynamically" ON)
option(BUILD_TESTS "Enable testing" ON)
option(BUILD_DOCS "Enable documentation" OFF)
option(USE_ASSERT "Enable assertions" ON)
option(USE_VALGRIND "Allow Valgrind for unit tests" OFF)
option(USE_ALLSCALECC "Use allscalecc as compiler" OFF)
option(ENABLE_PROFILING "Enable AllScale profiling support" OFF)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "CMake Build Type" FORCE)
endif()

if(BUILD_TESTS)
	enable_testing()
endif()

if(ENABLE_PROFILING)
	add_definitions(-DENABLE_PROFILING)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
	# C flags
	set(CMAKE_C_FLAGS "-Wall -Wextra -std=c99")
	set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -ggdb")
	set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")

	# C++ flags
	set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++14")
	set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -ggdb")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

	# Yo Dawg, I heard you like templates!
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-900")

	# \[T]/ Praise the sun. There be dragons ahead.
	if(NOT USE_ASSERT)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DNDEBUG")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG")
	endif()

	if(NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
		# allow arbitrary library linking order (in case `-as-needed` is default)
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
		set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -Wl,--no-as-needed")
		set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")
	endif()
elseif(MSVC)
	include(msvc_settings)
	set(USE_VALGRIND OFF)
else()
	message(FATAL_ERROR "Unhandled Compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Forward settings for external projects.
set(CMAKE_EXTERNALPROJECT_FORWARDS
	-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	-DCMAKE_C_COMPILER_ARG1=${CMAKE_C_COMPILER_ARG1}
	-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DCMAKE_CXX_COMPILER_ARG1=${CMAKE_CXX_COMPILER_ARG1}
	-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
	-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
)
