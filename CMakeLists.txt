cmake_minimum_required(VERSION 2.8)
project(PIPS-SBB)

# CMake nested project boilerplate taken from
# http://techminded.net/blog/modular-c-projects-with-cmake.html
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# End boilerplate

# Check to see if using shared libraries
if (BUILD_SHARED_LIBS)
   message(STATUS "Using shared libraries")
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif (BUILD_SHARED_LIBS)

# Macros for PIPS
add_definitions(-DTIMING -DSTOCH_TESTING)

# Include math libraries if needed; inspired by elemental-bgp
if(MATH_LIBS)
        set(NEED_MATH FALSE)
        message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
        set(NEED_MATH TRUE)
endif(MATH_LIBS)

if(NEED_MATH)
        if(IS_DIRECTORY $ENV{MKLROOT})
                set(MATH_LIBS "-Wl,--start-group  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a
$ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lgfortran")
        else(IS_DIRECTORY $ENV{MKLROOT})
                enable_language( Fortran )
                find_package(LAPACK REQUIRED)
                set(MATH_LIBS ${LAPACK_LIBRARIES};${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
        endif(IS_DIRECTORY $ENV{MKLROOT})
endif(NEED_MATH)

# Right now, PIPS-SBB is contained in PIPS (parent directory), so
# so add parent directory.
# TODO: Restructure repository so this ugly construction is no longer
# needed.
set(PIPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/PIPS)
add_subdirectory(${PIPS_DIR} ${CMAKE_BINARY_DIR}/PIPS)

set(SHARED_DIR ${PIPS_DIR}/SharedLibraries)

# include CoinUtils
set(COIN_DIR ${SHARED_DIR}/Cbc-2.7.6)

if(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")
        message(FATAL_ERROR "Please run \"./configure --enable-static; make install\" in PIPS/SharedLibraries/Cbc-2.7.6")
endif(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")

set(ENV{PKG_CONFIG_PATH} ${COIN_DIR}/lib/pkgconfig)
find_package(PkgConfig REQUIRED)
pkg_check_modules(COINUTILS REQUIRED coinutils)

set(COIN_LIBS "")

message(Status "COINUTILS LIBRARIES ${COINUTILS_LIBRARIES} ${ZLIB_LIBRARY}")
foreach(f ${COINUTILS_LIBRARIES})
        if(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack") # skip libm
                message(STATUS "looking after ${f}")
                find_library(X NAMES ${f} HINTS "${COINUTILS_LIBRARY_DIRS}" PATHS /usr/lib /usr/local/lib /usr/lib64/)
                message(STATUS "found ${X}")
                set(COIN_LIBS ${COIN_LIBS} ${X})
                unset(X CACHE)
        endif(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack")
endforeach(f)
#message(STATUS "Coin libraries: ${COIN_LIBS}")

find_package(Boost COMPONENTS system filesystem thread log log_setup)
include_directories(${Boost_INCLUDE_DIR})
if (${Boost_FOUND}$ MATCHES "FALSE")
  message("! A Boost version containing BOOST.LOG has NOT been found.")
  message("! Please specify the location of the Boost library containing 'log' package.")
  message("! For this create a 'Toolchain.cmake' file or edit the one from the root ")
  message("! directory and use 'set(BOOST_ROOT \"/path/to/boost\")'.")
  message("! Then run  'cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake ..' in the build directory.")
  message(FATAL_ERROR "Configuration aborted due to missing BOOST libraries!")
endif(${Boost_FOUND}$ MATCHES "FALSE")

#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
#set(CMAKE_BUILD_TYPE RELEASE)
#message("FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")

#IF(CMAKE_COMPILER_IS_GNUCXX)
#       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)

include_directories(${COINUTILS_INCLUDE_DIRS})

# Include directories for PIPS-SBB
include_directories(${PIPS_DIR}/PIPS-S/CoinBALPFactorization)
include_directories(${PIPS_DIR}/PIPS-S/Basic)
include_directories(${PIPS_DIR}/PIPS-S/Core)
include_directories(${PIPS_DIR}/PIPS-S/Drivers)
include_directories(${PIPS_DIR}/SolverInterface)
include_directories(${PIPS_DIR}/Input)


include_directories(BBSMPSBranching)
include_directories(BBSMPSTree)
include_directories(BBSMPSHeuristics)
include_directories(BBSMPSUtils)
include_directories(.)
add_library(pipssbb BBSMPSTree/BBSMPSTree.cpp 
	BBSMPSTree/BBSMPSNode.cpp 
	BBSMPSBranching/BBSMPSBranchingInfo.cpp
	BBSMPSBranching/BBSMPSBranchingRuleManager.cpp 
	BBSMPSBranching/BBSMPSBranchingRule.cpp 
	BBSMPSBranching/BBSMPSMaxFracBranchingRule.cpp 
	BBSMPSUtils/BBSMPSUtils.cpp
	BBSMPSUtils/BBSMPSLogging.cpp
	BBSMPSHeuristics/BBSMPSHeuristicsManager.cpp 
	BBSMPSHeuristics/BBSMPSHeuristic.cpp 
	BBSMPSHeuristics/BBSMPSHeuristicRounding.cpp
	BBSMPSHeuristics/BBSMPSHeuristicRINS.cpp
	BBSMPSHeuristics/BBSMPSHeuristicRENS.cpp
	BBSMPSHeuristics/BBSMPSHeuristicFixAndDive.cpp
	BBSMPSTree/BBSMPSSolver.cpp 
	BBSMPSTree/BBSMPSSolverState.cpp
	BBSMPSTree/BBSMPSSolution.cpp)
add_executable(pipssbbSMPS pipssbbSMPS.cpp)
target_link_libraries(pipssbbSMPS pipssbb pipss stochInput ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})	
