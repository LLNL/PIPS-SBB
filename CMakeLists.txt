cmake_minimum_required(VERSION 2.8)
project(PIPS-SBB)

# CMake nested project boilerplate taken from
# http://techminded.net/blog/modular-c-projects-with-cmake.html
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# End boilerplate

# Right now, PIPS-SBB is contained in PIPS (parent directory), so
# so add parent directory.
# TODO: Restructure repository so this ugly construction is no longer
# needed.
set(PIPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/PIPS)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PIPS_DIR}/cmake/Modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")
message(STATUS "CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
set(BUILD_ALL OFF CACHE BOOL "Build all of PIPS (PIPS-S, PIPS-IPM, PIPS-NLP)")
set(BUILD_PIPS_S ON CACHE BOOL "Build PIPS-S")

add_subdirectory(${PIPS_DIR} ${CMAKE_BINARY_DIR}/PIPS)

# include different "whole archive" linking options depending on compiler
if (CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
  set(WHOLE_ARCHIVE "-Wl,-all_load")
  set(NO_WHOLE_ARCHIVE "-Wl,-noall_load")
  message(STATUS "SETTING HAVE_CLANG")
  set(HAVE_CLANG TRUE)

  # Clang doesn't have a Fortran compiler in its suite (yet),
  # so detect libraries for gfortran; we need equivalents to
  # libgfortran and libquadmath, which are implicitly
  # linked by flags in CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
  find_package(GFortranLibs REQUIRED)

  # Add directory containing libgfortran and libquadmath to
  # linker. Should also contain libgomp, if not using
  # Intel OpenMP runtime
  link_directories(${GFORTRAN_LIBRARIES_DIR})
else ()
  set(HAVE_CLANG FALSE)
  set(WHOLE_ARCHIVE "-Wl,--whole-archive")
  set(NO_WHOLE_ARCHIVE "-Wl,--no-whole-archive")
endif (CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")

# Check to see if using shared libraries
if (BUILD_SHARED_LIBS)
   message(STATUS "Using shared libraries")
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif (BUILD_SHARED_LIBS)

# inspired by elemental-bgp
if(MATH_LIBS)
	message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
elseif(IS_DIRECTORY $ENV{MKLROOT})
    set(MATH_LIBS "-Wl,--start-group  $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_sequential.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lgfortran -ldl")
    message(STATUS "Using MKLROOT:$ENV{MKLROOT} MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
    message(STATUS "Searching system path for MATH_LIBS")
    enable_language(Fortran)
    message(STATUS "FC = ${FC}")
    message(STATUS "Fortran compiler name is ${CMAKE_Fortran_COMPILER}")
    message(STATUS "Fortran compiler ID is ${CMAKE_Fortran_COMPILER_ID}")
    find_package(LAPACK REQUIRED)
    message(STATUS " LAPACK_LIBRARIES:  ${LAPACK_LIBRARIES}")
    message(STATUS " CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES:  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    set(MATH_LIBS ${LAPACK_LIBRARIES} -ldl -l${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ${OpenMP_CXX_FLAGS})
    message(STATUS " MATH_LIBS:  ${MATH_LIBS}")
endif(MATH_LIBS)

set(SHARED_DIR ${PIPS_DIR}/ThirdPartyLibs)

#CBC
set(COIN_DIR ${SHARED_DIR}/CBC/src)
if(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")
	message(FATAL_ERROR "Please run \"./wgetCBC.sh\" in ${COIN_DIR}")
endif(NOT EXISTS "${COIN_DIR}/lib/pkgconfig/coinutils.pc")

set(ENV{PKG_CONFIG_PATH} ${COIN_DIR}/lib/pkgconfig)
find_package(PkgConfig REQUIRED)
pkg_check_modules(COINUTILS REQUIRED coinutils)
pkg_check_modules(CLP REQUIRED clp)
pkg_check_modules(CBC REQUIRED cbc)
find_library(CLP_LIB NAMES libClp.a Clp PATHS ${COIN_DIR}/lib)

message(STATUS "Find all COINUTILS_LIBRARIES:${COINUTILS_LIBRARIES}")
set(COIN_LIBS "")
foreach(f ${COINUTILS_LIBRARIES})
	if(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack") # skip libm
	    message(STATUS "find ${f}...")
		find_library(X NAMES ${f} HINTS "${COIN_DIR}/lib" PATHS /usr/lib /usr/local/lib /usr/lib64/ /usr/lib/x86_64-linux-gnu/)
		message(STATUS "found ${f}:=${X}")
		set(COIN_LIBS ${COIN_LIBS} ${X})
		unset(X CACHE)
	endif(NOT ${f} MATCHES "m" AND NOT ${f} MATCHES "blas" AND NOT ${f} MATCHES "lapack")
endforeach(f)
message(STATUS "COIN_LIBS:=${COIN_LIBS}")

message(STATUS "CBC_LIBRARY_DIRS: ${CBC_LIBRARY_DIRS}")
find_library(CBCSolver_LIB NAMES libCbcSolver.a CbcSolver HINTS ${CBC_LIBRARY_DIRS})
find_library(CBC_LIB NAMES libCbc.a Cbc HINTS ${CBC_LIBRARY_DIRS})
find_library(OSICLP_LIB NAMES libOsiClp.a OsiClp HINTS ${CBC_LIBRARY_DIRS})
find_library(OSI_LIB NAMES libOsi.a Osi HINTS ${CBC_LIBRARY_DIRS})
find_library(CGL_LIB NAMES libCgl.a Cgl HINTS ${CBC_LIBRARY_DIRS})

set(CBC_LIBS ${CBCSolver_LIB} ${CBC_LIB} ${CGL_LIB} ${OSICLP_LIB} ${OSI_LIB})
include_directories(${COINUTILS_INCLUDE_DIRS})
message(STATUS "CBC_LIBS:=${CBC_LIBS}")
message(STATUS "CBC_INC:=${COINUTILS_INCLUDE_DIRS}")

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

include_directories(BBSMPSCuttingPlanes)
include_directories(BBSMPSBranching)
include_directories(BBSMPSTree)
include_directories(BBSMPSHeuristics)
include_directories(BBSMPSUtils)
include_directories(.)
add_library(pipssbb
	BBSMPSCuttingPlanes/BBSMPSCuttingPlanesManager.cpp
	BBSMPSCuttingPlanes/BBSMPSCuttingPlaneGenerator.cpp
	BBSMPSCuttingPlanes/BBSMPSCuttingPlaneGeneratorMIR.cpp
	BBSMPSCuttingPlanes/BBSMPSCuttingPlaneGeneratorGMI.cpp
	BBSMPSCuttingPlanes/BBSMPSCuttingPlane.cpp
	BBSMPSTree/BBSMPSTree.cpp
	BBSMPSTree/BBSMPSNode.cpp
	BBSMPSBranching/BBSMPSBranchingInfo.cpp
	BBSMPSBranching/BBSMPSBranchingRuleManager.cpp
	BBSMPSBranching/BBSMPSBranchingRule.cpp
	BBSMPSBranching/BBSMPSMaxFracBranchingRule.cpp
	BBSMPSBranching/BBSMPSPseudoCostBranchingRule.cpp
	BBSMPSUtils/BBSMPSUtils.cpp
	BBSMPSUtils/BBSMPSLogging.cpp
	BBSMPSHeuristics/BBSMPSHeuristicsManager.cpp
	BBSMPSHeuristics/BBSMPSHeuristic.cpp
	BBSMPSHeuristics/BBSMPSHeuristicRounding.cpp
	BBSMPSHeuristics/BBSMPSHeuristicRINS.cpp
	BBSMPSHeuristics/BBSMPSHeuristicRENS.cpp
	BBSMPSHeuristics/BBSMPSHeuristicFixAndDive.cpp
	BBSMPSHeuristics/BBSMPSHeuristicFixAndDiveLocks.cpp
	BBSMPSHeuristics/BBSMPSHeuristicCrossover.cpp
	BBSMPSHeuristics/BBSMPSHeuristicLockRounding.cpp
	BBSMPSHeuristics/BBSMPSHeuristicMagic.cpp
	BBSMPSHeuristics/BBSMPSHeuristicSolutionRINS.cpp
	BBSMPSHeuristics/BBSMPSHeuristicBestRINSJump.cpp
	BBSMPSHeuristics/BBSMPSHeuristicSolutionPolishing.cpp
	BBSMPSTree/BBSMPSSolver.cpp
	BBSMPSTree/BBSMPSSolverState.cpp
	BBSMPSTree/BBSMPSSolution.cpp)
add_executable(pipssbbSMPS pipssbbSMPS.cpp)
target_link_libraries(pipssbbSMPS pipssbb pipss stochInput ${COIN_LIBS} ${MATH_LIBS} ${Boost_LIBRARIES})
