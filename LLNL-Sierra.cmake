SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_Fortran_COMPILER mpif90)

SET(Boost_NO_BOOST_CMAKE TRUE)
SET(BOOST_ROOT /usr/gapps/ppp/boost-nompi-1.57.0)
SET(Boost_LIBRARY_DIRS /usr/gapps/ppp/boost-nompi-1.57.0/lib)
SET(Boost_INCLUDE_DIR /usr/gapps/ppp/boost-nompi-1.57.0/include)