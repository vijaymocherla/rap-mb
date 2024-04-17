if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()

cmake_policy(SET CMP0048 NEW)  # project_VERSION* variables populated from project(... VERSION x.x.x) string

project(rap_mb
VERSION 0.1.0
LANGUAGES CXX)
set(rap_mb_AUTHORS      "Sai Vijay Bhaskar Mocherla")
set(rap_mb_DESCRIPTION  "")
set(rap_mb_URL          "https://github.com/vijaymocherla/rap-mb")
set(rap_mb_LICENSE      "MIT")

include_directories(.)
aux_source_directory(. SRC_LIST)

set(CMAKE_PREFIX_PATH ${CONDA_PREFIX})

# Uncomment the following lines to force cmake to build rap_mb with icpx using mkl
set(CMAKE_C_COMPILER icx)
set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_FC_COMPILER ifx)
set(OpenMP_FIND_COMPONENTS "CXX")



set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# Checking MKL
find_package(MKL)
if(MKL_FOUND)
include_directories(${MKL_INCLUDE_DIR})
add_definitions(-DEIGEN_USE_MKL=1)
set(EIGEN_USE_MKL_ALL 1) # This will go into config.h - it makes Eigen use MKL
else()
set(EIGEN_USE_MKL_ALL 0)
endif()
# Getting package dependencies
find_package(OpenMP REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)



if(CODE_COVERAGE)
  message("-- Code coverage enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
endif(CODE_COVERAGE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++1z-extensions -lstdc++ -lsundials_core -lsundials_arkode -lsundials_nvecserial -qopenmp -lstdc++ -qmkl") # avoid warnings for C++17


pybind11_add_module(_rap_mb ${SRC_LIST})
