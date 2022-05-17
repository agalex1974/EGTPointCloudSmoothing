cmake_minimum_required(VERSION 3.17)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 50)
endif()

project(EGTSmoothing CUDA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CGAL_DISABLE_GMP TRUE)
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

include(CheckLanguage)
check_language(CUDA)
find_package(OpenMP)

add_executable(EGTSmoothing main.cpp EGT.cu MCLSEXST.h Arrays.h CPointCloud.cpp CPointCloud.h CCreoPointCloud.cpp CCreoPointCloud.h EGT.cuh)

if(OpenMP_CXX_FOUND)
    target_link_libraries(EGTSmoothing PUBLIC OpenMP::OpenMP_CXX)
endif()

