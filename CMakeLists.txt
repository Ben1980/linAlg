cmake_minimum_required(VERSION 3.14)
project(linAlg)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCE main.cpp matrix.h utils.h matrixFactory.h luDecomposition.h pivotLUDecomposition.h)
add_executable(${PROJECT_NAME} ${SOURCE})

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
find_package(fmt CONFIG REQUIRED)
find_package(doctest CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt fmt::fmt-header-only doctest::doctest)
