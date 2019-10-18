cmake_minimum_required(VERSION 3.14)
project(linAlg)

set(SOURCE main.cpp matrix.h utils.h matrixFactory.h luDecomposition.h)

add_executable(${PROJECT_NAME} ${SOURCE})
#add_compile_definitions(_ASMATRIX)
add_compile_definitions(_ASARRAY)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
find_package(fmt CONFIG REQUIRED)
find_package(doctest CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt fmt::fmt-header-only doctest::doctest)
