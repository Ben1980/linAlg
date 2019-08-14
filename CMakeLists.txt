cmake_minimum_required(VERSION 3.13)
project(linAlg LANGUAGES CXX)

add_executable(linAlg main.cpp matrix.h)

find_package(fmt CONFIG REQUIRED)
find_package(doctest CONFIG REQUIRED)
target_link_libraries(linAlg PRIVATE fmt::fmt fmt::fmt-header-only doctest::doctest)

set_property(TARGET linAlg PROPERTY CXX_STANDARD 17)