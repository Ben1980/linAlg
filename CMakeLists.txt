cmake_minimum_required(VERSION 3.14)
project(linAlg)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#add_compile_options(-O3)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#add_compile_definitions(_ASMATRIX)
add_compile_definitions(_ASARRAY)

set(SOURCE main.cpp matrix.h utils.h matrixFactory.h luDecomposition.h pivotLUDecomposition.h)
add_executable(${PROJECT_NAME} ${SOURCE})

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
find_package(fmt CONFIG REQUIRED)
find_package(doctest CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE fmt::fmt fmt::fmt-header-only doctest::doctest)
