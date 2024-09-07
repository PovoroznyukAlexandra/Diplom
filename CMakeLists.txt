cmake_minimum_required(VERSION 3.20)
project(Diplom)

set(GCC_COVERAGE_COMPILE_FLAGS "-flto -Ofast -march=native -llapack")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

find_package(LAPACK REQUIRED)


add_executable(Diplom diplom_main.cpp diplom_func.h diplom_func.cpp)
