cmake_minimum_required(VERSION 2.8)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_C_STANDARD 11)

add_compile_options(-std=c++11 -Wextra -Werror -Wall -pedantic -O3 -fopenmp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CMAKE_SOURCE_DIR}/libs)
add_executable(particleSim main.cpp particle.cpp particle.hpp particles.cpp particles.hpp)
