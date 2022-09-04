cmake_minimum_required(VERSION 2.8)

# set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_C_STANDARD 11)

add_compile_options(-std=c++11 -Wall -Wextra -pedantic -Werror -O3)

include_directories(${CMAKE_SOURCE_DIR}/libs)
add_executable(particleSim main.cpp particle.cpp particle.hpp particles.cpp particles.hpp)
