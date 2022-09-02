cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 11)

add_compile_options(-Wall -Wextra -pedantic -Werror)

include_directories(${CMAKE_SOURCE_DIR}/libs)
add_executable(particleSim main.cpp particle.cpp particle.hpp particles.cpp particles.hpp)
