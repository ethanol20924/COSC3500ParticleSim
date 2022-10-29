cmake_minimum_required(VERSION 3.8)
project(particleSim LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

target_compile_features(particleSim PUBLIC cxx_std_11)
set_target_properties(particleSim PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(particleSim)

add_compile_options(-std=c++11 -Wextra -Werror -Wall -pedantic -O3)

include_directories(${CMAKE_SOURCE_DIR}/libs)
add_executable(particleSim main.cpp particle.cpp particle.hpp particles.cu particles.hpp)
