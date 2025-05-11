cmake_minimum_required(VERSION 3.29)
project(mpi_binomial_broadcast C)

set(CMAKE_C_STANDARD 17)

find_package(MPI REQUIRED)

add_executable(prog main.c)

target_link_libraries(prog PRIVATE MPI::MPI_C)