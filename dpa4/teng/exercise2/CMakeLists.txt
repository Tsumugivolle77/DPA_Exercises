cmake_minimum_required(VERSION 3.10.0)
project(dpa4-nr2 VERSION 0.1.0 LANGUAGES C)

add_executable(dpa4-nr2 main.c mmio.c)

add_compile_options(-fopenmp -O3 -Wall -Wextra -Wpedantic)

find_package(OpenMP REQUIRED)
if(OpenMP_C_FOUND)
    target_link_libraries(dpa4-nr2 PUBLIC OpenMP::OpenMP_C)
endif()