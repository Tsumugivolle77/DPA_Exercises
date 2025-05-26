cmake_minimum_required(VERSION 3.10.0)
project(parallel_algo_demo VERSION 0.1.0 LANGUAGES C)

set (TARGET_NAME "quick_sort_teng")

add_executable(${TARGET_NAME} p1.c)

add_compile_options(-fopenmp)

include_directories(/opt/homebrew/Caskroom/miniforge/base/include/)

find_package(OpenMP REQUIRED)
if(OpenMP_C_FOUND)
    target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_C)
endif()