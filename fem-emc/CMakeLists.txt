cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(FemGauntlet LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(BUILD_SHARED_LIBS ON)

include_directories(.)


# Add the directory with your submission
add_subdirectory(project)


add_executable(
    solver
    ./main.c
)

find_library(m m)
target_link_libraries(
    solver 
    PUBLIC femsolver 
)
target_link_libraries(solver PUBLIC m)
