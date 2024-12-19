# TRACCC tutorial for beginners
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the project.
cmake_minimum_required( VERSION 3.22 )
project( momentum_resolution LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 20 CACHE STRING "The (host) C++ standard to use" )

find_package( Boost 1.86.0 REQUIRED COMPONENTS program_options filesystem log)

# Include covfie and traccc
add_subdirectory(extern/covfie)
add_subdirectory(extern/traccc)

# Clusterization
add_executable( clusterization tutorials/clusterization.cpp )
target_link_libraries( clusterization traccc::core traccc::simulation )

# Seeding
add_executable( seeding tutorials/seeding.cpp )
target_link_libraries( seeding traccc::core )