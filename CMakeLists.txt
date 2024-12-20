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

# Include traccc
add_subdirectory(extern/traccc)

# Clusterization
add_executable( clusterization tutorials/clusterization.cpp )
target_link_libraries( clusterization traccc::core traccc::simulation )

# Seeding
add_executable( seeding tutorials/seeding.cpp )
target_link_libraries( seeding traccc::core )

# Track finding
add_executable( track_finding tutorials/track_finding.cpp )
target_link_libraries( track_finding traccc::core traccc::simulation )

# Track fitting
add_executable( track_fitting tutorials/track_fitting.cpp )
target_link_libraries( track_fitting traccc::core traccc::simulation )

# Telescope detector writer
add_executable( write_detector tutorials/write_detector.cpp )
target_link_libraries( write_detector traccc::core traccc::simulation traccc::io )