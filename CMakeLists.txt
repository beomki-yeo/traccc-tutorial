# Set up the project.
cmake_minimum_required( VERSION 3.22 )
project( momentum_resolution LANGUAGES CXX )

# Set up the used C++ standard(s).
set( CMAKE_CXX_STANDARD 20 CACHE STRING "The (host) C++ standard to use" )

find_package( Boost 1.86.0 REQUIRED COMPONENTS program_options filesystem log)

# Include covfie and traccc
add_subdirectory(extern/covfie)
add_subdirectory(extern/traccc)

# Track parameter estimation
add_executable( track_parameter_estimation tutorials/track_parameter_estimation.cpp )
target_link_libraries( track_parameter_estimation traccc::core )