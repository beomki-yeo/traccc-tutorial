# TRACCC tutorial for beginners
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required( VERSION 3.22 )
include( FetchContent )

# Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# Tell the user what's happening.
message( STATUS "Building Traccc as part of the project" )

# Declare where to get Traccc from.
set( TRACCC_SOURCE  
   "GIT_REPOSITORY;https://github.com/acts-project/traccc;GIT_TAG;c06d483980feeebc7f2c874255f8d89617771365"
   CACHE STRING "Source for Traccc, when built as part of this project" )

mark_as_advanced( TRACCC_SOURCE )

# Mark the import as a system library on modern CMake versions
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
   set(TRACCC_SOURCE_FULL "${TRACCC_SOURCE};SYSTEM")
else()
   set(TRACCC_SOURCE_FULL "${TRACCC_SOURCE}")
endif()
mark_as_advanced( TRACCC_SOURCE_FULL )

FetchContent_Declare( Traccc ${TRACCC_SOURCE_FULL} )

# Options used in the build of Detray.
set( TRACCC_CUSTOM_SCALARTYPE "double" CACHE STRING
   "Scalar type to use in the Detray code" )
set( DETRAY_CUSTOM_SCALARTYPE "double" CACHE STRING
   "Scalar type to use in the Detray code" )
set( TRACCC_BUILD_IO TRUE CACHE BOOL "Turn on the IO build" )
set( TRACCC_BUILD_TESTING FALSE CACHE BOOL "Turn off the Test build" )
set( TRACCC_BUILD_BENCHMARKS FALSE CACHE BOOL "Turn off benchmark build" )
set( TRACCC_USE_ROOT FALSE CACHE BOOL "Turn off ROOT" )
set( TRACCC_USE_ROOT FALSE CACHE BOOL "Turn off ROOT" )
set( TRACCC_BUILD_EXAMPLES TRUE CACHE BOOL "Turn on example build" )
set( TRACCC_SETUP_COVFIE FALSE CACHE BOOL
   "Do not set up covfie as part of Traccc" )

# Get it into the current directory.
FetchContent_MakeAvailable( Traccc )
