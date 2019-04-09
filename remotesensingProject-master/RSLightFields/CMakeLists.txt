cmake_minimum_required( VERSION 3.7 )

## RSLightFields

# The project is composed of a library and some tests that make use of this library.
# The library sources and headers are located under /src and /include.
# The test sources are located under /tests.


# Configure project
project( RSLightFields )
set( PROJECT_DESCRIPTION "Implements some depth map estimator from 3D light fields" )

# Set C++ compiler flags.
set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lstdc++fs -fopenmp" )

# Specify library content :
include_directories( include )
file( GLOB_RECURSE ${PROJECT_NAME}_HEADERS include/*.hpp )

# Build library
file( GLOB LIB_SOURCES src/*.cpp )
file( GLOB LIB_HEADERS include/*.hpp )
add_library( RSLightFields ${LIB_SOURCES} ${LIB_HEADERS} )

# Add and link OpenCV to the library
find_package( OpenCV REQUIRED )
target_link_libraries( RSLightFields ${OpenCV_LIBS} )
target_link_libraries( RSLightFields stdc++fs )

# Build tests and link them to the library
file( GLOB TEST_SOURCES tests/*.cpp )

foreach( TEST_SOURCE_FILE ${TEST_SOURCES} )

	get_filename_component( TEST_NAME ${TEST_SOURCE_FILE} NAME_WE )
	add_executable( ${TEST_NAME} ${TEST_SOURCE_FILE} )
	target_link_libraries( ${TEST_NAME} RSLightFields )

endforeach( TEST_SOURCE_FILE ${TEST_SOURCES} )


