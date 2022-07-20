cmake_minimum_required(VERSION 3.12)

project(AMIQS VERSION 1.0)


# compiler
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
# set(CMAKE_CXX_FLAGS "-L/usr/local/lib")
# with -Werror warning will be raised to error, for example unused variable will be an error..
# add_definitions(-Wall -Wextra -Werror)
# add_definitions(-Wall)



# Configuration.
find_package(Threads REQUIRED)
find_package(GSL REQUIRED)


# set executable directory
# could change the amiqs.exe directory to the amiss folder but then would need to change 
# the relative paths in the .cpp files
# set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)


# targets
add_library( lsoda ${CMAKE_CURRENT_SOURCE_DIR}/Source/LSODA.cpp)
include_directories(${GSL_INCLUDE_DIR})


#stef, amiqs
add_executable( amiqs.exe ${CMAKE_CURRENT_SOURCE_DIR}/amiqs.cpp)
target_link_libraries(amiqs.exe lsoda)
target_link_libraries(amiqs.exe GSL::gsl GSL::gslcblas)



# Tests and benchmark
# enable_testing()
#stef 
# add_test(NAME amiqs.exe COMMAND $<TARGET_FILE:amiqs>)

