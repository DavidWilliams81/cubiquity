####################################################################################################
# Cubiquity - A micro-voxel engine for games and other interactive applications                    #
#                                                                                                  #
# Written in 2019 by David Williams                                                                #
#                                                                                                  #
# To the extent possible under law, the author(s) have dedicated all copyright and related and     #
# neighboring rights to this software to the public domain worldwide. This software is distributed #
# without any warranty.                                                                            #
#                                                                                                  #
# You should have received a copy of the CC0 Public Domain Dedication along with this software.    #
# If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   #
####################################################################################################
cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

project (Cubiquity)

################################################################################
#
# All targets
#
################################################################################

# C++17 required for flags.h library.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Lots of warnings and all warnings as errors
if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# Additional include directories
include_directories(src/cubiquity)
include_directories(src/external)
include_directories(src/utility)

# SDL 2 is needed for the demos. See https://stackoverflow.com/a/44347594
# and https://trenki2.github.io/blog/2017/06/02/using-sdl2-with-cmake/
find_package(SDL2 REQUIRED)

# I'm not quite clear on why we need this for a standard library feature but it
# seems that we do (at least on GCC?). See https://stackoverflow.com/a/39547577
find_package(Threads REQUIRED)

# Temporary solution to get C++ filesystem support through CMake:
# See https://gitlab.kitware.com/cmake/cmake/-/issues/17834
if (NOT MSVC)
    link_libraries(stdc++fs)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
#
# Main library
#
################################################################################

file(GLOB CUBIQUITY_SRCS src/cubiquity/*)
add_library(cubiquity ${CUBIQUITY_SRCS})
target_link_libraries(cubiquity Threads::Threads)

################################################################################
#
# External dependencies for tests, tools, and demos.
#
################################################################################

file(GLOB EXTERNAL_SRCS src/external/*)
add_library(external ${EXTERNAL_SRCS})
target_link_libraries(external)

################################################################################
#
# Utility code for tests, tools, and demos.
#
################################################################################

file(GLOB UTILITY_SRCS src/utility/*)
add_library(utility ${UTILITY_SRCS})
target_link_libraries(utility cubiquity external ${SDL2_LIBRARIES} ${CMAKE_DL_LIBS})

################################################################################
#
# Instancing demo
#
################################################################################

file(GLOB INSTANCING_SRCS
     src/demos/instancing/*
     src/demos/instancing/framework/*)
add_executable(instancing ${INSTANCING_SRCS})
include_directories(instancing ${SDL2_INCLUDE_DIRS})
target_link_libraries(instancing cubiquity external utility)

################################################################################
#
# Pathtracing demo
#
################################################################################

add_executable(pathtracing src/demos/pathtracing/main.cpp)
include_directories(pathtracing ${SDL2_INCLUDE_DIRS})
target_link_libraries(pathtracing cubiquity external utility ${SDL2_LIBRARIES})

# Set startup project in Visual Studio: https://stackoverflow.com/a/37994396
set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT pathtracing)

################################################################################
#
# Tests
#
################################################################################

file(GLOB TEST_SRCS src/tests/*.cpp src/tests/*.h src/tests/*.inl)
add_executable(test ${TEST_SRCS})
target_link_libraries(test cubiquity external utility)

################################################################################
#
# Voxelise tool
#
################################################################################

add_executable(voxelize src/tools/voxelize/voxelize.cpp)
target_link_libraries(voxelize cubiquity)
