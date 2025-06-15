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

cmake_minimum_required(VERSION 3.15)

project (Cubiquity)

################################################################################
#
# All targets
#
################################################################################

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Lots of warnings and all warnings as errors
if (MSVC)
    add_compile_options(/W2) # Can increase this to W4 after I fix the most severe warnings.
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

if(APPLE OR MINGW)
    # MinGW does not properly support std::execution::par (it does compile, but
    # does not actually parallelise. Using PoolSTL is a temporary workaround.
    add_compile_definitions(CUBIQUITY_USE_POOLSTL)
endif()

# Temporary solution to get C++ filesystem support through CMake:
# See https://gitlab.kitware.com/cmake/cmake/-/issues/17834
if (NOT MSVC AND NOT APPLE)
    link_libraries(stdc++fs)
endif()

# Needed to group CMake predefined targets.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
#
# Main library
#
################################################################################

# Eventually we may want to allow the user to build shared or static version of Cubiquity?
# See https://cmake.org/cmake/help/latest/variable/BUILD_SHARED_LIBS.html
file(GLOB LIBRARY_SRCS src/library/*)
add_library(Library ${LIBRARY_SRCS})

# I'm not quite clear on why we need this for a standard library feature but it
# seems that we do (at least on GCC?). See https://stackoverflow.com/a/39547577
find_package(Threads REQUIRED)
target_link_libraries(Library Threads::Threads)

if(UNIX)
    find_package(TBB REQUIRED)
    # GCC parallel algorithms depend on Intel TBB on Linux. This is not easily
    # available on MinGW, but we use PoolSTL as a substitute for that.
    target_link_libraries(Library TBB::tbb)
endif()

if(MINGW)
    # Statically link (as Windows users are unlikely to have MinGW runtimes?).
    target_link_libraries(Library -static)
endif()

# Required to profile release builds.
if(MSVC)
    target_link_options(Library PUBLIC /PROFILE)
    target_compile_options(Library PUBLIC /Zi)
endif()

# Note: I would like to uncomment the line below, but doing so seems to break the  
# expected dependency between the library being updated and a relink being required.
# At least withing Visual Studio. More investigation required.
# This indicates it should work: https://discourse.cmake.org/t/how-to-use-the-same-name-for-a-target-and-a-library/5044
#set_target_properties(Library PROPERTIES OUTPUT_NAME cubiquity)

# For tidier folder layout in Visual Studio (mimics application target below)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src/library PREFIX "Source Code" FILES ${LIBRARY_SRCS})

################################################################################
#
# Frontend Application
#
################################################################################

# SDL 2 is needed for the view subcommand. Info on using SDL via CMake is here:
#   https://wiki.libsdl.org/SDL2/README/cmake
find_package(SDL2 CONFIG)

# Source code for most of the application
file(GLOB APPLICATION_SRCS
    src/application/*
    src/application/base/*
    src/application/commands/export/*
    src/application/commands/export/vox_writer/vox_writer.* # Skip example.cpp
    src/application/commands/generate/*
    src/application/commands/test/*
    src/application/commands/voxelize/*
    src/application/external/*)

# Source files for optional 'view' subcommand (requires SDL)
if (SDL2_FOUND)
    file(GLOB APPLICATION_VIEW_SRCS
        src/application/commands/view/*
        src/application/commands/view/framework/*)
endif()

# Configuration of base application
add_executable(Application ${APPLICATION_SRCS} ${APPLICATION_VIEW_SRCS})
include_directories(src/library src/application)
include_directories(src/application/external) # Needed so external deps can find own headers
target_link_libraries(Application Library ${CMAKE_DL_LIBS})
set_target_properties(Application PROPERTIES OUTPUT_NAME cubiquity)

# Extra config if we are building the 'view' subcommand.
if (SDL2_FOUND)
    add_compile_definitions(CUBIQUITY_APP_ENABLE_VIEW)
    include_directories(${SDL2_INCLUDE_DIRS})
    target_link_libraries(Application SDL2::SDL2main SDL2::SDL2-static)
    if (WIN32)
        # '-mconsole' is needed to get stdout/stderr output on Windows
        # See https://stackoverflow.com/a/54079540
        target_link_libraries(Application -mconsole)
    endif()
endif()

# For tidier folder layout in Visual Studio
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src/application PREFIX "Source Code" FILES ${APPLICATION_SRCS} ${APPLICATION_VIEW_SRCS})

install(TARGETS Application DESTINATION .)
install(FILES COPYING.txt DESTINATION .)
install(DIRECTORY src/application/commands/view/glsl DESTINATION .)
