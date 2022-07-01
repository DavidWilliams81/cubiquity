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
    add_compile_options(/W2) # Can increase this to W4 after I fix the most severe warnings.
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# SDL 2 is needed for the demos. See https://stackoverflow.com/a/44347594
# and https://trenki2.github.io/blog/2017/06/02/using-sdl2-with-cmake/
find_package(SDL2)

# I'm not quite clear on why we need this for a standard library feature but it
# seems that we do (at least on GCC?). See https://stackoverflow.com/a/39547577
find_package(Threads REQUIRED)

# Temporary solution to get C++ filesystem support through CMake:
# See https://gitlab.kitware.com/cmake/cmake/-/issues/17834
if (NOT MSVC)
    link_libraries(stdc++fs)
endif()

# Needed to group CMake predefined targets.
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

################################################################################
#
# Main library
#
################################################################################

file(GLOB LIBRARY_SRCS src/library/*)
add_library(Library ${LIBRARY_SRCS})
target_link_libraries(Library Threads::Threads)

# Required to profile release bulds.
if(MSVC)
    target_link_options(Library PUBLIC /PROFILE)
    target_compile_options(Library PUBLIC /Zi)
endif()

# Note: I would like to uncomment the line below, but doing so seems to break the  
# expected dependancy between the library being updated and a relink being required. 
# At least withing Visual Studio. More investigation required. 
#set_target_properties(Library PROPERTIES OUTPUT_NAME cubiquity)

# For tidier folder layout in Visual Studio (mimics application target below)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src/library PREFIX "Source Code" FILES ${LIBRARY_SRCS})

################################################################################
#
# Frontend Application
#
################################################################################

if (SDL2_FOUND)
    file(GLOB APPLICATION_SRCS
        src/application/*
        src/application/base/*
        src/application/commands/export/*
        src/application/commands/generate/*
        src/application/commands/test/*
        src/application/commands/view/*
        src/application/commands/view/framework/*
        src/application/commands/voxelise/*
        src/application/external/*)
    add_executable(Application ${APPLICATION_SRCS})
    include_directories(src/library src/application ${SDL2_INCLUDE_DIRS})
    include_directories(src/application/external) # Needed so external deps can find own headers
    target_link_libraries(Application Library ${SDL2_LIBRARIES} ${CMAKE_DL_LIBS})
    set_target_properties(Application PROPERTIES OUTPUT_NAME cubiquity)

    # For tidier folder layout in Visual Studio
    source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src/application PREFIX "Source Code" FILES ${APPLICATION_SRCS})
endif()
