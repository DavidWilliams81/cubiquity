Build notes
===========

Prerequisites
-------------
Install [CMake](https://cmake.org/) and [SDL2](https://libsdl.org/). On Linux these can typically be found in your package manager.

Linux with GCC
--------------
```
mkdir build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release ..
make
```

Windows with Visual Studio
--------------------------
* Download and install CMake
* Download and build SDL
* Run CMake Gui, run configure, point it at SDL install folder, generate.
* Build in Visual Studio.

Note: Cubiquity defaults to using static SDL which is not actually present in the SDL prebuilts on GitHub. Hence is is necessary to build our own version. It must also be installed (via the INSTALL target in Visual Studio) in order to be useable from the Cubiquity CMake scripts.

Cross-compile to Windows from Ubuntu with GCC
--------------------------------------------
Tested on Ubuntu 24.04. Note that there are significant performance problem with this build which appear to be at least partly due to Intel Threaded Building Blocks for MinGW not being installed. This needs a review.

Run:

    sudo apt install mingw-w64*

Save the following snippet as 'mingw-w64-x86_64.cmake' (taken from [here](https://gist.github.com/peterspackman/8cf73f7f12ba270aa8192d6911972fe8#file-mingw-w64-x86_64-cmake))

```
# Sample toolchain file for building for Windows from an Ubuntu Linux system.
#
# Typical usage:
#    *) install cross compiler: `sudo apt-get install mingw-w64`
#    *) cd build
#    *) cmake -DCMAKE_TOOLCHAIN_FILE=~/mingw-w64-x86_64.cmake ..
# This is free and unencumbered software released into the public domain.

set(CMAKE_SYSTEM_NAME Windows)
set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)

# cross compilers to use for C, C++ and Fortran
set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
set(CMAKE_Fortran_COMPILER ${TOOLCHAIN_PREFIX}-gfortran)
set(CMAKE_RC_COMPILER ${TOOLCHAIN_PREFIX}-windres)

# target environment on the build host system
set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})

# modify default behavior of FIND_XXX() commands
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
```

Download and extract MinGW prebuilt SDL. This is presumably contains Windows files but it works in the Linux cross-compilation toolchain. Note: The MinGW build of SDL 2.30 has broken CMake support. It is (or will be) fixed in 2.31. See this bug report: https://github.com/libsdl-org/SDL/issues/11795

Make a build folder within Cubiquity:

```
mkdir build-mingw
cd build-mingw/
```

Run CMake, passing it the location of the toolchain and SDL. E.g.:

    cmake -DCMAKE_TOOLCHAIN_FILE=~/Downloads/mingw-w64-x86_64.cmake -DSDL2_DIR=/home/davidw/Downloads/SDL2-2.30.10/cmake/ ..

Run ```make```

Further reading on cross-compiling:
* https://cmake.org/cmake/help/book/mastering-cmake/chapter/Cross%20Compiling%20With%20CMake.html
* https://stackoverflow.com/questions/63178407/cmake-compile-in-linux-execute-in-windows
* https://stackoverflow.com/questions/73612230/cross-compiling-c-with-sdl-on-linux
