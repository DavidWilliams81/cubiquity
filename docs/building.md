Build notes
===========
Cubiquity is straight-forward to build (especially on Linux) due to it's small size and minimal dependencies. It uses [CMake](https://cmake.org/) to generate a platform-specific build system such as a Makefile under Linux or a Visual Studio solution under Windows. CMake can be found in the package manager of most Linux distributions, or a Windows installer can be found on the CMake homepage.

Dependencies
------------
Cubiquity consists of a library implementing the main data structures and algorithms, and a wrapper application which exposes some of the functionality via a command line interface. The core library has no dependencies beyond the C++ standard library, while the application has a number or lightweight dependencies which are included in the source tree and hence do not have to be installed.

The only external dependency for the wrapper application is the [Simple Direct Media Layer (SDL, Version 2)](https://libsdl.org/) which is used to display the window used by the Cubiquity 'view' command. This is an optional dependency - if it is not available then Cubiquity still be built but the 'view' command will not work. Other commands such as 'voxelize' and 'export' still work without SDL2 being present

Linux with G++
--------------
These instructions are tested on Ubuntu 24.10. You may need to adapt them for you own system.

1. Install the required dependencies:

    `sudo apt install g++ make cmake libsdl2-dev`

2. Create a build directory in the root of the Cubiquity repository and change into it:

    ```
    mkdir build
    cd build/
    ```

3. Run CMake to generate a Makefile.

    `cmake -DCMAKE_BUILD_TYPE=Release ..`

4. Make Cubiquity

    `make`

This should result in an executable called 'cubiquity' in the current ('build') directory.

Windows with Visual Studio
--------------------------
The steps required for Windows are broadly similar to those for Linux (above) but are typically performed via the CMake GUI and Visual Studio IDE rather than from the command line (though the command line is still possible).

These steps were tested on Windows 11 using Visual Studio 2022 and CMake version 3.31.4.

1. Install Visual Studio (making sure to include the C++ components) and CMake.

2. Obtain SDL (version 2) if you wish to include the 'view' command in the resulting Cubiquity application.

   * **Note:** SDL2 Currently needs to be obtained by building it from source. This is because Cubiquity is configured to statically link against it, and the SDL2 prebuilt binaries do not include a static version. Full details are too much to include here (see the [SDL page on CMake](https://wiki.libsdl.org/SDL2/README/cmake)) but basically you should download the SDL source code and build it using CMake. It must also be installed (via the INSTALL target in Visual Studio) in order to be usable from the Cubiquity CMake scripts.

3. Use the CMake GUI to generate a Visual Studio solution:

   * Run the CMake GUI and provide it with the path to the local Cubiquity repository. You also need to provide it with a build folder, which should be the same path again but with a 'build' subdirectory. Click 'Configure', and let CMake create the build folder as it does not yet exist.

   * Accept the default options for creating a Visual Studio solution

   * CMake will probably fail to find SDL. You need to edit the SDL2_DIR to point at the 'cmake' subdirectory within the SDL2 install location (not the location of the SDL2 source code, which also has a 'cmake' folder!).

   * If you want to statically link against the Visual C++ runtime then add a new CMake variable called 'CMAKE_MSVC_RUNTIME_LIBRARY' and set it to 'MultiThreaded' (as opposed to 'MultiThreadedDLL').

   * Click 'Configure' again after making the above changes.

   * Click 'Generate' to create the Visual Studio solution.

   * Click 'Open project' to open the solution file which was just generated.

4. Build Cubiquity in Visual Studio (you may wish to select the 'Release' configuration).


Cross-compile to Windows from Linux with G++
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

Download and extract MinGW prebuilt SDL. This is presumably contains Windows files but it works in the Linux cross-compilation toolchain. Note: The MinGW build of SDL 2.30.10 has broken CMake support. It is (or will be) fixed in later versions. See this bug report: https://github.com/libsdl-org/SDL/issues/11795

Make a build folder within Cubiquity:

```
mkdir build-mingw
cd build-mingw/
```

Run CMake, passing it the location of the toolchain and SDL. E.g.:

    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=~/Downloads/mingw-w64-x86_64.cmake -DSDL2_DIR=/home/davidw/Downloads/SDL2-2.31.0/cmake/ ..

Run ```make```

Further reading on cross-compiling:
* https://cmake.org/cmake/help/book/mastering-cmake/chapter/Cross%20Compiling%20With%20CMake.html
* https://stackoverflow.com/questions/63178407/cmake-compile-in-linux-execute-in-windows
* https://stackoverflow.com/questions/73612230/cross-compiling-c-with-sdl-on-linux
