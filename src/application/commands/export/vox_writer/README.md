# vox_writer - Writing high-resolution MagicaVoxel files from C++

This library provides a simple way to write high-resolution [MagicaVoxel](https://ephtracy.github.io/) files from a C++ program. You can use it to wrap your existing data structure for exporting to MagicaVoxel, or you can generate voxel scenes procedurally.

The library is part of the [Cubiquity Voxel Engine](https://github.com/DavidWilliams81/cubiquity) but can be used completely independently. The source code is released into the public domain under the [CC0 dedication](https://creativecommons.org/publicdomain/zero/1.0/).

<p align="center">
	<a href="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/header.png"><img src="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/header.png" alt="Screenshots"></a>
	<i>This scene was modelled using assets from <a href="https://quaternius.com/">Quaternius</a>, voxelized using <a href="https://github.com/DavidWilliams81/cubiquity">Cubiquity</a>, and exported to the <a href="https://ephtracy.github.io/">MagicaVoxel</a> format using this library.</i>
</p>

## Contents
* [Features](#features)
* [Limitations](#limitations)
* [Usage](#usage)
* [Viewing high-resolution files](#viewing-high-resolution-files)

## Features
**High-resolution:** There are no limits on the dimensions of the volume which can be written (though see [the note below](#limitations) on max *file size*). The library takes care of splitting the volume into a number of models and building a scene graph to accommodate the constraints of the .vox format.

**Performant:** Data can be read directly from your own data structure and is written straight to disk. There are no intermediate memory structures or allocations required. 

**Multithreaded** The writer can operate in multithreaded mode if your underlying data structure or algorithm supports concurrent access to voxels.

**Easy to use:** Provided as a single pair of .cpp/.h files which you can add directly to your project. No knowledge of the .vox format is required. Just subclass the vox_writer and implement a couple of functions to specify the bounds and contents of the volume.

**Cross-platform:** Written in C++17 and has no dependencies apart from the C++ standard library.

**Public domain:** The source code is released into the public domain via the [CC0 dedication](https://creativecommons.org/publicdomain/zero/1.0/). It is completely free for any purpose.

## Limitations
The library only supports writing out a single large volume. It does not support animation or explicit scene graph management. These are outside the intended scope of the library and there is no plan to add them.

The library only supports an indexed RGBA color per voxel. Advanced MagicaVoxel materials (metal, emissive, etc) are not currently supported though it is possible they could be added in the future.

Although this library has no limit on the dimensions of the volume, the .vox format has a 2 GB file size limit. Furthermore, MagicaVoxel can struggle to render large volumes. See the section ['Viewing high-resolution files'](#Viewing-high-resolution-files) below for tips on viewing large volumes.

## Usage
It is easy to integrate this library with your application. It is provided as a single pair of .cpp/.h files which you can add directly to your project or makefile. You can then subclass `vox_writer` to implement a couple of functions as shown below (from [example.cpp](example.cpp)):

```c++
#include "vox_writer.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// Subclass the vox_writer to wrap your own data structure or implement
// procedural generation (shown here). You need to override the bounds()
// and voxel() functions, and optionally the on_progress() function.
class example_vox_writer : public vox_writer
{
public:
	// Perform any required initialisation in the constructor
	example_vox_writer(int radius) : m_radius(radius)
	{
		// Set up the palette (slot 0 is empty space and can't be set).
		set_palette(1, { 0xf3, 0xfa, 0xe1, 0xff }); // Beige
		set_palette(2, { 0xf7, 0xf6, 0xc5, 0xff }); // Cream
		set_palette(3, { 0xfa, 0xb2, 0xb8, 0xff }); // Light pink
		set_palette(4, { 0xfc, 0x6d, 0xab, 0xff }); // Darker pink
		set_palette(5, { 0xde, 0x5d, 0xd4, 0xff }); // Towards purple
		set_palette(6, { 0xc0, 0x4c, 0xfd, 0xff }); // Purple
		set_palette(7, { 0x8f, 0x3c, 0xfe, 0xff }); // Towards indigo
		set_palette(8, { 0x5e, 0x2b, 0xff, 0xff }); // Bright indigo
	}

	// Return the inclusive bounds of your volume data.
	box bounds() override
	{
		// These bounds are chosen to cut the sphere in half.
		return { {-m_radius, -m_radius, 0},          // Lower bound
		         { m_radius,  m_radius, m_radius} }; // Upper bound
	}

	// Get the color index for the voxel at the specified position. This
	// is called for every voxel inside the region defined by bounds().
	uint8_t voxel(const vec3i& position) override
	{
		// Raise the sphere up to sit on ground plane
		int new_z = position.z - m_radius;

		// Choose a color index based on distance from center
		const int color_bands = 8;
		float dist = sqrt((position.x * position.x) + 
		                  (position.y * position.y) +
		                  (new_z * new_z));
		dist /= m_radius;  // Normalize (0.0 is center, 1.0 is edge)
		dist = 1.0 - dist; // Invert (1.0 is center, 0.0 is edge)
		int index = static_cast<int>(dist * color_bands + 1);
		return static_cast<uint8_t>(std::clamp(index, 0, color_bands));
	}

	// Override this function to monitor progress.
	void on_progress(int done, int total) override
	{
		std::cout << "Done " << done << " of " << total << std::endl;
	}

private:
	const int m_radius;
};

int main(void)
{
	std::filesystem::path filename = "example.vox";
	example_vox_writer writer(200); // Create the writer
	writer.write(filename);         // Write the file to disk
	std::cout << "Finished writing " << filename << std::endl;
	return 0;
}
```
The resulting example.vox looks as shown below when opened in MagicaVoxel:

<p align="center"><a href="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/example.png"><img src="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/example.png" alt="Output from the above example program"></a></p>

### Building the example
#### Linux
On Linux you can build this example with GCC by running:

```
make
```

Internally this runs:

```
g++ example.cpp vox_writer.cpp -std=c++17 -lstdc++ -ltbb -o example
```

Note that C++17 support is required for std::filesystem and std::for_each(). The GCC implementation of parallel algorithms also depends on Intel Threaded Building Blocks (see the bottom of [this page](https://gcc.gnu.org/onlinedocs/libstdc++/manual/using.html)), hence the need for the -ltbb flag above. You may be able to do without this on other compilers.

#### Windows
There is no makefile for building the example under Windows. The code is so simple that you can just make your own e.g. Visual Studio project if you want to test it. To use the library itself you just add the vox_writer.cpp/.h to your project.

### Running the example
You can run the example to generate the output volume:

```
./example
```
gives:

```
Done 0 of 4
Done 1 of 4
Done 2 of 4
Done 3 of 4
Done 4 of 4
Finshed writing "example.vox"
```

## Viewing high-resolution files
There are a couple of options for viewing the high-resolution .vox files which can be generated by this library:

### MagicaVoxel
[MagicaVoxel](https://ephtracy.github.io/) is probably the most well-known voxel editor and the program for which the .vox format was originally created. However, despite its popularity it is limited in the size of scene it can open and render. Therefore you might find that when you open a scene only a portion of the volume is actually displayed. There are two things you can do to improve this:

#### 1. Increase render buffer size.
Open the following file:
```
...\path\to\MagicaVoxel\config\config.txt
```
And set the following:
```
render :
{
    dense_buffer   : '1024'   // [16, 1024] MB
    sparse_buffer  : '1024'   // [16, 1024] MB
}
```
This maximizes the number of voxels which can be loaded. 

#### 2. Switch to 'sparse' mode for rendering.
This can be found by switching to the 'Render' tab, finding the 'Sample' settings, scrolling down to the 'Geometry' group, and enabling the 'Sparse' toggle. See the screenshot below to make it clearer:

<p align="center"><a href="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/sparse.png"><img src="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/sparse.png" alt="Sparse setting"></a></p>

Note that enabling sparse mode can trigger a significant pause (perhaps a minute or more) while MagicaVoxel recalculates the data.

### Avoyd
[Avoyd](https://www.avoyd.com/) is a powerful voxel editor with a free version available. It supports very high resolution volumes (up to 256k voxels on each side) so can easily handle the files output from vox_writer.

The .vox files can be opened via `File > Open` or `File > Import > MagicaVoxel (.vox)`, or by dragging and dropping a .vox file into Avoyd. Further information is available [here](https://www.avoyd.com/avoyd-voxel-editor-documentation.html#ImportVox).

<p align="center"><a href="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/avoyd.png"><img src="https://s3.eu-west-2.amazonaws.com/dpw81.public/vox_writer/avoyd.png" alt="MagicaVoxel file imported into Avoyd"></a></p>
