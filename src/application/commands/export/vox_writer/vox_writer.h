/*******************************************************************************
* vox_writer - A library for writing volume data in the MagicaVoxel format.    *
*                                                                              *
* Written in 2024 by David Williams                                            *
*                                                                              *
* To the extent possible under law, the author(s) have dedicated all copyright *
* and related and neighboring rights to this software to the public domain     *
* worldwide. This software is distributed without any warranty.                *
*                                                                              *
* You should have received a copy of the CC0 Public Domain Dedication along    *
* with this software.                                                          *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.               *
*******************************************************************************/

// Basic Usage
// -----------
// Begin by subclassing vox_writer to wrap your own data structure (or your
// procedural generation algorithm) and by implementing the bounds() and voxel()
// functions. You can then instantiate your subclass and call the write()
// function. The libary will call your provided voxel() function for each
// position within the range returned by bounds(), in order to obtain the data
// which it writes out.
// 
// Each voxel in the .vox format is an 8-bit index into a color palette (note
// that more advanced MagicaVoxel materials are not currently supported). You
// can specify the RGBA color for relevent palette slots using set_palette(),
// and you are responsible for mapping your voxel data to an appropriate 8-bit 
// index within your implementation of voxel().
//
// You can monitor progress by overriding the on_progress() function.
//
// Further details
// ---------------
// A MagicalVoxel model has a maximum size of 256 voxels along each side,
// but multiple models can be placed in a scene. This vox_writer works by
// building a 3D grid of such models to cover the occupied parts of your volume.
//
// The vox_writer can operate in single-threaded (the default) or multi-threaded
// mode. When running in multi-threaded mode both voxel() and on_progress() can
// be called concurrently from multiple threads. Therefore you need to ensure
// that your implementations of these are thread safe. When using this library
// it is likly that your voxel() function will be the bottleneck, so enabling
// multi-threaded operation will probably give a significant performance boost.
//
// Limitations
// -----------
// This library is intended solely for writing out individual large volues
// and does not expose the full range of MagicaVoxel capabilities. Amoung
// other things it is is missing:
//     - Key-frame animation
//     - Scene hierarchy management
//     - Advanced materials
//
// Known issues
// ------------
// The bounding boxes of each model are aligned to multiple of 256 regardless of
// the contents of the model. This means they are often larger than they need to
// be. This can look strange when viewing the models in edit mode but does not
// affect the final rendered image.
// 
// Future work
// -----------
//     - Tighten model bounding boxes
//     - Key-frame animation support
//     - Split large files into multiple parts
//     - Support for advanced material parameters
//     - Let user specify model size (always power-of-two)
//     - Alert caller if scene bounds exceed MagicaVoxel limits

#ifndef CUBIQUITY_VOX_WRITER_H
#define CUBIQUITY_VOX_WRITER_H

#include <array>
#include <filesystem>
#include <fstream>
#include <vector>

// Subclass this class to implement your MagicaVoxel file writer.
class vox_writer
{
public:
	// Simple utility structures
	struct vec3i { int32_t x; int32_t y; int32_t z; };
	struct box   { vec3i lower; vec3i upper; };
	struct color { uint8_t r; uint8_t g; uint8_t b; uint8_t a; };

	// Set the RGBA color for the specified index (1-255).
	void set_palette(int index, const color& col);

	// Write the volume data to the specified file. This function calls voxel()
	// for each position with bounds() to obtain the data to write, and splits
	// it into as many models as necessary. If the multithreaded flag is set
	// then calls to voxel() and on_progress() can occur from multiple threads
	// concurrently, otherwise all calls occur from the current thread.
	void write(const std::filesystem::path& filename,
		       bool multithreaded = false);

protected:

	// Protected constructor and destructor as this class
	// cannot be instantiated directly (must be subclassed).
	vox_writer();
	virtual ~vox_writer() {};

	// You must override this function to return the inclusive bounds of the
	// volume. It is called from the same thread as write(). The voxel()
	// function will be called for each position within these bounds.
	virtual box     bounds() = 0;
 
	// You must override this function to return the 8-bit colour index for the
	// voxel at the specified position. Note that it may be called concurrently
	// from multiple threads in multithreaded mode.
	virtual uint8_t voxel(const vec3i& position) = 0;

	// Optionally override this function to monitor progress. Note that it may
	// be called concurrently from multiple threads in multithreaded mode.
	virtual void    on_progress(int done, int total);

private:
	// Write basic types
	void write_char (char data);
	void write_int32(int32_t data);
	void write_string(const std::string& data);

	// Write standard MagicaVoxel chunks
	void write_id(const std::array<char, 4>& id);
	void write_size(const box& model_bounds);
	void write_xyzi(const std::vector<int32_t>& voxels);
	void write_rgba();
	void write_main();

	// Write extended MagicaVoxel chunks (for model hierarchy)
	void write_ntrn(int32_t node_id, int32_t child_id, const vec3i& pos);
	void write_ngrp(int32_t node_id, const std::vector<int32_t>& child_ids);
	void write_nshp(int32_t node_id, int32_t model_id);

	// Writer configuration
	std::ofstream m_file;
	bool m_multithreaded = false;
	std::array<color, 256> m_palette;
};

#endif // CUBIQUITY_VOX_WRITER_H
