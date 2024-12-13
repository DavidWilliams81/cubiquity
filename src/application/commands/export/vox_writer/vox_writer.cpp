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

// Based on specs here:
// 
// https://github.com/ephtracy/voxel-model/blob/master/
//
// Specs give version number as 150 but latest MagicaVoxel (0.99.7.1) uses 200.
// Some things are not clearly specified:
//     - Is a layer chunk needed?
//     - Does the order of chunks matter? At least it determines the model id.
//     - Does the transform tree need to be the same as the chank layout
//       (should chunks be nested)? As far as I can tell, no. The  hierarchy is
//       only determined by the node id, child id, model, id, etc.
//     - Does node zero have to be the root, or is the root the first transform
//       node in the file?
//     - Although chunks have a chunk size and child size, the main chunk
//       appears to be the only one which has children. All other chunks are a
//       flat list within this.
// Given some lack of clarity I have also made test scenes in MagicaVoxel and
// viewed them in a hex editor. When I have doubts about the spec I just try to
// do what MagicaVoxel does.

#include "vox_writer.h"

#include <algorithm>
#include <execution>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <vector>

// Maximum size for MagicaVoxel models is 256.
const int model_size_power = 8;
const int model_size = 1 << model_size_power;

// The vox_writer creates a 3D grid of models overlaid on the volume data to be
// written. Each model starts at a multiple of model_size and extends to the
// next multiple of model_size minus 1. This means that models can extend beyond
// the volume bounds, and this can sometimes present as bounding boxes which 
// appear excessively large in MagicaVoxel.
// 
// In contrast, the clipped bounds represent the model bounds clipped to the
// volume bounds. These are used to ensure we never actually sample the volume
// outside of the region specified by the user. Ideally I might make the clipped
// bounds show up as the the bounding boxes in MagicaVoxel but I've found this
// suprisingly tricky to achieve, and it doesn't really matter.
struct model
{
	int id;
	vox_writer::box bounds;
	vox_writer::box clipped_bounds;
};

// Implementation of FNV-1a hash function
uint64_t fnv1a_hash_64bit(const uint8_t* data, size_t size)
{
	uint64_t hash = UINT64_C(0xcbf29ce484222325);
	for (size_t i = 0; i < size; i++) {
		hash ^= data[i];
		hash *= UINT64_C(0x00000100000001B3);
	}
	return hash;
}

vox_writer::vox_writer()
{
	// Not too bright for UI but still visible if used by mistake.
	const color unused_col = { 102, 102, 102, 255 }; // 40% grey
	std::fill(m_palette.begin(), m_palette.end(), unused_col);
}

void vox_writer::set_palette(int index, const color& col)
{
	if (index < 1 || index > 255) {
		// Note that .vox does not store a palette entry for index zero.
		throw std::out_of_range("palette index out of range [1-255]");
	}

	// .vox stores palette offset by one
	m_palette[index - 1] = col;
}

void vox_writer::on_progress(int /*done*/, int /*total*/)
{
	// Empty - can optionally be implemented by subclasses.
}

void vox_writer::write_char(char data)
{
	m_file.write(&data, sizeof(data));
}

void vox_writer::write_int32(int32_t data)
{
	m_file.write(reinterpret_cast<const char*>(&data), sizeof(data));
}

void vox_writer::write_string(const std::string& data)
{
	write_int32(data.size());        // String length
	m_file.write(reinterpret_cast<const char*>(
		data.c_str()), data.size()); // String data
}

void vox_writer::write_id(const std::array<char, 4>& id)
{
	write_char(id[0]);
	write_char(id[1]);
	write_char(id[2]);
	write_char(id[3]);
}

void vox_writer::write_size(const box& model_bounds)
{
	write_id({ 'S', 'I', 'Z', 'E' });
	write_int32(12); // Chunk bytes
	write_int32(0);  // Child bytes

	write_int32(model_bounds.upper.x - model_bounds.lower.x + 1); // X size
	write_int32(model_bounds.upper.y - model_bounds.lower.y + 1); // Y size
	write_int32(model_bounds.upper.z - model_bounds.lower.z + 1); // Z size
}

void vox_writer::write_xyzi(const std::vector<int32_t>& voxels)
{
	write_id({ 'X', 'Y', 'Z', 'I' });

	const int voxel_size = sizeof(voxels[0]); // 4
	const int32_t num_voxels = voxels.size();
	write_int32(sizeof(num_voxels) + (num_voxels * voxel_size)); // Chunk bytes
	write_int32(0);                                              // Child bytes

	write_int32(num_voxels);
	m_file.write(reinterpret_cast<const char*>(
		voxels.data()), voxels.size() * voxel_size);
}

void vox_writer::write_rgba()
{
	write_id({ 'R', 'G', 'B', 'A' });

	write_int32(sizeof(m_palette)); // Chunk bytes
	write_int32(0);                 // Child bytes

	// No need to apply offset by one as this was done in set_palette().
	m_file.write(reinterpret_cast<const char*>(
		m_palette.data()), sizeof(m_palette));
}

void vox_writer::write_ntrn(int32_t node_id, int32_t child_id, const vec3i& pos)
{
	write_id({ 'n', 'T', 'R', 'N' });

	auto chunk_bytes_pos = m_file.tellp();
	write_int32(0);  // Placeholder for chunk bytes
	write_int32(0);  // Child bytes

	auto data_pos = m_file.tellp();
	write_int32(node_id);
	write_int32(0);  // Empty node attributes dict
	write_int32(child_id);
	write_int32(-1); // Reserved
	write_int32(0);  // Layer id (layer 0 seems to exist automatically?)
	write_int32(1);  // Number of frames

	write_int32(1); // Translation is only entry
	write_string("_t");
	std::stringstream ss;
	ss << pos.x << " " << pos.y << " " << pos.z;
	write_string(ss.str());

	auto chunk_bytes = m_file.tellp() - data_pos;

	// Go back and fill in chunk bytes
	m_file.seekp(chunk_bytes_pos);
	write_int32(chunk_bytes);
	m_file.seekp(0, std::ios::end);
}

void vox_writer::write_ngrp(int32_t node_id,
	                        const std::vector<int32_t>& child_ids)
{
	write_id({ 'n', 'G', 'R', 'P' });

	write_int32(12 + child_ids.size() * 4);     // Chunk bytes
	write_int32(0);                             // Child bytes

	write_int32(node_id);
	write_int32(0);                             // Empty node attributes dict

	write_int32(child_ids.size());              // No of children
	m_file.write(reinterpret_cast<const char*>( // Actual ids
		child_ids.data()), child_ids.size() * sizeof(child_ids[0]));
}

void  vox_writer::write_nshp(int32_t node_id, int32_t model_id)
{
	write_id({ 'n', 'S', 'H', 'P' });

	write_int32(20); // Chunk bytes
	write_int32(0);  // Child bytes

	write_int32(node_id);
	write_int32(0);  // Empty node attributes dict
	write_int32(1);  // No of models

	write_int32(model_id);
	write_int32(0);  // Empty model attributes dict
}

void vox_writer::write_main()
{
	write_id({ 'M', 'A', 'I', 'N' });
	write_int32(0); // Chunk bytes

	auto child_bytes_pos = m_file.tellp();
	write_int32(0); // Placeholder for child bytes

	auto data_pos = m_file.tellp();

	// ====== Build a 3D grid of models covering the volume bounds ======
	box vol_bounds = bounds();
	int32_t mask = ~(model_size - 1); // Converts position to model origin
	vec3i first_origin = { vol_bounds.lower.x & mask, 
		                   vol_bounds.lower.y & mask, 
		                   vol_bounds.lower.z & mask };
	vec3i last_origin  = { vol_bounds.upper.x & mask,
		                   vol_bounds.upper.y & mask,
		                   vol_bounds.upper.z & mask };

	std::vector<model> models;
	for (int z = first_origin.z; z <= last_origin.z; z += model_size) {
		for (int y = first_origin.y; y <= last_origin.y; y += model_size) {
			for (int x = first_origin.x; x <= last_origin.x; x += model_size) {

				model mdl;
				mdl.bounds.lower = { x, y, z };
				mdl.bounds.upper = { x + model_size - 1,
					                 y + model_size - 1,
					                 z + model_size - 1 };

				mdl.clipped_bounds.lower.x =
					std::max(mdl.bounds.lower.x, vol_bounds.lower.x);
				mdl.clipped_bounds.lower.y =
					std::max(mdl.bounds.lower.y, vol_bounds.lower.y);
				mdl.clipped_bounds.lower.z =
					std::max(mdl.bounds.lower.z, vol_bounds.lower.z);
				mdl.clipped_bounds.upper.x =
					std::min(mdl.bounds.upper.x, vol_bounds.upper.x);
				mdl.clipped_bounds.upper.y =
					std::min(mdl.bounds.upper.y, vol_bounds.upper.y);
				mdl.clipped_bounds.upper.z =
					std::min(mdl.bounds.upper.z, vol_bounds.upper.z);

				mdl.id = -1; // Indicates invalid, may be updated later.
				models.push_back(mdl);
			}
		}
	}

	// ====== Write the models to the file (possibly in parallel) ======

	// Initialise progress prior to processing
	int done = 0;
	on_progress(done, models.size());

	// The main model writing code is a lambda function which
	// can be called sequentialy or in parallel for each model.
	std::mutex mut;
	int model_id = 0;
	std::unordered_map<uint64_t, int> models_ids; // Map hash to model id.
	auto write_model = [&](model& mdl)       // Called later for_each() model
	{
		// ---------------- Begin parallel code ----------------

		// Get a list of voxels for the model
		std::vector<int32_t> voxels;
		const auto& lower = mdl.clipped_bounds.lower;
		const auto& upper = mdl.clipped_bounds.upper;
		for (int z = lower.z; z <= upper.z; z++) {
			for (int y = lower.y; y <= upper.y; y++) {
				for (int x = lower.x; x <= upper.x; x++) {

					int i = voxel({ x, y, z });
					if (i > 0) {
						voxels.push_back((x - mdl.bounds.lower.x) |
							((y - mdl.bounds.lower.y) << 8 ) |
							((z - mdl.bounds.lower.z) << 16) |
							 (i << 24));
					}
				}
			}
		}

		std::lock_guard<std::mutex> guard(mut);
		// ----------------  End parallel code  ----------------

		// If the model contained any voxels then write it out.
		if (voxels.size() > 0)
		{
			// Check whether the model is the same as one we have previously
			// written out. If so we can reuse the previous model to save space.
			// This is particularly effective for filled high-resolution objects
			// which may have large homoneous regions.
			// We identify matching models by their hash and do not handle
			// collisions. For a few thousand models and a 64-bit hash the
			// chances of a collision are approx one in a trillion. See:
			// https://preshing.com/20110504/hash-collision-probabilities/
			uint64_t hash = fnv1a_hash_64bit(reinterpret_cast<uint8_t*>(
				voxels.data()), sizeof(voxels[0]) * voxels.size());

			if (models_ids.find(hash) != models_ids.end()) {
				mdl.id = models_ids[hash];
			} else {
				write_size(mdl.bounds);
				write_xyzi(voxels);
				mdl.id = model_id;
				models_ids[hash] = model_id;
				model_id++;
			}
		}

		// Update progress
		done++;
		on_progress(done, models.size());
	};

	// Voxel access can be done in parallel if user allows
	if (m_multithreaded) {
		std::for_each(std::execution::par,
			models.begin(), models.end(), write_model);
	} else {
		std::for_each(std::execution::seq,
			models.begin(), models.end(), write_model);
	}

	// ======  Write the transform chunks to position the model ======

	int id = 2; // Offset accounts for root transform and group nodes.
	std::vector<int> transform_ids;
	std::vector<int> shape_ids;
	for (auto& m : models) {
		if (m.id != -1) {
			transform_ids.push_back(id++); // Gets 2,4,6,...
			shape_ids.push_back(id++);     // Gets 3,5,7,...
		}
	}

	// Write root transform and group (nodes 0 and 1). We have to
	// offset by half the model size to align with the MagicaVoxel
	// origin, which I find a little strange but it seesm to work.
	write_ntrn(0, 1, { model_size / 2, model_size / 2, model_size / 2 });
	write_ngrp(1, transform_ids);

	// Write out a transform and shape for each model.
	unsigned int i = 0;
	for (auto& m : models) {
		if (m.id != -1) {
			write_ntrn(transform_ids[i], shape_ids[i], m.bounds.lower);
			write_nshp(shape_ids[i], m.id);
			i++;
		}
	}

	write_rgba();
	auto child_bytes = m_file.tellp() - data_pos;

	// MagicaVoxel uses 32-bit signed ints to represent chunk/child sizes.
	// This limits the size to 2 GB but only the main chunk ever gets so big.
	if (child_bytes > std::numeric_limits<int32_t>::max()) {
		throw std::runtime_error("child data in main chunk exceeds 2GB limit");
	}

	// Go back and fill in child bytes
	m_file.seekp(child_bytes_pos);
	write_int32(child_bytes);
	m_file.seekp(0, std::ios::end);
}

void vox_writer::write(const std::filesystem::path& filename,
	                   bool multithreaded)
{
	try {
		m_file = std::ofstream(filename, std::ios::out | std::ios::binary);
		m_multithreaded = multithreaded;

		write_id({ 'V', 'O', 'X', ' ' });
		write_int32(150); // Version
		write_main();
		m_file.close(); // Explicit, as m_file stays in scope
	} catch (std::exception& e) {
		// Close and reopen to truncate to zero bytes.
		m_file.close();
		m_file = std::ofstream(filename, std::ios::out | std::ios::binary);
		m_file.close();
		throw e; // Rethrow
	}
}
