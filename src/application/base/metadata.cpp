#include "metadata.h"

#include "logging.h"

// The TOML parser is a fairly large single-file include which
// could potentially bloat compilation times. We only include it
// from this .cpp file and not any of our header files.
#include "toml.hpp"

#include <fstream>

// Useful typedef for mapping between vec and toml types
using array3i = std::array<int, 3>;
using array3f = std::array<float, 3>;

vec3 DefaultColor = { 0.8f, 0.8f, 0.8f };

const std::string root_comment =
R"(# Overview
# --------
# This text file stores volume bounds and materials using the TOML markup format
# (see https://toml.io for specs and parsers). The actual voxels are stored in a
# separate binary file with a '.bin' or '.dag' extension.
#
# The format of a .bin file is very simple. The voxels are stored as a raw 3D
# array of 8-bit unsigned integers. The array dimensions can be determined from
# the volume bounds (see below). The array is laid out with the x-axis changing
# most quickly and the z-axis changing most slowly. Each voxel represents an
# index (0-255) into the material array (also below) describing the voxel.
#
# A .dag file is a more complex format but it can be converted to a .bin file
# using the Cubiquity voxel engine. See https://cubiquity.net for details.)";

const std::string bounds_comment =
R"(# Bounds
# ------
# These are *inclusive* bounds, so add one when computing the dimensions:
# 
#     dimensions[0] = (upper_bound[0] - lower_bound[0]) + 1
#     dimensions[1] = (upper_bound[1] - lower_bound[1]) + 1
#     dimensions[2] = (upper_bound[2] - lower_bound[2]) + 1
#
#     number_of_voxels = dimensions[0] * dimensions[1] * dimensions[2]
#
# The number of voxels should match the size of the .bin file.)";

const std::string materials_comment =
R"(# Materials
# ---------
# This array of tables contains the material properties. Each voxel is an 8-bit
# unsigned integer and so there can be up to 256 materials in the array below.)";

Metadata::Metadata()
{
	// This toml::value exists for the lifetime of the Metadata, though it's
	// contents may be reassigned. It is only deleted in the destructor.
	m_toml_data = new toml::ordered_value(toml::ordered_table());
	m_toml_data->comments().push_back(root_comment);
}

Metadata::~Metadata()
{
	delete m_toml_data;
}

/* ============================= Serialisation ============================== */

void Metadata::load(const std::filesystem::path& path)
{
	*(m_toml_data) = toml::parse<toml::ordered_type_config>(path);
}

void Metadata::save(const std::filesystem::path& path) const
{
	std::ofstream file(path);
	file.exceptions(~std::ofstream::goodbit); // Enable exceptions
	file << toml::format(*(m_toml_data));
}

/* ================================= Bounds ================================= */

void Metadata::ensure_bounds_exists()
{
	// This is a bit of a hack! The issue is that the 'toml11' library we are
	// using can store values in a toml::table (in which the order gets
	// randomised) ot a toml::ordered_table (in which the order is preserved)
	// but it does not allow explicit control over the ordering. The ordered
	// table can only have items added at the end, not inserted elsewhere.
	// 
	// However, for purely aesthetic reasons, I would like to ensure that the
	// bounds appear before the (potentially long) list of materials as this
	// makes the files more readable. But I don't want to force myself to always
	// set the bounds first. In fact, it is likely I want to set them last as I
	// may not know then until the volume has been populated.
	// 
	// The (dirty) solution here is to create a new TOML document, set the
	// bounds on it first, and then copy across any values from the original
	// document. Note that we only do this if the bounds don't already exist,
	// i.e. we avoid tampering with any existing documents in which the bounds
	// might have been manually placed at the end of the file by the user.

	bool has_bounds_table = false;
	if (m_toml_data->contains("bounds")) {
		if (m_toml_data->at("bounds").is_table()) {
			has_bounds_table = true;
		}
		else {
			// Presumably this could occur in the case of a broken input file.
			log_warning("Found 'bounds' but it is not a table. "
			            "It will be discarded and/or overwritten");
		}
	}

	if (has_bounds_table == false) {
		// Create new document with empty 'bounds' table
		toml::ordered_table new_data = toml::ordered_table();

		new_data["bounds"] = toml::ordered_table();
		new_data["bounds"].comments().push_back(bounds_comment);

		// Copy over any other values
		for (auto val : m_toml_data->as_table()) {
			new_data.insert(val);
		}

		// Replace the old document with the new one.
		*m_toml_data = new_data;
	}
}

std::optional<std::array<int, 3>> Metadata::lower_bound() const
{
	// Return std::optional as no sensible default.
	return toml::find<std::optional<std::array<int, 3>>>(
		*m_toml_data, "bounds", "lower_bound");
}

void Metadata::set_lower_bound(const std::array<int, 3>& lower_bound)
{
	ensure_bounds_exists();
	(*m_toml_data)["bounds"]["lower_bound"] = lower_bound;
}

std::optional<std::array<int, 3>> Metadata::upper_bound() const
{
	// Return std::optional as no sensible default.
	return toml::find<std::optional<std::array<int, 3>>>(
		*m_toml_data, "bounds", "upper_bound");
}

void Metadata::set_upper_bound(const std::array<int, 3>& upper_bound)
{
	ensure_bounds_exists();
	(*m_toml_data)["bounds"]["upper_bound"] = upper_bound;
}

/* =============================== Materials ================================ */

int Metadata::material_count() const
{
	// If the materials array doesn't exist then the
	// returned default has the expected size of zero.
	return toml::find_or(
		*(m_toml_data), "materials", toml::ordered_array()).size();
}

void Metadata::clear_all_materials()
{
	resize_materials(0);
}

void Metadata::resize_materials(int size)
{
	// Ensure there is a material array to resize. It does not exist by
	// default, and may not have been present in a loaded TOML file.
	bool has_material_array = false;
	if (m_toml_data->contains("materials")) {
		if (m_toml_data->at("materials").is_array()) {
			has_material_array = true;
		} else {
			// Presumably this could occur in the case of a broken input file.
			log_warning("Found 'materials' but it is not an array. "
			            "It will be discarded and/or overwritten");
		}
	}

	// Create the array if needed.
	if (has_material_array == false) {
		(*(m_toml_data))["materials"] = toml::ordered_array();
		(*(m_toml_data))["materials"].as_array_fmt().fmt =
			toml::array_format::array_of_tables;
	}

	// Set to our requested size
	toml::ordered_array& materials_array =
		m_toml_data->at("materials").as_array();
	materials_array.resize(size);

	// Any newly-created elements are in the 'empty' state,
	// whereas our materials should be an array of *tables*.
	for (int i = 0; i < materials_array.size(); i++) {
		if (materials_array[i].is_empty()) {
			materials_array[i] = toml::ordered_table();

			// We have a top-level comment describing the materials in general
			// as well as a per-material comment identifying its index. Ideally
			// we would attach the general description to the 'materials' array
			// of tables, but doing so seems to break the formatting (possibly
			// a bug in the TOML library we are using?). Therefore our (dirty)
			// solution is instead to prepend our general description to the
			// comment for the zeroth material.
			if (i == 0) {
				materials_array[i].comments().push_back(materials_comment);
				materials_array[i].comments().push_back(" "); // Blank line
			}
			materials_array[i].comments().push_back(
			    std::string(" Material index ") + std::to_string(i));
		}
	}
}

void Metadata::ensure_material_exists(int index)
{
	// When populating a material array I find it useful to be able to write to
	// any material index and assume it is valid. However, I don't want to
	// preallocate all 256 valid materials because in many cases most will be
	// empty, and so I want to avoid cluttering up the output file with empty
	// entries. My solution is simply to expand the materials array on demand.
	if (index >= material_count()) {
		resize_materials(index + 1);
	}
}

void Metadata::clear_material(int index)
{
	ensure_material_exists(index);

	// Create an empty table for the material
	(*(m_toml_data))["materials"][index] = toml::ordered_table();
}

void Metadata::set_material(int index, std::string name, vec3 base_color)
{
	clear_material(index); // Also ensures it exists
	set_material_name(index, name);
	set_material_base_color(index, base_color);
}

std::string Metadata::material_name(int index) const
{
	return toml::find_or< std::string >(
		*(m_toml_data), "materials", index, "name", "");
}

void Metadata::set_material_name(int index, std::string name)
{
	ensure_material_exists(index);
	(*(m_toml_data))["materials"][index]["name"] = name;
}

vec3 Metadata::material_base_color(int index) const
{
	return toml::find_or< array3f >(
		*(m_toml_data), "materials", index, "base_color", DefaultColor);
}

void Metadata::set_material_base_color(int index, vec3 color)
{
	ensure_material_exists(index);
	(*(m_toml_data))["materials"][index]["base_color"] = array3f(color);
}

void Metadata::set_material_to_empty_space(int index)
{
	// Note: Unlike 'set_material_to_default()' (below) we do not simply call 
	// 'set_material()' because we have no base color to set for empty space.
	clear_material(index); // Also ensures it exists
	set_material_name(index, "EmptySpace");
}

void Metadata::set_material_to_default(int index)
{
	set_material(index, "Default", DefaultColor);
}