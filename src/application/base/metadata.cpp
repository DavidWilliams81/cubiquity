#include "metadata.h"

#include "logging.h"

// The TOML parser is a fairly large single-file include which
// could potentially bloat compilation times. We only include it
// from this .cpp file and not any of our header files.
#include "toml.hpp"

#include <fmt/ostream.h>

#include <array>
#include <fstream>

const Material Material::Default = { "Default", DefaultColor };
const Material Material::EmptySpace = { "EmptySpace", std::nullopt };

const std::string root_comment =
R"(# Overview
# --------
# This text file stores volume metadata such as dimensions and materials, using 
# the TOML markup format (see https://toml.io for specs and parsers). The actual
# voxels are stored in a separate binary file with a '.bin' or '.dag' extension:
#
#   - The format of a .bin file is very simple. The voxels are stored as a raw 
#     3D array of 8-bit unsigned integers. The array dimensions are given below.
#     The array is laid out with the 'x' component changing most quickly and the
#     'z' component changing most slowly. Each voxel represents an index (0-255)
#     into the material array (also given below) describing the voxel.
#
#   - A .dag file is more complex but it can be converted to a .bin file using
#     the Cubiquity voxel engine. See https://cubiquity.net for details.
#)";

const std::string dimensions_comment =
R"(# Dimensions
# ----------
# These are dimensions of the volume (and array) along the x, y and z axis.
# Multiplying these together should match the file size of the .bin file.)";

const std::string materials_comment =
R"(# Materials
# ---------
# This array of tables contains the material properties. Each voxel is an 8-bit
# unsigned integer and so there can be up to 256 materials in the array below.
#
# Note that each material begins with *double* square brackets [[...]] as this
# is the TOML syntax for an array-of-tables.)";

void Metadata::load(const std::filesystem::path& path)
{
	// Our TOML library knows how to handle std::array but not our vec class,
	// while the compiler can implicitely convert from one to the other. Hence
	// we work with std::array while parsing, and typedefs make it tidier.
	using array3i = std::array<int, 3>;
	using array3f = std::array<float, 3>;

	toml::value toml_data = toml::parse(path);

	dimensions = toml::find<std::optional<array3i>>(toml_data, "dimensions");

	if (toml_data.contains("materials")) {
		if (toml_data.at("materials").is_array()) {
			toml::array& toml_materials_array =
				toml_data.at("materials").as_array();

			materials.resize(toml_materials_array.size());

			for (int i = 0; i < toml_materials_array.size(); i++) {
				materials[i].set_name(toml::find<std::optional<std::string>>(
					toml_data, "materials", i, "name"));

				materials[i].set_base_color(toml::find<std::optional<array3f>>(
					toml_data, "materials", i, "base_color"));
			}
		}
	}
}

void Metadata::save(const std::filesystem::path& path) const
{
	std::ofstream file(path);
	fmt::print(file, "{}\n", root_comment);

	// Write out dimensions if present
	if (dimensions.has_value())	{
		fmt::print(file, "{}\n\n", dimensions_comment);
		fmt::print(file, "dimensions = [{}, {}, {}]\n\n",
		    dimensions->x, dimensions->y, dimensions->z);
	}

	// Write out any materials
	if (materials.empty() == false) {
		fmt::print(file, "{}\n\n", materials_comment);
	}

	for (int i = 0; i < materials.size(); i++)
	{
		fmt::print(file, "[[materials]] # Material index {}\n", i);

		if (materials[i].m_name.has_value()) {
			std::string name = materials[i].m_name.value();
			//file << "name = \"" << name << "\"" << std::endl;
			fmt::print(file, "name = \"{}\"\n", name);
		}

		// Fixed precision formatting prevents scientific notation (which
		// can look messy) and ensures a decimal point is always present.
		// This ensures a TOML parser reads it back as a float (not an int).
		if (materials[i].m_base_color.has_value()) {
			vec3 base_color = materials[i].m_base_color.value();
			fmt::print(file, "base_color = [{:.4f}, {:.4f}, {:.4f}]\n",
				base_color.x, base_color.y, base_color.z);
		}

		fmt::print(file, "\n");
	}
}
