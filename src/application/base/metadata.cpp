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

	// We don't read the header - it's just a comment string so not someting
	// we can access easily. But we don't need to pass it trhrough anyway.
	header = "";

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

void Metadata::save(std::ostream& os) const
{
	// Write out the header if present
	if (!header.empty()) {
		fmt::print(os, "{}\n", header);
	}

	// Write out dimensions if present
	if (dimensions.has_value())	{
		fmt::print(os, "{}\n\n", dimensions_comment);
		fmt::print(os, "dimensions = [{}, {}, {}]\n\n",
		    dimensions->x, dimensions->y, dimensions->z);
	}

	// Write out any materials
	if (materials.empty() == false) {
		fmt::print(os, "{}\n\n", materials_comment);
	}

	for (int i = 0; i < materials.size(); i++)
	{
		fmt::print(os, "[[materials]] # Material index {}\n", i);

		if (materials[i].m_name.has_value()) {
			std::string name = materials[i].m_name.value();
			//os << "name = \"" << name << "\"" << std::endl;
			fmt::print(os, "name = \"{}\"\n", name);
		}

		// Fixed precision formatting prevents scientific notation (which
		// can look messy) and ensures a decimal point is always present.
		// This ensures a TOML parser reads it back as a float (not an int).
		if (materials[i].m_base_color.has_value()) {
			vec3 base_color = materials[i].m_base_color.value();
			fmt::print(os, "base_color = [{:.4f}, {:.4f}, {:.4f}]\n",
				base_color.x, base_color.y, base_color.z);
		}

		fmt::print(os, "\n");
	}
}
