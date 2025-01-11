#include "volume_vox_writer.h"

#include "base/logging.h"

#include "base.h"

#include "cubiquity.h"

using namespace Cubiquity;

volume_vox_writer::volume_vox_writer(Volume& vol, const Metadata& metadata)
	:m_vol(vol), m_metadata(metadata)
{
	// MagicaVoxel assumes that palette index 0 is empty space, so we
	// can only write a valid .vox file if our volume does the same.
	if (metadata.materials[0].name != Metadata::EmptySpace.name) {
		log_warning("Material '{}' found in slot zero (it should be empty space "
		            "and will be treated as such)", metadata.materials[0].name);
	}

	for (int i = 1; i < 256; i++)
	{
		if (i < metadata.materials.size())
		{
			Col base_color = metadata.materials.at(i).base_color;

			float gamma = 1.0f / 2.2f;
			base_color[0] = pow(base_color[0], gamma);
			base_color[1] = pow(base_color[1], gamma);
			base_color[2] = pow(base_color[2], gamma);

			uint8 r = std::clamp(std::lround(base_color[0] * 255.0f), 0L, 255L);
			uint8 g = std::clamp(std::lround(base_color[1] * 255.0f), 0L, 255L);
			uint8 b = std::clamp(std::lround(base_color[2] * 255.0f), 0L, 255L);

			set_palette(i, { r, g, b, 255 });
		}
	}
}

vox_writer::box volume_vox_writer::bounds()
{
	uint8 outside_material;
	int32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&m_vol, &outside_material,
		                      &lower_x, &lower_y, &lower_z,
		                      &upper_x, &upper_y, &upper_z);

	// Including a border lets the user see the objects have not been
	// clipped, and if the outside  material is solid then it makes
	// sense to include some of it as part of the scene anyway.
	const int border = 1;
	lower_x -= border; lower_y -= border; lower_z -= border;
	upper_x += border; upper_y += border; upper_z += border;

	return{ { lower_x, lower_y, lower_z }, { upper_x, upper_y, upper_z } };
}

uint8_t volume_vox_writer::voxel(const vec3i& position)
{
	return m_vol.voxel(position.x, position.y, position.z);
}

void volume_vox_writer::on_progress(int done, int total)
{
	if(done == 0) {
		log_info_no_newline("Exporting as .vox file");
	}
	log_info_no_newline(".");
	// FIXME - We would benefit from log_info_if() here.
	if(done == total) {
		log_info(" Done");
	}
}
