#include "volume_vox_writer.h"

#include "base/logging.h"

#include "base.h"

#include "cubiquity.h"

using Cubiquity::Volume;

volume_vox_writer::volume_vox_writer(Volume& vol, const Metadata& metadata)
	:m_vol(vol), m_metadata(metadata)
{
	// FIXME - Restore this test? Check opacity is zero?
	// MagicaVoxel assumes that palette index 0 is empty space, so we
	// can only write a valid .vox file if our volume does the same.
	//if (metadata.name(0) != Metadata::EmptySpace.name) {
	//	log_warning("Material '{}' found in slot zero (it should be empty space "
	//	            "and will be treated as such)", metadata.name(0));
	//}

	for (int i = 1; i < 256; i++)
	{
		if (i < metadata.material_count())
		{
			vec3 base_color = metadata.find_material_base_color(i);

			float gamma = 1.0f / 2.2f;
			base_color[0] = pow(base_color[0], gamma);
			base_color[1] = pow(base_color[1], gamma);
			base_color[2] = pow(base_color[2], gamma);

			u8 r = std::clamp(std::lround(base_color[0] * 255.0f), 0L, 255L);
			u8 g = std::clamp(std::lround(base_color[1] * 255.0f), 0L, 255L);
			u8 b = std::clamp(std::lround(base_color[2] * 255.0f), 0L, 255L);

			set_palette(i, { r, g, b, 255 });
		}
	}
}

vox_writer::box volume_vox_writer::bounds()
{
	ivec3 lower_bound = m_metadata.find_lower_bound();
	ivec3 upper_bound = m_metadata.find_upper_bound();

	return { { lower_bound.x, lower_bound.y, lower_bound.z }, 
	         { upper_bound.x, upper_bound.y, upper_bound.z } };
}

u8 volume_vox_writer::voxel(const vec3i& position)
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
