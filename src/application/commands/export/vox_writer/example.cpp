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
	std::cout << "Finshed writing " << filename << std::endl;
	return 0;
}
