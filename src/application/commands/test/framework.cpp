#include "framework.h"

#include "utility.h"

#include "stb_image_write.h"

using namespace Cubiquity;
using namespace Cubiquity::Internals;

void saveVisibilityMaskAsImage(VisibilityMask& visMask, const std::string& filename)
{
	std::vector<uint8> imageData;
	for (uint32_t y = 0; y < visMask.height(); y++)
	{
		for (uint32_t x = 0; x < visMask.width(); x++)
		{
			imageData.push_back(visMask.testPixel(x, y) ? 255 : 0);
		}
	}

	int result = stbi_write_png(filename.c_str(), visMask.width(), visMask.height(), 1, imageData.data(), visMask.width());
	if (result == 0)
	{
		log_error("Failed to write PNG image");
	}
}
