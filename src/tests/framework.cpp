#include "framework.h"

#include "utility.h"

using namespace Cubiquity;
using namespace Cubiquity::Internals;

void saveVisibilityMaskAsImage(VisibilityMask& visMask, const std::string& filename)
{
	Image image(visMask.width(), visMask.height());
	for (uint32_t y = 0; y < visMask.height(); y++)
	{
		for (uint32_t x = 0; x < visMask.width(); x++)
		{
			image.setPixel(x, y, visMask.testPixel(x, y) ? 255 : 0);
		}
	}

	image.save(filename);
}
