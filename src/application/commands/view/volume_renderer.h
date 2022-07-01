#ifndef VOLUME_RENDERER_H_1F2E659B
#define VOLUME_RENDERER_H_1F2E659B

#include "camera.h"
#include "instance_list.h"
#include "framework/program.h"

#include "base/materials.h"

#include "storage.h"
#include "rendering.h"

enum class GlyphType
{
	Cube,
	Point
};

class VolumeRenderer
{
public:
	VolumeRenderer(const MaterialSet& materials);
	virtual ~VolumeRenderer();

	virtual void render(const Camera& camera);

	Program program;

	GlyphType mGlyphType;

	const Cubiquity::Volume* mVolume = nullptr;
	Cubiquity::CameraData cameraData;

	GLuint modelMatrixID;
	GLuint viewMatrixID;
	GLuint projMatrixID;
	GLuint cameraPosID;

	GLuint Texture;
	GLuint TextureID;

	InstanceList instanceList;

	Cubiquity::VisibilityCalculator* mVisibilityCalculator = nullptr;

	Cubiquity::Glyph* mGlyphs;
	bool mDoGlyphUpdates = true;
};

#endif // VOLUME_RENDERER_H_1F2E659B
