#ifndef CUBIQUITY_GPU_PATHTRACING_VIEWER_H
#define CUBIQUITY_GPU_PATHTRACING_VIEWER_H

#include "opengl_viewer.h"

#include "rendering.h"
#include "utility.h"


class GPUPathtracingViewer : public OpenGLViewer
{
public:

	GPUPathtracingViewer(const std::string& filename)
		: OpenGLViewer(filename) {}

	void onInitialise() override;
	void onUpdate(float deltaTime) override;
	void onShutdown() override;

	void onKeyDown(const SDL_KeyboardEvent& event);

	void onCameraModified();

	void clear();

	GLuint mainProgram;
	GLuint screenQuadProgram;

	GLuint mMaterialsTexture;

	unsigned int framebuffer;
	unsigned int textureColorbuffer;

	Cubiquity::Timer mTimer;

	unsigned int randSeed = 42;
};

#endif // CUBIQUITY_GPU_PATHTRACING_VIEWER_H
