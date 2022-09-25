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


	GLuint screenQuadProgram;
	GLuint mMaterialsTexture;

	Cubiquity::Timer mTimer;
};

#endif // CUBIQUITY_GPU_PATHTRACING_VIEWER_H
