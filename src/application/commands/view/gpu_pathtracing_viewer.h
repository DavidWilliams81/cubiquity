#ifndef CUBIQUITY_GPU_PATHTRACING_VIEWER_H
#define CUBIQUITY_GPU_PATHTRACING_VIEWER_H

#include "opengl_viewer.h"

#include "visibility.h"
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

	void onMouseButtonDown(const SDL_MouseButtonEvent& event) override;
	void onMouseButtonUp(const SDL_MouseButtonEvent& event) override;

	void onCameraModified();

	void clear();

	GLuint previewProgram;
	GLuint progressiveProgram;
	GLuint copyProgram;
	GLuint hBlurProgram;
	GLuint vBlurProgram;

	GLuint mMaterialsTexture;

	unsigned int accFramebuffer;
	unsigned int accTexture;

	unsigned int blurFramebuffer;
	unsigned int blurTexture;

	Cubiquity::Timer mTimer;

	unsigned int frameId = 0;
	unsigned int staticFrameCount = 0;
};

#endif // CUBIQUITY_GPU_PATHTRACING_VIEWER_H
