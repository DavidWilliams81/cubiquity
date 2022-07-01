#ifndef CUBIQUITY_INSTANCING_DEMO_H
#define CUBIQUITY_INSTANCING_DEMO_H

#include "demo.h"
#include "camera.h"
#include "instance_list.h"
#include "volume_renderer.h"

#include "framework/program.h"

#include "storage.h"
#include "rendering.h"

class InstancingDemo : public Demo
{
public:

	InstancingDemo(const std::string& filename)
		: Demo(filename, WindowType::OpenGL) {}

	VolumeRenderer* mVolumeRenderer;

	void onInitialise() override;
	void onUpdate(float deltaTime) override;
	void onShutdown() override;

	void onKeyDown(const SDL_KeyboardEvent& event);
};

#endif // CUBIQUITY_INSTANCING_DEMO_H
