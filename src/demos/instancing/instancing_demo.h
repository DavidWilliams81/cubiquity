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

	InstancingDemo(int argc, char** argv)
		: Demo(argc, argv, WindowType::OpenGL) {}

	VolumeRenderer* mVolumeRenderer;

	void onInitialise() override;
	void onUpdate(float deltaTime) override;
	void onShutdown() override;
};

#endif // CUBIQUITY_INSTANCING_DEMO_H
