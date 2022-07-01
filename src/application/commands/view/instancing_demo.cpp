#include "instancing_demo.h"

#include "framework/program.h"


#include "volume_renderer.h"

#include "utility.h"

#include <iostream>

using namespace Cubiquity;

void InstancingDemo::onInitialise()
{
	Demo::onInitialise();

	mVolumeRenderer = new VolumeRenderer(materials());

	mVolumeRenderer->mVolume = &(volume());

	// Black background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	// Depth test disabled for points as it interferes with blending when the glyphs are scaled so that they overlap.
	mVolumeRenderer->mGlyphType == GlyphType::Point ? glDisable(GL_DEPTH_TEST) : glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
}

void InstancingDemo::onUpdate(float deltaTime)
{
	Demo::onUpdate(deltaTime);

	mVolumeRenderer->render(camera());
}

void InstancingDemo::onShutdown()
{
	Demo::onShutdown();

	delete mVolumeRenderer;
}

void InstancingDemo::onKeyDown(const SDL_KeyboardEvent & event)
{
	Demo::onKeyUp(event);

	if (event.keysym.sym == SDLK_SPACE)
	{
		mVolumeRenderer->mDoGlyphUpdates = !(mVolumeRenderer->mDoGlyphUpdates);
	}
}
