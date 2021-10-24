/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#ifndef CUBIQUITY_PATHTRACING_DEMO_H
#define CUBIQUITY_PATHTRACING_DEMO_H

#include "demo.h"

// Cubiquity
#include "geometry.h"
#include "rendering.h"
#include "utility.h"
#include "storage.h"

// Utility
#include "camera.h"

// SDL
#define SDL_MAIN_HANDLED
#include <SDL.h>

// Standard library
#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>

using namespace Cubiquity;

class Pathtracer
{
public:

	Pathtracer();
	void setResolution(uint width, uint height);
	void clear();
	void raytrace(const Volume& volume, const MaterialSet& materials, const Camera& camera, uint timeoutMs = 0); // Zero for no timeout

	uint mWidth;
	uint mHeight;
	Vector3u8* renderSurface = nullptr;

	uint bounces = 1;
	uint samples = 1;
	bool includeSun = true;
	bool includeSky = true;
	bool addNoise = true;

	// Iterates over pixels in a pseudorandom order (used for progressive rendering).
	ShuffledSequence gPixelSequence;
	uint pixelCount;
};

class PathtracingDemo : public Demo
{
public:

	PathtracingDemo(const std::string& filename)
		: Demo(filename, WindowType::Software) {}

	void setPreviewMode(bool previewMode);

protected:

	void onUpdate(float deltaTime) override;

	void onKeyDown(const SDL_KeyboardEvent& event) override;

	void onMouseButtonDown(const SDL_MouseButtonEvent& event) override;
	void onMouseButtonUp(const SDL_MouseButtonEvent& event) override;
	void onMouseMotion(const SDL_MouseMotionEvent& event) override;

	void onWindowSizeChanged(int width, int height) override;

	void onCameraModified() override;
	void onVolumeModified() override;

	Pathtracer mPathtracer;
	bool mPreviewMode = false;
};

#endif // CUBIQUITY_PATHTRACING_DEMO_H
