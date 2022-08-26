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

#include "pathtracing_demo.h"

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

Vector3d randomPointInUnitSphere()
{
	static std::mt19937 gen;
	static std::uniform_real_distribution<double> dist(-1.0, 1.0);

	Vector3d result;
	do
	{
		result[0] = dist(gen);
		result[1] = dist(gen);
		result[2] = dist(gen);
	} while (dot(result, result) >= 1.0);

	return result;
}

// Return a small positive noise value
double positionBasedNoise(const Vector3d& position)
{
	// Because the intersectionl lies exactly between two voxels a proper round to
	// nearest suffers from floating point problems. Therefore we apply a tiny offset.
	Vector3i roundedIntesectionPosition = static_cast<Vector3i>(position + Vector3d::filled(0.499));
	uint32 hash = murmurHash3(&roundedIntesectionPosition, sizeof(roundedIntesectionPosition));
	return (hash & 0xff) / 255.0f; // 0.0 to 1.0
}

Vector3d gatherLighting(Vector3d position, Vector3d normal, const Volume& volume,
	bool includeSun = true, bool includeSky = true)
{
	const Vector3d offset = normal * 0.001;
	double intensity = 0.0;

	if (includeSun)
	{
		double sunIntensity = 2.0;
		Vector3d sunDir(normalize(Vector3d({ 1.0, -2.0, 10.0 })));

		Ray3d sunShadowRay(position + offset, sunDir);
		RayVolumeIntersection sunShadowIntersection = ray_parameter(volume, sunShadowRay);
		if (!sunShadowIntersection)
		{
			intensity += std::max(dot(sunDir, normal), 0.0) * sunIntensity;
		}
	}

	if (includeSky)
	{
		double skyIntensity = 1.0;
		const Vector3d skyDir = normalize(normal + randomPointInUnitSphere());
		Ray3d skyShadowRay(position + offset, skyDir);
		RayVolumeIntersection skyShadowIntersection = ray_parameter(volume, skyShadowRay);
		if (!skyShadowIntersection)
		{
			intensity += skyIntensity;
		}
	}

	return Vector3d::filled(intensity);
}

Vector3d traceSingleRay(const Ray3d& ray, const Volume& volume, const MaterialSet& materials, uint bounces, bool includeSun, bool includeSky, bool addNoise, uint depth)
{
	if (depth > bounces) { return Vector3d::filled(0); }

	Vector3d pixelColour = { 0.0f, 0.0f, 0.0f };

	RayVolumeIntersection intersection = ray_parameter(volume, ray);
	if (intersection)
	{
		pixelColour = Vector3d({ materials[intersection.material][0], materials[intersection.material][1], materials[intersection.material][2] });

		if (addNoise)
		{
			// Noise is applied multiplicatively as this avoids creating overshoots
			// and undershoots for saturated or dark surface colours respectively.
			float noise = positionBasedNoise(intersection.position);
			noise = (noise * 0.1) + 0.9; // Map 0.0 - 1.0 to range 0.9 - 1.0.
			pixelColour *= Vector3d::filled(noise);
		}

		if (includeSun || includeSky)
		{
			Vector3d directLighting = gatherLighting(intersection.position, intersection.normal,
				volume, includeSun, includeSky);

			const Vector3d reflectedDir = normalize(intersection.normal + randomPointInUnitSphere());
			const Ray3d reflectedRay(intersection.position + (intersection.normal * 0.01), reflectedDir);

			Vector3d indirectLighting = traceSingleRay(reflectedRay, volume, materials, bounces,
				includeSun, includeSky, addNoise, depth + 1);

			pixelColour *= (directLighting + indirectLighting);
		}
	}

	return pixelColour;
}

Pathtracer::Pathtracer()
	:gPixelSequence(12345)
{

}

void Pathtracer::setResolution(uint width, uint height)
{
	mWidth = width;
	mHeight = height;

	delete[] renderSurface;
	renderSurface = new Vector3u8[width * height];
	gPixelSequence = ShuffledSequence(width * height);

	clear();
}

void Pathtracer::clear()
{
	memset(renderSurface, 0, mWidth * mHeight * sizeof(renderSurface[0]));
	pixelCount = 0;
}

void Pathtracer::raytrace(const Volume& volume, const MaterialSet& materials, const Camera& camera, uint timeoutMs)
{
	Timer timer;
	while(pixelCount < mWidth * mHeight)
	{
		uint32 index = gPixelSequence.state();
		gPixelSequence.next();

		uint y = index / mWidth;
		uint x = index % mWidth;


		//Uint8* start = ((Uint8*)surface->pixels);
		//start += surface->pitch * y;
		//Uint32* pixel = ((Uint32*)start) + x;

		Ray3d ray = camera.rayFromViewportPos(x, y, mWidth, mHeight);

		Vector3f pixel = { 0, 0, 0 };
		for (uint sample = 0; sample < samples; sample++)
		{
			pixel += static_cast<Vector3f>(traceSingleRay(ray, volume, materials, bounces, includeSun, includeSky, addNoise, 0));
		}
		pixel /= samples;

		Vector3u8* start = ((Vector3u8*)renderSurface);
		start += mWidth * y;
		Vector3u8* pixelColour = ((Vector3u8*)start) + x;

		// Clamp to valid range for conversion.
		pixel = max(pixel, Vector3f::filled(0.0f));
		pixel = min(pixel, Vector3f::filled(1.0f));

		uint8 red = static_cast<uint8>(pixel.x() * 255.0f);
		uint8 green = static_cast<uint8>(pixel.y() * 255.0f);
		uint8 blue = static_cast<uint8>(pixel.z() * 255.0f);

		*pixelColour = Vector3u8({ red, green, blue });

		pixelCount++;

		if (timeoutMs && timer.elapsedTimeInMilliSeconds() > static_cast<float>(timeoutMs))
		{
			return;
		}
	}
}

void PathtracingDemo::onUpdate(float deltaTime)
{
	Viewer::onUpdate(deltaTime);

	// No timeout needed in preview mode because it is fast anyway, and having
	// it appear immeditely looks nicer than seeing the progressive loading.
	mPathtracer.raytrace(volume(), materials(), camera(), mPreviewMode ? 0 : 50);

	SDL_Surface* renderedSurface = SDL_CreateRGBSurfaceFrom(mPathtracer.renderSurface,
		mPathtracer.mWidth, mPathtracer.mHeight, 24,
		mPathtracer.mWidth * sizeof(mPathtracer.renderSurface[0]),
		0x000000ff,	0x0000ff00,	0x00ff0000,	0);

	SDL_BlitScaled(renderedSurface, NULL, surface(), NULL);

	SDL_FreeSurface(renderedSurface);
}

void PathtracingDemo::setPreviewMode(bool previewMode)
{
	mPreviewMode = previewMode;
	if (previewMode)
	{
		mPathtracer.samples = 1;
		mPathtracer.bounces = 0;
		mPathtracer.addNoise = false;
		mPathtracer.includeSun = false;
		mPathtracer.includeSky = false;
		mPathtracer.setResolution(surface()->w / 4, surface()->h / 4);
	}
	else
	{
		mPathtracer.samples = 5;
		mPathtracer.bounces = 1;
		mPathtracer.addNoise = true;
		mPathtracer.includeSun = true;
		mPathtracer.includeSky = true;
		mPathtracer.setResolution(surface()->w, surface()->h);
	}
}

void PathtracingDemo::onKeyDown(const SDL_KeyboardEvent& event)
{
	Viewer::onKeyDown(event);

	if (event.keysym.sym == SDLK_F1)
	{
		if (mPathtracer.bounces > 0) { mPathtracer.bounces--; }
		mPathtracer.clear();
		std::cout << "Bounces = " << mPathtracer.bounces << std::endl;
	}
	if (event.keysym.sym == SDLK_F2)
	{
		if (mPathtracer.bounces < 5) { mPathtracer.bounces++; }
		mPathtracer.clear();
		std::cout << "Bounces = " << mPathtracer.bounces << std::endl;
	}
	if (event.keysym.sym == SDLK_F3)
	{
		mPathtracer.includeSun = !mPathtracer.includeSun;
		mPathtracer.clear();
		std::cout << "includeSun = " << mPathtracer.includeSun << std::endl;
	}
	if (event.keysym.sym == SDLK_F4)
	{
		mPathtracer.includeSky = !mPathtracer.includeSky;
		mPathtracer.clear();
		std::cout << "includeSky = " << mPathtracer.includeSky << std::endl;
	}
	if (event.keysym.sym == SDLK_F5)
	{
		mPathtracer.addNoise = !mPathtracer.addNoise;
		
		std::cout << "addNoise = " << mPathtracer.addNoise << std::endl;
	}
}

void PathtracingDemo::onMouseButtonDown(const SDL_MouseButtonEvent& event)
{
	Viewer::onMouseButtonDown(event);
	setPreviewMode(true);
}

void PathtracingDemo::onMouseButtonUp(const SDL_MouseButtonEvent& event)
{
	Viewer::onMouseButtonUp(event);
	setPreviewMode(false);
}

void PathtracingDemo::onMouseMotion(const SDL_MouseMotionEvent& event)
{
	Viewer::onMouseMotion(event);

	if (mouseButtonState(SDL_BUTTON_RIGHT) == MouseButtonState::Down)
	{
		mPathtracer.clear();
	}
}

void PathtracingDemo::onWindowSizeChanged(int width, int height)
{
	mPathtracer.setResolution(width, height);
}

void PathtracingDemo::onCameraModified()
{
	mPathtracer.clear();
}

void PathtracingDemo::onVolumeModified()
{
	mPathtracer.clear();
}
