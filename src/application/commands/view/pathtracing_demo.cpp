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

uint32_t nextPointInUnitSphere = 17;
const uint pointsInUnitSphereCount = 97;
vec3 pointsInUnitSphere[pointsInUnitSphereCount] = {
{-0.383666, -0.804919, 0.0944412}, {-0.443004, -0.623236, 0.093763}, {0.596212, 0.600561, -0.405941},
{0.00732541, 0.311481, 0.595857}, {0.35747, -0.202523, 0.51548}, {0.481295, 0.486265, -0.0504826},
{-0.215546, -0.155825, 0.310956}, {-0.366899, -0.446154, 0.744858}, {0.389657, -0.749635, -0.365801},
{-0.931108, 0.327211, -0.122511}, {-0.748207, -0.236883, -0.579582}, {-0.0204712, -0.0840217, -0.108828},
{-0.0248622, 0.292626, 0.58795}, {0.615062, -0.44795, 0.411548}, {0.421408, -0.674777, 0.287922},
{-0.762005, -0.0879344, -0.00327188}, {-0.319229, 0.753515, 0.170535}, {0.502534, 0.642492, -0.48981},
{0.398153, -0.174667, 0.781806}, {-0.15367, 0.918583, 0.161913}, {0.094431, -0.683885, -0.722751},
{0.62857, -0.400337, -0.51295}, {0.232089, 0.557795, -0.0534223}, {0.661657, 0.587195, 0.170528},
{0.189007, 0.0994473, 0.465597}, {0.35964, 0.5144, -0.215359}, {0.507458, 0.123115, -0.239108},
{-0.583864, 0.135643, 0.0547429}, {-0.294475, 0.0615951, 0.185648}, {0.137647, -0.210184, -0.0612187},
{-0.325755, -0.22286, -0.675635}, {-0.37757, 0.725356, 0.0570663}, {0.203964, -0.0560864, -0.474057},
{-0.319561, 0.308158, 0.0596839}, {0.378429, 0.432201, 0.496303}, {0.0109971, 0.826675, 0.116538},
{-0.695244, 0.00638008, 0.651634}, {-0.0750516, 0.0766848, 0.0931839}, {0.708902, -0.114643, 0.208463},
{-0.273627, 0.737389, 0.359039}, {-0.831128, -0.307533, -0.200435}, {0.600137, 0.320239, -0.137172},
{-0.844886, 0.159409, 0.254769}, {0.360574, 0.706062, 0.0678662}, {0.24411, -0.122666, -0.298095},
{-0.600898, 0.026499, -0.723997}, {-0.196384, -0.235334, -0.848067}, {0.00954199, -0.520095, 0.64521},
{-0.504305, -0.0214947, 0.0881122}, {0.754728, -0.220522, -0.579396}, {-0.516617, -0.454494, -0.192176},
{-0.015116, -0.807091, -0.55316}, {-0.129509, -0.53044, 0.105762}, {0.618274, 0.298231, -0.483408},
{0.463445, -0.740307, 0.295492}, {-0.0133463, -0.0981526, -0.242999}, {0.0940177, 0.43694, -0.407358},
{-0.42439, 0.489386, 0.246871}, {-0.62209, 0.623748, 0.373551}, {-0.374984, -0.632978, -0.228579},
{-0.263031, -0.312141, 0.251237}, {0.631538, 0.560455, 0.33857}, {-0.0830062, 0.551425, -0.392772},
{-0.0264167, 0.82113, -0.128283}, {-0.0227646, -0.106432, 0.51158}, {-0.387301, -0.783876, 0.0170174},
{-0.21282, 0.0215431, 0.740373}, {0.635255, -0.278279, 0.589663}, {0.834236, 0.288636, 0.087611},
{-0.242781, -0.719712, 0.623161}, {0.100313, -0.22277, 0.24495}, {0.559379, 0.174089, -0.0467238},
{-0.584515, 0.080276, -0.397507}, {0.771996, -0.610471, -0.117553}, {-0.520995, -0.544671, 0.599306},
{-0.128603, -0.0539699, -0.377795}, {-0.139585, 0.266127, -0.630367}, {0.168764, 0.809762, 0.453309},
{0.360813, -0.777762, 0.414643}, {-0.483871, -0.675343, -0.18256}, {-0.732527, 0.189792, -0.100888},
{0.594728, 0.422432, -0.664888}, {-0.350073, -0.406648, 0.31156}, {-0.362443, -0.15962, -0.151666},
{0.563818, 0.0157166, -0.773615}, {-0.022782, -0.311344, 0.15705}, {0.211856, 0.0936115, -0.546897},
{0.0422716, -0.620189, -0.536811}, {-0.196857, -0.0222045, -0.173419}, {0.24812, -0.659695, 0.358271},
{-0.401896, -0.20897, 0.272937}, {-0.344404, -0.329286, -0.741868}, {0.359456, -0.416078, -0.726894},
{0.694219, 0.442455, -0.538471}, {-0.786476, 0.00394881, 0.307515}, {0.558103, -0.426827, 0.430074},
{0.781845, -0.289303, -0.331674}
};

// Return a small positive noise value
float PathtracingDemo::positionBasedNoise(const vec3& position)
{
	// Because the intersectionl lies exactly between two voxels a proper round to
	// nearest suffers from floating point problems. Therefore we apply a tiny offset.
	Vector3i roundedIntesectionPosition = static_cast<Vector3i>(position + vec3::filled(0.499));
	uint32 hash = murmurHash3(&roundedIntesectionPosition, sizeof(roundedIntesectionPosition));
	return (hash & 0xff) / 255.0f; // 0.0 to 1.0
}

vec3 PathtracingDemo::surfaceColour(const RayVolumeIntersection& intersection)
{
	vec3 colour = vec3({ materials()[intersection.material][0], materials()[intersection.material][1], materials()[intersection.material][2] });

	if (addNoise)
	{
		// Noise is applied multiplicatively as this avoids creating overshoots
		// and undershoots for saturated or dark surface colours respectively.
		float noise = positionBasedNoise(intersection.position);
		noise = (noise * 0.1) + 0.9; // Map 0.0 - 1.0 to range 0.9 - 1.0.
		colour *= vec3::filled(noise);
	}

	return colour;
}

const vec3& PathtracingDemo::randomPointInUnitSphere()
{
	/*static std::mt19937 gen;
	static std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

	vec3 result;
	for (int i = 0; i < 97; i++)
	{
		do
		{
			result[0] = dist(gen);
			result[1] = dist(gen);
			result[2] = dist(gen);
		} while (dot(result, result) >= 1.0f);
	}

	return result;*/

	// Note that we visit random points in pseudorandom order rather then sequentially. This is
	// because the period does need to be long, otherwise successive frames slip into synchronisation
	// (by chance) after a few hundred pixels and then never deviate again, resulting in successive
	// frames being mostly identical and hence no use for averaging together.
	nextPointInUnitSphere = mix(nextPointInUnitSphere);
	vec3 result = pointsInUnitSphere[nextPointInUnitSphere % pointsInUnitSphereCount];

	return result;
}

vec3 PathtracingDemo::gatherLighting(vec3 position, vec3 normal)
{
	const vec3 offset = normal * 0.001f;
	float intensity = 0.0f;

	if (includeSun)
	{
		float sunIntensity = 2.0f;
		vec3 sunDir(normalize(vec3({ 1.0, -2.0, 10.0 })));

		Ray3f sunShadowRay(position + offset, sunDir);
		RayVolumeIntersection sunShadowIntersection = intersectVolume(volume(), subDAGs, sunShadowRay, maxFootprint);
		if (!sunShadowIntersection)
		{
			intensity += std::max(dot(sunDir, normal), 0.0f) * sunIntensity;
		}
	}

	if (includeSky)
	{
		float skyIntensity = 1.0f;
		const vec3 skyDir = normalize(normal + randomPointInUnitSphere());
		Ray3f skyShadowRay(position + offset, skyDir);
		RayVolumeIntersection skyShadowIntersection = intersectVolume(volume(), subDAGs, skyShadowRay, maxFootprint);
		if (!skyShadowIntersection)
		{
			intensity += skyIntensity;
		}
	}

	return vec3::filled(intensity);
}

vec3 PathtracingDemo::traceSingleRayRecurse(const Ray3f& ray, uint depth)
{
	if (depth > bounces) { return vec3::filled(0); }

	vec3 pixelColour = { 0.0f, 0.0f, 0.0f };

	RayVolumeIntersection intersection = intersectVolume(volume(), subDAGs, ray, maxFootprint);
	if (intersection)
	{
		pixelColour = surfaceColour(intersection);

		vec3 directLighting = gatherLighting(intersection.position, intersection.normal);

		const vec3 reflectedDir = normalize(intersection.normal + randomPointInUnitSphere());
		const Ray3f reflectedRay(intersection.position + (intersection.normal * 0.01), reflectedDir);

		vec3 indirectLighting = traceSingleRayRecurse(reflectedRay, depth + 1);

		pixelColour *= (directLighting + indirectLighting);

	}

	return pixelColour;
}

vec3 PathtracingDemo::traceSingleRay(const Ray3f& ray, uint depth)
{
	vec3 pixelColour = { 0.0f, 0.0f, 0.0f };

	RayVolumeIntersection intersection0 = intersectVolume(volume(), subDAGs, ray, maxFootprint);
	if (intersection0)
	{
		vec3 surfCol0 = surfaceColour(intersection0);
		vec3 directLighting0 = gatherLighting(intersection0.position, intersection0.normal);

		const vec3 reflectedDir = normalize(intersection0.normal + randomPointInUnitSphere());
		const Ray3f reflectedRay(intersection0.position + (intersection0.normal * 0.01), reflectedDir);

		vec3 indirectLighting0 = {0.0f, 0.0f, 0.0f};
		RayVolumeIntersection intersection1 = intersectVolume(volume(), subDAGs, reflectedRay, maxFootprint);
		if (intersection1)
		{
			vec3 surfCol1 = surfaceColour(intersection1);
			vec3 directLighting1 = gatherLighting(intersection1.position, intersection1.normal);
			indirectLighting0 = surfCol1 * directLighting1;
		}

		pixelColour = surfCol0 * (directLighting0 + indirectLighting0);

	}

	return pixelColour;
}

void PathtracingDemo::updateResolution()
{
	int divisor = mPreviewMode ? 4 : 1;

	mImageWidth = width() / divisor;
	mImageHeight = height() / divisor;
	mImage.resize(mImageWidth * mImageHeight);

	SDL_FreeSurface(mRgbSurface);
	mRgbSurface = SDL_CreateRGBSurface(
		0, mImageWidth, mImageHeight, 24,
		0x000000ff, 0x0000ff00, 0x00ff0000, 0);

	clear();
}

void PathtracingDemo::clear()
{
	std::fill(mImage.begin(), mImage.end(), Vector3f::filled(0.0f));
	mAccumulatedFrameCount = 0;
}

void PathtracingDemo::raytrace(const Camera& camera)
{
	for (uint y = 0; y < mImageHeight; y++)
	{
		for (uint x = 0; x < mImageWidth; x++)
		{
			Ray3f ray = static_cast<Ray3f>(camera.rayFromViewportPos(x, y, mImageWidth, mImageHeight));

			Vector3f pixel = traceSingleRay(ray, 0);

			// Consider whether we need these.
			pixel = max(pixel, Vector3f::filled(0.0f));
			pixel = min(pixel, Vector3f::filled(1.0f));

			mImage[y * mImageWidth + x] += pixel;
		}
	}

	mAccumulatedFrameCount++;
}

void PathtracingDemo::onUpdate(float deltaTime)
{
	Viewer::onUpdate(deltaTime);

	// Update the pathtraced image
	Timer timer;
	raytrace(camera());
	std::cout << "Rendered frame in " << timer.elapsedTimeInMilliSeconds() << "ms" << std::endl;

	// Copy from floating point image to 24-bit RGB SDL surface
	float* src = &(mImage[0][0]);
	Uint8* dst = static_cast<Uint8*>(mRgbSurface->pixels);
	int dstActivePitch = mRgbSurface->w * mRgbSurface->format->BytesPerPixel;
	const float scaleFactor = 255.0f / mAccumulatedFrameCount;
	for (uint y = 0; y < mImageHeight; y++)
	{
		Uint8* dstBegin = dst + (y * mRgbSurface->pitch);
		Uint8* dstEnd = dstBegin + dstActivePitch;
		for (Uint8* dst = dstBegin; dst < dstEnd; ++src, ++dst)
		{
			*dst = *src * scaleFactor;
		}
	}

	// Blit to display with scaling
	SDL_BlitScaled(mRgbSurface, NULL, surface(), NULL);
}

void PathtracingDemo::setPreviewMode(bool previewMode)
{
	mPreviewMode = previewMode;
	updateResolution();
}

void PathtracingDemo::onKeyDown(const SDL_KeyboardEvent& event)
{
	Viewer::onKeyDown(event);

	if (event.keysym.sym == SDLK_F1)
	{
		if (bounces > 0) { bounces--; }
		clear();
		std::cout << "Bounces = " << bounces << std::endl;
	}
	if (event.keysym.sym == SDLK_F2)
	{
		if (bounces < 5) { bounces++; }
		clear();
		std::cout << "Bounces = " << bounces << std::endl;
	}
	if (event.keysym.sym == SDLK_F3)
	{
		includeSun = !includeSun;
		clear();
		std::cout << "includeSun = " << includeSun << std::endl;
	}
	if (event.keysym.sym == SDLK_F4)
	{
		includeSky = !includeSky;
		clear();
		std::cout << "includeSky = " << includeSky << std::endl;
	}
	if (event.keysym.sym == SDLK_F5)
	{
		addNoise = !addNoise;
		
		std::cout << "addNoise = " << addNoise << std::endl;
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
		clear();
	}
}

void PathtracingDemo::onWindowSizeChanged(int width, int height)
{
	updateResolution();
}

void PathtracingDemo::onCameraModified()
{
	clear();
}

void PathtracingDemo::onVolumeModified()
{
	subDAGs = findSubDAGs(
		Internals::getNodes(volume()).nodes(), getRootNodeIndex(volume()));

	clear();
}
