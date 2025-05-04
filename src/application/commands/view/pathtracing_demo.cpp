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

#include "base/logging.h"

// Cubiquity
#include "extraction.h"
#include "raytracing.h"
#include "utility.h"
#include "storage.h"

// Utility
#include "camera.h"

// Standard library
#include <cstdlib>
#include <cmath>
#include <functional>
#include <random>

uint32_t nextPointInUnitSphere = 17;

// Return a small positive noise value
float PathtracingDemo::positionBasedNoise(const vec3& position)
{
	// Because the intersectionl lies exactly between two voxels a proper round to
	// nearest suffers from floating point problems. Therefore we apply a tiny offset.
	ivec3 roundedIntesectionPosition = static_cast<ivec3>(position + vec3(0.499));
	uint32_t hash = Cubiquity::murmurHash3(&roundedIntesectionPosition, sizeof(roundedIntesectionPosition));
	return (hash & 0xff) / 255.0f; // 0.0 to 1.0
}

vec3 PathtracingDemo::surfaceColour(const Cubiquity::RayVolumeIntersection& intersection)
{
	vec3 colour = vec3({ colours()[intersection.material][0], colours()[intersection.material][1], colours()[intersection.material][2] });

	if (addNoise)
	{
		// Noise is applied multiplicatively as this avoids creating overshoots
		// and undershoots for saturated or dark surface colours respectively.
		vec3 position(intersection.position.x, intersection.position.y, intersection.position.z);
		float noise = positionBasedNoise(position);
		noise = (noise * 0.1) + 0.9; // Map 0.0 - 1.0 to range 0.9 - 1.0.
		colour *= vec3(noise);
	}

	return colour;
}

vec3 PathtracingDemo::randomPointInUnitSphere()
{
	vec3 result;
	do
	{
		nextPointInUnitSphere = Cubiquity::mixBits(nextPointInUnitSphere);
		result[0] = nextPointInUnitSphere & 0x3FF;
		result[1] = (nextPointInUnitSphere >> 10) & 0x3FF;
		result[2] = (nextPointInUnitSphere >> 20) & 0x3FF;

		vec3 offset = { 511.5f, 511.5f, 511.5f };
		result = result - offset;
		result /= 511.5f;

	} while (dot(result, result) >= 1.0f);

	return result;
}

vec3 PathtracingDemo::gatherLighting(vec3 position, vec3 normal)
{
	const vec3 offset = normal * 0.001f;
	vec3 intensity = vec3(0.0f);

	if (includeSun)
	{
		const vec3 sunColour = vec3({ 0.1, 0.1, 0.1 });
		vec3 sunDir(normalize(vec3({ 1.0, -2.0, 10.0 })));

		Ray3f sunShadowRay(position + offset, sunDir);
		Cubiquity::RayVolumeIntersection sunShadowIntersection = intersectVolume(volume(), subDAGs,
			sunShadowRay.mOrigin.x, sunShadowRay.mOrigin.y, sunShadowRay.mOrigin.z,
			sunShadowRay.mDir.x, sunShadowRay.mDir.y, sunShadowRay.mDir.z,
			false, maxFootprint);
		if (sunShadowIntersection.hit == false)
		{
			intensity += sunColour * std::max(dot(sunDir, normal), 0.0f);
		}
	}

	if (includeSky)
	{
		const vec3 skyColour = vec3({ 1.5f, 1.5f, 1.5f });
		const vec3 skyDir = normalize(normal + randomPointInUnitSphere());
		Ray3f skyShadowRay(position + offset, skyDir);
		Cubiquity::RayVolumeIntersection skyShadowIntersection = intersectVolume(volume(), subDAGs,
			skyShadowRay.mOrigin.x, skyShadowRay.mOrigin.y, skyShadowRay.mOrigin.z,
			skyShadowRay.mDir.x, skyShadowRay.mDir.y, skyShadowRay.mDir.z,
			false, maxFootprint);
		if (skyShadowIntersection.hit == false)
		{
			intensity += skyColour;
		}
	}

	return intensity;
}

vec3 PathtracingDemo::traceSingleRayRecurse(const Ray3f& ray, uint depth)
{
	if (depth > bounces) { return vec3(0, 0, 0); }

	vec3 pixelColour = { 0.8f, 0.8f, 1.0f }; // Light blue background

	Cubiquity::RayVolumeIntersection intersection = intersectVolume(volume(), subDAGs,
		ray.mOrigin.x, ray.mOrigin.y, ray.mOrigin.z,
		ray.mDir.x, ray.mDir.y, ray.mDir.z,
		true, maxFootprint);
	if (intersection.hit)
	{
		pixelColour = surfaceColour(intersection);

		vec3 position(intersection.position.x, intersection.position.y, intersection.position.z);
		vec3 normal(intersection.normal.x, intersection.normal.y, intersection.normal.z);
		vec3 directLighting = gatherLighting(position, normal);

		const vec3 reflectedDir = normalize(normal + randomPointInUnitSphere());
		const Ray3f reflectedRay(position + (normal * 0.01f), reflectedDir);

		vec3 indirectLighting = traceSingleRayRecurse(reflectedRay, depth + 1);

		pixelColour *= (directLighting + indirectLighting);

	}

	return pixelColour;
}

vec3 PathtracingDemo::traceSingleRay(const Ray3f& ray, uint depth)
{
	vec3 pixelColour = { 0.8f, 0.8f, 1.0f }; // Light blue background

	Cubiquity::RayVolumeIntersection intersection0 = intersectVolume(volume(), subDAGs,
		ray.mOrigin.x, ray.mOrigin.y, ray.mOrigin.z,
		ray.mDir.x, ray.mDir.y, ray.mDir.z,
		true, maxFootprint);
	if (intersection0.hit)
	{
		vec3 surfCol0 = surfaceColour(intersection0);
		vec3 position0(intersection0.position.x, intersection0.position.y, intersection0.position.z);
		vec3 normal0(intersection0.normal.x, intersection0.normal.y, intersection0.normal.z);
		vec3 directLighting0 = gatherLighting(position0, normal0);

		const vec3 reflectedDir = normalize(normal0 + randomPointInUnitSphere());
		const Ray3f reflectedRay(position0 + (normal0 * 0.01f), reflectedDir);

		vec3 indirectLighting0 = {0.0f, 0.0f, 0.0f};
		Cubiquity::RayVolumeIntersection intersection1 = intersectVolume(volume(), subDAGs,
			reflectedRay.mOrigin.x, reflectedRay.mOrigin.y, reflectedRay.mOrigin.z,
			reflectedRay.mDir.x, reflectedRay.mDir.y, reflectedRay.mDir.z,
			true, maxFootprint);
		if (intersection1.hit)
		{
			vec3 surfCol1 = surfaceColour(intersection1);
			vec3 position1(intersection1.position.x, intersection1.position.y, intersection1.position.z);
			vec3 normal1(intersection1.normal.x, intersection1.normal.y, intersection1.normal.z);
			vec3 directLighting1 = gatherLighting(position1, normal1);
			indirectLighting0 = surfCol1 * directLighting1;
		}

		pixelColour = surfCol0 * (directLighting0 + indirectLighting0);

		float gamma = 1.0 / 2.2;
		pixelColour = pow(pixelColour, gamma);

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
	std::fill(mImage.begin(), mImage.end(), vec3(0.0f));
	mAccumulatedFrameCount = 0;
}

void PathtracingDemo::raytrace(const Camera& camera)
{
	for (uint y = 0; y < mImageHeight; y++)
	{
		for (uint x = 0; x < mImageWidth; x++)
		{
			Ray3f ray = static_cast<Ray3f>(camera.rayFromViewportPos(x, y, mImageWidth, mImageHeight));

			vec3 pixel = traceSingleRay(ray, 0);

			mImage[y * mImageWidth + x] += pixel;
		}
	}

	mAccumulatedFrameCount++;
}

void PathtracingDemo::onUpdate(float deltaTime)
{
	Viewer::onUpdate(deltaTime);

	// Update the pathtraced image
	Cubiquity::Timer timer;
	raytrace(camera());
	log_info("Rendered frame in {}ms", timer.elapsedTimeInMilliSeconds());

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
			float fDst = *src * scaleFactor;
			fDst = std::max(fDst, 0.0f);
			fDst = std::min(fDst, 255.0f);
			*dst = static_cast<Uint8>(fDst);
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
		log_info("Bounces = {}", bounces);
	}
	if (event.keysym.sym == SDLK_F2)
	{
		if (bounces < 5) { bounces++; }
		clear();
		log_info("Bounces = {}", bounces);
	}
	if (event.keysym.sym == SDLK_F3)
	{
		includeSun = !includeSun;
		clear();
		log_info("includeSun = {}", includeSun);
	}
	if (event.keysym.sym == SDLK_F4)
	{
		includeSky = !includeSky;
		clear();
		log_info("includeSky = {}", includeSky);
	}
	if (event.keysym.sym == SDLK_F5)
	{
		addNoise = !addNoise;
		log_info("addNoise = {}", addNoise);
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
	subDAGs = Cubiquity::findSubDAGs(
		Cubiquity::Internals::getNodes(volume()).nodes(), Cubiquity::getRootNodeIndex(volume()));

	clear();
}
