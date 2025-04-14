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

#include "viewer.h"

// Cubiquity
#include "geometry.h"
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

using namespace Cubiquity;

typedef vec3f vec3;
typedef vec3i ivec3;
typedef vec4i ivec4;

class PathtracingDemo : public Viewer
{
public:

	PathtracingDemo(const std::string& filename)
		: Viewer(filename, WindowType::Software) {}

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

private:
	void updateResolution();
	void clear();
	void raytrace(const Camera& camera);

	float positionBasedNoise(const vec3& position);
	vec3 surfaceColour(const RayVolumeIntersection& intersection);
	vec3 randomPointInUnitSphere();
	vec3 gatherLighting(vec3 position, vec3 normal);
	vec3 traceSingleRayRecurse(const Ray3f& ray, uint depth);
	vec3 traceSingleRay(const Ray3f& ray, uint depth);

	// Floating point target for path tracing
	uint mImageWidth;
	uint mImageHeight;
	std::vector<vec3f> mImage;

	// Intermediate SDL target for blitting with scaling.
	SDL_Surface* mRgbSurface = nullptr;

	// Pathtracing properties
	uint bounces = 1;
	bool addNoise = true;
	bool includeSun = true;
	bool includeSky = true;
	float maxFootprint = 0.0035; // For LOD

	// General properties
	uint mAccumulatedFrameCount = 0;
	bool mPreviewMode = false;

	SubDAGArray subDAGs;
};

#endif // CUBIQUITY_PATHTRACING_DEMO_H
