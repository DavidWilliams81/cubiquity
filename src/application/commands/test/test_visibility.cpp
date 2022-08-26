#include "test_visibility.h"

#include "framework.h"

#include "utility.h"
#include "geometry.h"
#include "fractal_noise.h"
#include "storage.h"
#include "rendering.h"

#include <functional>
#include <iostream>
#include <random>

using namespace Cubiquity;
using namespace std;

template <typename Function>
void fillVolume(Volume* volume, const Box3i& bounds, Function function)
{

	for (uint32_t z = bounds.lower().z(); z <= bounds.upper().z(); z++)
	{
		std::cout << "Generating slice " << z << std::endl;

		for (uint32_t y = bounds.lower().y(); y <= bounds.upper().y(); y++)
		{
			for (uint32_t x = bounds.lower().x(); x <= bounds.upper().x(); x++)
			{
				volume->setVoxel(x, y, z, function(x, y, z));
			}
		}
	}
}

bool testVisibility()
{
	bool uniResult = testVisibilityUnidirectional();
	return uniResult;
}

bool testVisibilityUnidirectional()
{
	std::cout << "Running unidirectional test..." << std::endl;

    Volume* mVolume = nullptr;
    bool loadFromDisk = true;
    if (loadFromDisk)
    {
		mVolume = new Volume("testVisibilityUnidirectional.vol");
    }
    else
    {
		std::cout << "Creating volume..." << std::endl;
		mVolume = new Volume;

		FractalNoise fractalNoise(9);
		Box3i filledBounds(Vector3i::filled(0), Vector3i::filled(511));
		fillVolume(mVolume, filledBounds, fractalNoise);

		std::cout << "Saving volume...";
		mVolume->save("testVisibilityUnidirectional.vol");
		std::cout << "done." << std::endl;
    }

    const float fovyInDegrees = 60.0f;
	// I think the linalg.h perspective_matrix takes it's input in radians, though this isn't documented.
	const float fovyInRadians = fovyInDegrees / 57.2958f;
	/*linalg::mat<float,4,4> projMat = linalg::perspective_matrix<float>(fovyInRadians, 1.0f / 1.0f , 0.1f , 5000.0f);

	linalg::mat<float,4,4> trans = linalg::translation_matrix<float>(linalg::vec<float,3>(0.0, 0.0, -500));
	linalg::mat<float,4,4> rotx = linalg::rotation_matrix<float>(linalg::vec<float,4>(0.392699f, 1.0f, 0.0f, 0.0f));
	linalg::mat<float,4,4> roty = linalg::rotation_matrix<float>(linalg::vec<float,4>(0.392699f, 0.0f, 1.0f, 0.0f));
	linalg::mat<float,4,4> rotz = linalg::rotation_matrix<float>(linalg::vec<float,4>(0.392699f, 0.0f, 0.0f, 1.0f));

	linalg::mat<float,4,4> viewMat = linalg::identity;
	viewMat = mul(viewMat, trans);
	viewMat = mul(viewMat, rotx);
	viewMat = mul(viewMat, roty);
	viewMat = mul(viewMat, rotz);

	std::cout << viewMat << std::endl;

	CameraData cameraData(viewMat, projMat);*/

	CameraData cameraData(Vector3d({ -110, -100, -90 }), Vector3d({ 0, 0, 0 }), Vector3d({ 0.1, 0.2, 0.8 }), fovyInRadians, 1.0);

    VisibilityCalculator visCalc;
	visCalc.mMaxFootprintSize = 0.0055f;

	/*PolygonVertexArray vertices{
		Vector2i(888, 216),
		Vector2i(889, 215),
		Vector2i(885, 216),
		Vector2i(886, 215),
		Vector2i(888, 219),
		Vector2i(888, 218),
		Vector2i(885, 219),
		Vector2i(886, 218) };

	VisibilityMask m(1024, 1024);
	//m.drawQuads(vertices);

	for (int i = 0; i < 8; i++)
	{
		m.drawPixel(vertices[i][0], vertices[i][1]);
	}

	std::cout << "Hash = " << m.hash();
	saveVisibilityMaskAsImage(m, "mask.png");
	exit(1);*/

	Timer timer;
	const int iterations = 100;
	const uint32_t maxGlyphCount = 1000000;
	Glyph* glyphs = new Glyph[maxGlyphCount];
	uint32_t glyphCount = 0;
	for (int ct = 0; ct < iterations; ct++)
	{
		glyphCount = visCalc.findVisibleOctreeNodes(mVolume, &cameraData, NormalEstimation::None, false, glyphs, maxGlyphCount);
	}
	
	float elapsedTime = timer.elapsedTimeInMilliSeconds();
	cout << "\tTime elapsed = " << elapsedTime << "ms." << endl;

	std::cout << "\tFound " << glyphCount << " glyphs." << std::endl;

	saveVisibilityMaskAsImage(*(visCalc.mVisMask), "PerspectiveMask.png");

	delete mVolume;

	uint32_t hash = visCalc.mVisMask->hash();
	std::cout << "\tHash = " << hash << std::endl;

	const size_t expectedGlyphCount = 62117;
	// Tile size affects memory layout and hence hash.
	uint32_t expectedHash = 0;
	if (VisibilityMask::TileSize == 4)
	{
		expectedHash = 4184501084;
	}
	else if (VisibilityMask::TileSize == 8)
	{
		expectedHash = 1810787742;
	}

	check(glyphCount, expectedGlyphCount);
	check(hash, expectedHash);

	bool result = (glyphCount == expectedGlyphCount) && (hash == expectedHash);

    return result;
}
