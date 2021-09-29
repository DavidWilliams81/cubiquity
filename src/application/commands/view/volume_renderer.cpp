#include "volume_renderer.h"

#include "rendering.h"

#include <iostream>

using namespace Cubiquity;

const uint32_t maxGlyphCount = 1000000;

VolumeRenderer::VolumeRenderer()
{
	mGlyphType = GlyphType::Cube;

	// Should fix this, but I'll thank myself for the warning in the meantime...
	std::cout << "Warning - Using hard-coded paths to shaders in ../src/demos" << std::endl;

	// Create and compile our GLSL program from the shaders
	switch (mGlyphType)
	{
	case GlyphType::Cube:
		program.LoadShaders("../src/application/commands/view/cube_vertex_program.glsl", "../src/application/commands/view/cube_fragment_program.glsl");
		break;
	case GlyphType::Point:
		program.LoadShaders("../src/application/commands/view/point_vertex_program.glsl", "../src/application/commands/view/point_fragment_program.glsl");
		break;
	default:
		std::cerr << "Invalid glyph type" << std::endl;
	}

	modelMatrixID = glGetUniformLocation(program.programID, "modelMatrix");
	viewMatrixID = glGetUniformLocation(program.programID, "viewMatrix");
	projMatrixID = glGetUniformLocation(program.programID, "projectionMatrix");
	cameraPosID = glGetUniformLocation(program.programID, "cameraPos");

	mVisibilityCalculator = new VisibilityCalculator;
	mVisibilityCalculator->mMaxNodeDistance = 3000.0f;
	instanceList.initialise();

	mGlyphs = new Glyph[maxGlyphCount];
}

VolumeRenderer::~VolumeRenderer()
{
	delete[] mGlyphs;

	glDeleteProgram(program.programID);
}

void VolumeRenderer::render(const Camera& camera)
{
	cameraData = CameraData(camera.position, camera.position + camera.forward(), camera.up(), camera.fovInDegrees / 57.2958f, camera.aspect);

	const Vector3f volumeCentre = Vector3f(0.0f, 0.0f, 0.0f);

	mVisibilityCalculator->mMaxFootprintSize = 0.01f;
	uint32_t glyphCount = mVisibilityCalculator->findVisibleOctreeNodesPerspective(&(cameraData), mVolume, mGlyphs, maxGlyphCount);

	// We render back-to-front in order to support alpha blending for points, but for cubes we render
	// front to back and use the depth buffer to discard fragments (probably not needed, to be tested).
	if (mGlyphType == GlyphType::Point)
	{
		std::reverse(mGlyphs, mGlyphs + glyphCount);
	}

	instanceList.setPerInstanceData(mGlyphs, glyphCount);

	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	Matrix4x4f ProjectionMatrix = camera.projectionMatrix();
	Matrix4x4f ViewMatrix = camera.viewMatrix();

	// Use our shader
	glUseProgram(program.programID);

	glUniformMatrix4fv(viewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
	glUniformMatrix4fv(projMatrixID, 1, GL_FALSE, &ProjectionMatrix[0][0]);
	glUniform3f(cameraPosID, camera.position.x(), camera.position.y(), camera.position.z());

	if (instanceList.ParticlesCount > 0)
	{
		Cubiquity::Matrix4x4f ModelMatrix = translation_matrix(volumeCentre);
		glUniformMatrix4fv(modelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		instanceList.render();
	}
}
