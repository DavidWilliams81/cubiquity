#include "volume_renderer.h"

#include "rendering.h"
#include "utility.h"

#include <iostream>

using namespace Cubiquity;

const uint32_t maxGlyphCount = 1000000;

VolumeRenderer::VolumeRenderer(const MaterialSet& materials)
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

	// Create one OpenGL texture
	glGenTextures(1, &Texture);

	// "Bind" the newly created texture : all future texture functions will modify this texture
	glBindTexture(GL_TEXTURE_1D, Texture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, materials.data().size(), 0, GL_RGB, GL_FLOAT, materials.data().data());

	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	// Get a handle for our "materials" uniform
	TextureID = glGetUniformLocation(program.programID, "materials");

	mVisibilityCalculator = new VisibilityCalculator;
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
	// Settings for debugging a difficult angle on the 'building' map.
	// 'camera' needs to be changed to a non-const non-ref to use these.
	//std::cout << camera.position << " : " << camera.pitch << " : " << camera.yaw << std::endl;
	//camera.position = Vector3f(262.409, 320.604, 187.974);
	//camera.pitch = -0.00539806;
	//camera.yaw = -0.92;
	
	cameraData = CameraData(camera.position, camera.position + camera.forward(), camera.up(), camera.fovInDegrees / 57.2958f, camera.aspect);

	const Vector3f volumeCentre = Vector3f({ 0.0f, 0.0f, 0.0f });

	mVisibilityCalculator->mMaxFootprintSize = 0.01f;
	Timer timer;
	uint32_t glyphCount = mVisibilityCalculator->findVisibleOctreeNodes(&(cameraData), mVolume, mGlyphs, maxGlyphCount);
	std::cout << "Found " << glyphCount << " glyphs in " << timer.elapsedTimeInMilliSeconds() << "ms" << std::endl;

	// We render back-to-front in order to support alpha blending for points, but for cubes we render
	// front to back and use the depth buffer to discard fragments (probably not needed, to be tested).
	if (mGlyphType == GlyphType::Point)
	{
		std::reverse(mGlyphs, mGlyphs + glyphCount);
	}

	if (mDoGlyphUpdates)
	{
		instanceList.setPerInstanceData(mGlyphs, glyphCount);
	}

	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	Matrix4x4f ProjectionMatrix = static_cast<Matrix4x4f>(camera.projectionMatrix());
	Matrix4x4f ViewMatrix = static_cast<Matrix4x4f>(camera.viewMatrix());

	// Use our shader
	glUseProgram(program.programID);

	glUniformMatrix4fv(viewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
	glUniformMatrix4fv(projMatrixID, 1, GL_FALSE, &ProjectionMatrix[0][0]);
	glUniform3f(cameraPosID, camera.position.x(), camera.position.y(), camera.position.z());

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_1D, Texture);
	// Set our "myTextureSampler" sampler to use Texture Unit 0
	glUniform1i(TextureID, 0);

	if (instanceList.ParticlesCount > 0)
	{
		Cubiquity::Matrix4x4f ModelMatrix = translation_matrix(volumeCentre);
		glUniformMatrix4fv(modelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		instanceList.render();
	}
}
