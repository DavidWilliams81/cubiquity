#include "instancing_demo.h"

#include "utility.h"

#include <iostream>

#include "rendering.h"

using namespace Cubiquity;

const uint32_t MaxGlyphCount = 1000000;

static const GLfloat gCubeVertices[] =
{
	-0.5f, -0.5f, -0.5f, -0.5f, -0.5f,  0.5f, -0.5f,  0.5f,  0.5f,
	 0.5f,  0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,  0.5f, -0.5f,
	 0.5f, -0.5f,  0.5f, -0.5f, -0.5f, -0.5f,  0.5f, -0.5f, -0.5f,
	 0.5f,  0.5f, -0.5f,  0.5f, -0.5f, -0.5f, -0.5f, -0.5f, -0.5f,
	-0.5f, -0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f,  0.5f, -0.5f,
	 0.5f, -0.5f,  0.5f, -0.5f, -0.5f,  0.5f, -0.5f, -0.5f, -0.5f,
	-0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f,  0.5f,
	 0.5f,  0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f, -0.5f,
	 0.5f, -0.5f, -0.5f,  0.5f,  0.5f,  0.5f,  0.5f, -0.5f,  0.5f,
	 0.5f,  0.5f,  0.5f,  0.5f,  0.5f, -0.5f, -0.5f,  0.5f, -0.5f,
	 0.5f,  0.5f,  0.5f, -0.5f,  0.5f, -0.5f, -0.5f,  0.5f,  0.5f,
	 0.5f,  0.5f,  0.5f, -0.5f,  0.5f,  0.5f,  0.5f, -0.5f,  0.5f
};

void InstancingDemo::onInitialise()
{
	OpenGLViewer::onInitialise();

	mVisibilityCalculator = new VisibilityCalculator;

	// Create and compile our GLSL program from the shaders
	std::cout << "Warning - Using hard-coded paths to shaders in ../src/application/commands/view" << std::endl;
	glyphProgram = loadProgram("../src/application/commands/view/glyph_vertex_program.glsl", "../src/application/commands/view/glyph_fragment_program.glsl");
	screenQuadProgram = loadProgram("../src/application/commands/view/screen_quad_vertex_program.glsl", "../src/application/commands/view/screen_quad_fragment_program.glsl");

	// Set up access to uniforms
	modelMatrixID = glGetUniformLocation(glyphProgram, "modelMatrix");
	viewMatrixID = glGetUniformLocation(glyphProgram, "viewMatrix");
	projMatrixID = glGetUniformLocation(glyphProgram, "projectionMatrix");
	modeID = glGetUniformLocation(glyphProgram, "mode");
	cameraPosID = glGetUniformLocation(glyphProgram, "cameraPos");

	// Set up the texture which holds our voxel materials
	glGenTextures(1, &materialsTexture);
	glBindTexture(GL_TEXTURE_1D, materialsTexture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, materials().data().size(), 0, GL_RGB, GL_FLOAT, materials().data().data());
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	materialsTextureID = glGetUniformLocation(glyphProgram, "materials");


	//Init instance list
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	//mPerInstanceData = new GLfloat[MaxParticles * 8];

	glGenBuffers(1, &mCubeVerticesBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, mCubeVerticesBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(gCubeVertices), gCubeVertices, GL_STATIC_DRAW);

	// The VBO containing the positions and sizes of the particles
	glGenBuffers(1, &mPerInstanceDataBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, mPerInstanceDataBuffer);

	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, MaxGlyphCount * 8 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);




	mGlyphs = new Glyph[MaxGlyphCount];

	glGenFramebuffers(1, &framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

	// generate texture
	glGenTextures(1, &textureColorbuffer);
	glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 1600, 1200, 0, GL_RGBA, GL_FLOAT, NULL); // FIXME - Use correct dimensions
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	// attach it to currently bound framebuffer object
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);

	unsigned int rbo;
	glGenRenderbuffers(1, &rbo);
	glBindRenderbuffer(GL_RENDERBUFFER, rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, 1600, 1200); // use a single renderbuffer object for both a depth AND stencil buffer.
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo); // now actually attach it

	glBindFramebuffer(GL_FRAMEBUFFER, 0);



	// Black background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
}

void InstancingDemo::onUpdate(float deltaTime)
{
	OpenGLViewer::onUpdate(deltaTime);

	CameraData cameraData(camera().position, camera().position + camera().forward(), camera().up(), camera().fovInDegrees / 57.2958f, camera().aspect);

	const Vector3f volumeCentre = Vector3f({ 0.0f, 0.0f, 0.0f });

	mVisibilityCalculator->mMaxFootprintSize = 0.007f;
	Timer timer;

	NormalEstimation normalEstimation = mGlyphType == GlyphType::Disc ?
		NormalEstimation::FromChildren : NormalEstimation::None;

	// Disc glyphs looks poor if they are too large, so subdivide.
	// Likewise for per-glyph normals, even if they are applied to cubic glyphs.
	bool subdivideMaterialNodes = (mGlyphType == GlyphType::Disc) || (normalEstimation != NormalEstimation::None);

	if (mDoGlyphUpdates)
	{
		mGlyphCount = mVisibilityCalculator->findVisibleOctreeNodes(&(volume()), &(cameraData), normalEstimation, subdivideMaterialNodes, mGlyphs, MaxGlyphCount);
		std::cout << "Found " << mGlyphCount << " glyphs in " << timer.elapsedTimeInMilliSeconds() << "ms" << std::endl;
		assert(mGlyphCount <= MaxGlyphCount);

		// Update the buffers that OpenGL uses for rendering.
		// See http://www.opengl.org/wiki/Buffer_Object_Streaming for other approaches
		glBindBuffer(GL_ARRAY_BUFFER, mPerInstanceDataBuffer);
		glBufferData(GL_ARRAY_BUFFER, MaxGlyphCount * sizeof(Glyph), NULL, GL_STREAM_DRAW); // Buffer orphaning
		glBufferSubData(GL_ARRAY_BUFFER, 0, mGlyphCount * sizeof(Glyph), mGlyphs);
	}

	Matrix4x4f ProjectionMatrix = static_cast<Matrix4x4f>(camera().projectionMatrix());
	Matrix4x4f ViewMatrix = static_cast<Matrix4x4f>(camera().viewMatrix());

	// Use our shader
	glUseProgram(glyphProgram);

	glUniformMatrix4fv(viewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
	glUniformMatrix4fv(projMatrixID, 1, GL_FALSE, &ProjectionMatrix[0][0]);
	glUniform1ui(modeID, (unsigned int)mGlyphType);
	glUniform3f(cameraPosID, camera().position.x(), camera().position.y(), camera().position.z());

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_1D, materialsTexture);
	// Set our "myTextureSampler" sampler to use Texture Unit 0
	glUniform1i(materialsTextureID, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// We draw all glyphs in a single drawcall and OpenGL does not guaentee the ordering
	// of these. Therefore we need to enable the depth buffer for a correct result.
	// https://community.khronos.org/t/rendering-order-within-a-single-draw-call/66591/4
	// This only applie to cube glyphs (as discs are simply additively blended).
	mGlyphType == GlyphType::Cube ? glEnable(GL_DEPTH_TEST) : glDisable(GL_DEPTH_TEST);

	Matrix4x4f ModelMatrix = translation_matrix(volumeCentre);
	glUniformMatrix4fv(modelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

	if (mGlyphType == GlyphType::Disc)
	{
		glEnable(GL_BLEND);
	}
	glBlendFunc(GL_ONE, GL_ONE); // Additive blending

	glBindVertexArray(VertexArrayID);

	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, mCubeVerticesBuffer);
	glVertexAttribPointer(
		0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
		3,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
	);

	// 2nd attribute buffer : positions of particles' centers
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, mPerInstanceDataBuffer);
	glVertexAttribPointer(
		1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
		4,                                // size : x + y + z + size => 4
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		sizeof(GLfloat) * 8,              // stride
		(void*)0                          // array buffer offset
	);

	// 3rd attribute buffer : particles' colors
	glEnableVertexAttribArray(2);
	glBindBuffer(GL_ARRAY_BUFFER, mPerInstanceDataBuffer);
	glVertexAttribPointer(
		2,                                // attribute. No particular reason for 1, but must match the layout in the shader.
		4,                                // size : r + g + b + a => 4
		GL_FLOAT,						  // type
		GL_FALSE,                         // normalized?    *** YES, this means that the unsigned char[4] will be accessible with a vec4 (floats) in the shader ***
		sizeof(GLfloat) * 8,              // stride
		(void*)16                         // array buffer offset
	);

	// These functions are specific to glDrawArrays*Instanced*.
	// The first parameter is the attribute buffer we're talking about.
	// The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
	// http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
	glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
	glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
	glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1

									// Draw the particules !
									// This draws many times a small triangle_strip (which looks like a quad).
									// This is equivalent to :
									// for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4), 
									// but faster.
	glDrawArraysInstanced(GL_TRIANGLES, 0, 12 * 3, mGlyphCount);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);

	glDisable(GL_BLEND);

	glBindFramebuffer(GL_FRAMEBUFFER, 0); // back to default

	//screenShader.use();
	glUseProgram(screenQuadProgram);
	glUniform1i(glGetUniformLocation(screenQuadProgram, "screenTexture"), 0);
	
	glBindTexture(GL_TEXTURE_2D, textureColorbuffer);	// use the color attachment texture as the texture of the quad plane
	glDisable(GL_DEPTH_TEST);
	drawScreenAlignedQuad();

	glCheckError();
}

void InstancingDemo::onShutdown()
{
	OpenGLViewer::onShutdown();

	delete[] mGlyphs;

	glDeleteVertexArrays(1, &VertexArrayID);

	glDeleteProgram(glyphProgram);

	//delete mVolumeRenderer;
}

void InstancingDemo::onKeyDown(const SDL_KeyboardEvent & event)
{
	OpenGLViewer::onKeyUp(event);

	if (event.keysym.sym == SDLK_SPACE)
	{
		mDoGlyphUpdates = !(mDoGlyphUpdates);
	}

	if (event.keysym.sym == SDLK_t)
	{
		mGlyphType = mGlyphType == GlyphType::Cube ? GlyphType::Disc : GlyphType::Cube;
	}
}
