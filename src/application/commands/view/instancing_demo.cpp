#include "instancing_demo.h"

#include "base/logging.h"

#include "utility.h"

#include "extraction.h"

using namespace Cubiquity;

const uint32_t MaxGlyphCount = 10000000;

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

	std::string shader_path = getShaderPath();
	log_debug("Using shader path'{}", shader_path);

	// Create and compile our GLSL program from the shaders
	glyphProgram = loadProgram(
		shader_path + std::string("glyph.vert"),
		shader_path + std::string("glyph.frag"));

	// Set up access to uniforms
	modelMatrixID = glGetUniformLocation(glyphProgram, "modelMatrix");
	viewMatrixID = glGetUniformLocation(glyphProgram, "viewMatrix");
	projMatrixID = glGetUniformLocation(glyphProgram, "projectionMatrix");
	cameraPosID = glGetUniformLocation(glyphProgram, "cameraPos");

	// Set up the texture which holds our voxel materials
	glGenTextures(1, &materialsTexture);
	glBindTexture(GL_TEXTURE_1D, materialsTexture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, colours().size(), 0, GL_RGB, GL_FLOAT, colours().data());
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
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

	// Light blue background
	glClearColor(0.8f, 0.8f, 1.0f, 0.0f);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
}

void InstancingDemo::onUpdate(float deltaTime)
{
	OpenGLViewer::onUpdate(deltaTime);

	const vec3f volumeCentre = vec3f({ 0.0f, 0.0f, 0.0f });

	Timer timer;

	NormalEstimation normalEstimation = NormalEstimation::None;

	// Per-glyph normals looks poor if they are too large, so subdivide.
	bool subdivideMaterialNodes = normalEstimation != NormalEstimation::None;

	if (mNeedsUpdate)
	{
		//mGlyphCount = mVisibilityCalculator->findVisibleOctreeNodes(&(volume()), camera().position, normalEstimation, subdivideMaterialNodes, mGlyphs, MaxGlyphCount);

		mGlyphCount = extractGlyphs(const_cast<Volume&>(volume()), subdivideMaterialNodes, mGlyphs, MaxGlyphCount);

		log_debug("Found {} glyphs in {}ms", mGlyphCount, timer.elapsedTimeInMilliSeconds());
		assert(mGlyphCount <= MaxGlyphCount);

		// Update the buffers that OpenGL uses for rendering.
		// See http://www.opengl.org/wiki/Buffer_Object_Streaming for other approaches
		glBindBuffer(GL_ARRAY_BUFFER, mPerInstanceDataBuffer);
		glBufferData(GL_ARRAY_BUFFER, MaxGlyphCount * sizeof(Glyph), NULL, GL_STREAM_DRAW); // Buffer orphaning
		glBufferSubData(GL_ARRAY_BUFFER, 0, mGlyphCount * sizeof(Glyph), mGlyphs);

		mNeedsUpdate = false;
	}

	Matrix4x4f ProjectionMatrix = static_cast<Matrix4x4f>(camera().projectionMatrix());
	Matrix4x4f ViewMatrix = static_cast<Matrix4x4f>(camera().viewMatrix());

	// Use our shader
	glUseProgram(glyphProgram);

	glUniformMatrix4fv(viewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
	glUniformMatrix4fv(projMatrixID, 1, GL_FALSE, &ProjectionMatrix[0][0]);
	glUniform3f(cameraPosID, camera().position.x, camera().position.y, camera().position.z);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_1D, materialsTexture);
	// Set our "myTextureSampler" sampler to use Texture Unit 0
	glUniform1i(materialsTextureID, 0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// We draw all glyphs in a single drawcall and OpenGL does not guaentee the ordering
	// of these. Therefore we need to enable the depth buffer for a correct result.
	// https://community.khronos.org/t/rendering-order-within-a-single-draw-call/66591/4
	glEnable(GL_DEPTH_TEST);

	Matrix4x4f ModelMatrix = translation_matrix(volumeCentre);
	glUniformMatrix4fv(modelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

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
	OpenGLViewer::onKeyDown(event);

	if (event.keysym.sym == SDLK_SPACE)
	{
		mNeedsUpdate = !(mNeedsUpdate);
	}
}

void InstancingDemo::onVolumeModified()
{
	//subDAGs = findSubDAGs(
	//	Internals::getNodes(volume()).nodes(), getRootNodeIndex(volume()));

	//clear();
	mNeedsUpdate = true;
}