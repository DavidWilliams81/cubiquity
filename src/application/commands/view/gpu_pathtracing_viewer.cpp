#include "gpu_pathtracing_viewer.h"

#include "utility.h"

#include <iostream>

#include "rendering.h"

using namespace Cubiquity;

void GPUPathtracingViewer::onInitialise()
{
	OpenGLViewer::onInitialise();

	// Create and compile our GLSL program from the shaders
	std::cout << "Warning - Using hard-coded paths to shaders in ../src/application/commands/view" << std::endl;
	mainProgram = loadProgram(
		"../src/application/commands/view/glsl/gpu_pathtracing_main_vertex_program.glsl",
		"../src/application/commands/view/glsl/gpu_pathtracing_main_fragment_program.glsl");
	screenQuadProgram = loadProgram(
		"../src/application/commands/view/glsl/gpu_pathtracing_screen_quad_vertex_program.glsl",
		"../src/application/commands/view/glsl/gpu_pathtracing_screen_quad_fragment_program.glsl");

	// Black background
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	const Internals::NodeStore& nodeStore = Internals::getNodes(volume()).nodes();
	SubDAGArray subDAGs = findSubDAGs(nodeStore, getRootNodeIndex(volume()));

	uint32 nodeCount = Internals::getNodes(volume()).bakedNodesEnd();

	GLuint dagData;
	glGenBuffers(1, &dagData);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, dagData);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(Node) * nodeCount, nodeStore.data(), GL_DYNAMIC_COPY);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, dagData);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // unbind

	GLuint subDagData;
	glGenBuffers(1, &subDagData);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, subDagData);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(subDAGs), &subDAGs, GL_DYNAMIC_COPY);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, subDagData);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // unbind

	// Set up the texture which holds our voxel materials
	glGenTextures(1, &mMaterialsTexture);
	glBindTexture(GL_TEXTURE_1D, mMaterialsTexture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, materials().data().size(), 0, GL_RGB, GL_FLOAT, materials().data().data());
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


	glGenFramebuffers(1, &framebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);

	// generate texture
	glGenTextures(1, &textureColorbuffer);
	glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width(), height(), 0, GL_RGBA, GL_FLOAT, NULL); // FIXME - Use correct dimensions
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	// attach it to currently bound framebuffer object
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);

	unsigned int rbo;
	glGenRenderbuffers(1, &rbo);
	glBindRenderbuffer(GL_RENDERBUFFER, rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width(), height()); // use a single renderbuffer object for both a depth AND stencil buffer.
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo); // now actually attach it

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	mTimer.start();
}

void GPUPathtracingViewer::onUpdate(float deltaTime)
{
	OpenGLViewer::onUpdate(deltaTime);

	Matrix4x4f PInv = static_cast<Matrix4x4f>(inverse(camera().projectionMatrix()));
	Matrix4x4f VInv = static_cast<Matrix4x4f>(inverse(camera().viewMatrix()));

	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(mainProgram);
	glUniformMatrix4fv(glGetUniformLocation(mainProgram, "VInv"), 1, GL_FALSE, &VInv[0][0]);
	glUniformMatrix4fv(glGetUniformLocation(mainProgram, "PInv"), 1, GL_FALSE, &PInv[0][0]);
	glUniform3f(glGetUniformLocation(mainProgram, "cameraPos"), camera().position.x(), camera().position.y(), camera().position.z());
	glUniform1ui(glGetUniformLocation(mainProgram, "randSeed"), randSeed);
	randSeed = mix(randSeed);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_1D, mMaterialsTexture);
	// Set our texture sampler sampler to use Texture Unit 0
	glUniform1i(glGetUniformLocation(mainProgram, "materials"), 0);

	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE); // Additive blending
	
	glDisable(GL_DEPTH_TEST);
	drawScreenAlignedQuad();

	glDisable(GL_BLEND);

	glBindFramebuffer(GL_FRAMEBUFFER, 0); // back to default
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//screenShader.use();
	glUseProgram(screenQuadProgram);
	glUniform1i(glGetUniformLocation(screenQuadProgram, "screenTexture"), 0);

	glBindTexture(GL_TEXTURE_2D, textureColorbuffer);	// use the color attachment texture as the texture of the quad plane
	glDisable(GL_DEPTH_TEST);
	drawScreenAlignedQuad();

	glCheckError();

	const int groupSize = 100;
	if (frameNumber() % groupSize == 0)
	{
		std::cout << "Frame time = " << mTimer.elapsedTimeInMilliSeconds() / groupSize << "ms" << std::endl;
		mTimer.start(); // Reset
	}
}

void GPUPathtracingViewer::onShutdown()
{
	OpenGLViewer::onShutdown();
}

void GPUPathtracingViewer::onKeyDown(const SDL_KeyboardEvent & event)
{
	OpenGLViewer::onKeyDown(event);
}

void GPUPathtracingViewer::onCameraModified()
{
	clear();
}

void GPUPathtracingViewer::clear()
{
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}