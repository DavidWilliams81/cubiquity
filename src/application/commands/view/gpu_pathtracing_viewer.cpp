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
	screenQuadProgram = loadProgram("../src/application/commands/view/gpu_pathtracing_screen_quad_vertex_program.glsl", "../src/application/commands/view/gpu_pathtracing_screen_quad_fragment_program.glsl");

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

	mTimer.start();
}

void GPUPathtracingViewer::onUpdate(float deltaTime)
{
	OpenGLViewer::onUpdate(deltaTime);

	Matrix4x4f PInv = static_cast<Matrix4x4f>(inverse(camera().projectionMatrix()));
	Matrix4x4f VInv = static_cast<Matrix4x4f>(inverse(camera().viewMatrix()));

	glUseProgram(screenQuadProgram);
	glUniformMatrix4fv(glGetUniformLocation(screenQuadProgram, "VInv"), 1, GL_FALSE, &VInv[0][0]);
	glUniformMatrix4fv(glGetUniformLocation(screenQuadProgram, "PInv"), 1, GL_FALSE, &PInv[0][0]);
	glUniform3f(glGetUniformLocation(screenQuadProgram, "cameraPos"), camera().position.x(), camera().position.y(), camera().position.z());

	// Bind our texture in Texture Unit 1
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_1D, mMaterialsTexture);
	// Set our texture sampler sampler to use Texture Unit 1
	glUniform1i(glGetUniformLocation(screenQuadProgram, "materials"), 1);
	
	glDisable(GL_DEPTH_TEST);
	drawScreenAlignedQuad();

	glCheckError();

	std::cout << "Frame time = " << mTimer.elapsedTimeInMilliSeconds() << "ms" << std::endl;
	mTimer.start(); // Reset
}

void GPUPathtracingViewer::onShutdown()
{
	OpenGLViewer::onShutdown();
}

void GPUPathtracingViewer::onKeyDown(const SDL_KeyboardEvent & event)
{
	OpenGLViewer::onKeyDown(event);
}
