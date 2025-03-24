#include "gpu_pathtracing_viewer.h"

#include "base/logging.h"

#include "utility.h"

#include "visibility.h"
#include "raytracing.h"

using namespace Cubiquity;

void GPUPathtracingViewer::onInitialise()
{
	OpenGLViewer::onInitialise();

	std::string shader_path = getShaderPath();
	log_debug("Using shader path'{}", shader_path);

	// Create and compile our GLSL program from the shaders
	progressiveProgram = loadProgram(
		shader_path + std::string("screen_aligned_quad.vert"),
		shader_path + std::string("pathtracing.frag"),
		"#define PROGRESSIVE");
	previewProgram = loadProgram(
		shader_path + std::string("screen_aligned_quad.vert"),
		shader_path + std::string("pathtracing.frag"));
	copyProgram = loadProgram(
		shader_path + std::string("screen_aligned_quad.vert"),
		shader_path + std::string("normalise.frag"));
	hBlurProgram = loadProgram(
		shader_path + std::string("screen_aligned_quad.vert"),
		shader_path + std::string("horz_blur.frag"));
	vBlurProgram = loadProgram(
		shader_path + std::string("screen_aligned_quad.vert"),
		shader_path + std::string("vert_blur.frag"));

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
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(Node) * nodeCount, nodeStore.data(), GL_STATIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, dagData);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // unbind

	GLuint subDagData;
	glGenBuffers(1, &subDagData);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, subDagData);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(subDAGs), &subDAGs, GL_STATIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, subDagData);
	//glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0); // unbind

	// Set up the texture which holds our voxel materials
	glGenTextures(1, &mMaterialsTexture);
	glBindTexture(GL_TEXTURE_1D, mMaterialsTexture);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, colours().size(), 0, GL_RGB, GL_FLOAT, colours().data());
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);


	glGenFramebuffers(1, &accFramebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, accFramebuffer);

	// generate texture
	glGenTextures(1, &accTexture);
	glBindTexture(GL_TEXTURE_2D, accTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width(), height(), 0, GL_RGBA, GL_FLOAT, NULL); // FIXME - Use correct dimensions
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	// attach it to currently bound framebuffer object
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, accTexture, 0);

	unsigned int rbo;
	glGenRenderbuffers(1, &rbo);
	glBindRenderbuffer(GL_RENDERBUFFER, rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width(), height()); // use a single renderbuffer object for both a depth AND stencil buffer.
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo); // now actually attach it

	glBindFramebuffer(GL_FRAMEBUFFER, 0);




	glGenFramebuffers(1, &blurFramebuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, blurFramebuffer);

	// generate texture
	glGenTextures(1, &blurTexture);
	glBindTexture(GL_TEXTURE_2D, blurTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width(), height(), 0, GL_RGBA, GL_FLOAT, NULL); // FIXME - Use correct dimensions
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);

	// attach it to currently bound framebuffer object
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, blurTexture, 0);

	unsigned int rbo2;
	glGenRenderbuffers(1, &rbo2);
	glBindRenderbuffer(GL_RENDERBUFFER, rbo2);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width(), height()); // use a single renderbuffer object for both a depth AND stencil buffer.
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo2); // now actually attach it

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	mTimer.start();
}

void GPUPathtracingViewer::onUpdate(float deltaTime)
{
	OpenGLViewer::onUpdate(deltaTime);

	// Switch to progressive mode for a high quality render once the camera is static.
	// The camera needs to be static because we do not clear the window's framebuffer
	// before drawing, instead we slowly draw the progressive version over the 
	// preview. This is nice because then we do not see black holes where the
	// progressive version has not yet drawn. I *think* the threshold here needs to
	// be set to '2' because we have double buffering enabled and need to make sure
	// both buffers have the same static preview image drawn into them before we do
	// the progressive render over the top. An alternative might be to post-process
	// into a 'single-buffered' target seperate from the window buffer, and then do
	// the copy to the double-buffered window buffer from there.
	bool progressive = staticFrameCount >= 2;

	Matrix4x4f PInv = static_cast<Matrix4x4f>(inverse(camera().projectionMatrix()));
	Matrix4x4f VInv = static_cast<Matrix4x4f>(inverse(camera().viewMatrix()));

	glBindFramebuffer(GL_FRAMEBUFFER, accFramebuffer);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	GLuint pathtracingProgram = progressive ? progressiveProgram : previewProgram;
	glUseProgram(pathtracingProgram);
	glUniformMatrix4fv(glGetUniformLocation(pathtracingProgram, "VInv"), 1, GL_FALSE, &VInv[0][0]);
	glUniformMatrix4fv(glGetUniformLocation(pathtracingProgram, "PInv"), 1, GL_FALSE, &PInv[0][0]);
	glUniform3f(glGetUniformLocation(pathtracingProgram, "cameraPos"), camera().position.x, camera().position.y, camera().position.z);
	glUniform1ui(glGetUniformLocation(pathtracingProgram, "frameId"), ++frameId);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_1D, mMaterialsTexture);
	// Set our texture sampler sampler to use Texture Unit 0
	glUniform1i(glGetUniformLocation(pathtracingProgram, "materials"), 0);

	if (progressive)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE, GL_ONE); // Additive blending
	}
	else
	{
		glDisable(GL_BLEND);
	}

	glDisable(GL_DEPTH_TEST);
	drawScreenAlignedQuad();

	/*int blurPasses = 0;
	for (int i = 0; i < blurPasses; i++)
	{
		glBindFramebuffer(GL_FRAMEBUFFER, blurFramebuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//screenShader.use();
		glUseProgram(hBlurProgram);
		glUniform1i(glGetUniformLocation(hBlurProgram, "screenTexture"), 0);

		glBindTexture(GL_TEXTURE_2D, accTexture);	// use the color attachment texture as the texture of the quad plane
		glDisable(GL_DEPTH_TEST);
		drawScreenAlignedQuad();

		glBindFramebuffer(GL_FRAMEBUFFER, accFramebuffer);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//screenShader.use();
		glUseProgram(vBlurProgram);
		glUniform1i(glGetUniformLocation(vBlurProgram, "screenTexture"), 0);

		glBindTexture(GL_TEXTURE_2D, blurTexture);	// use the color attachment texture as the texture of the quad plane
		glDisable(GL_DEPTH_TEST);
		drawScreenAlignedQuad();
	}*/

	glBindFramebuffer(GL_FRAMEBUFFER, 0); // back to default
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//screenShader.use();
	glUseProgram(copyProgram);
	glUniform1i(glGetUniformLocation(copyProgram, "screenTexture"), 0);

	glBindTexture(GL_TEXTURE_2D, accTexture);	// use the color attachment texture as the texture of the quad plane
	glDisable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	drawScreenAlignedQuad();
	glDisable(GL_BLEND);

	// In non-progressive mode we have now finished with the accumulation buffer.
	// We clear it because the next frame might be progressive and we don't want
	// to add to the preview image (we should replace it).
	if (!progressive)
	{
		glBindFramebuffer(GL_FRAMEBUFFER, accFramebuffer);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}

	glCheckError();

	staticFrameCount++;

	const int groupSize = 100;
	if (frameNumber() % groupSize == 0)
	{
		log_info("Frame time = {} ms", mTimer.elapsedTimeInMilliSeconds() / groupSize);
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

void GPUPathtracingViewer::onMouseButtonDown(const SDL_MouseButtonEvent& event)
{
	OpenGLViewer::onMouseButtonDown(event);
}

void GPUPathtracingViewer::onMouseButtonUp(const SDL_MouseButtonEvent& event)
{
	OpenGLViewer::onMouseButtonUp(event);
}

void GPUPathtracingViewer::onCameraModified()
{
	//cameraModified = true;
	staticFrameCount = 0;
	//clear();
}

void GPUPathtracingViewer::clear()
{
	//glBindFramebuffer(GL_FRAMEBUFFER, accFramebuffer);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glBindFramebuffer(GL_FRAMEBUFFER, 0);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
