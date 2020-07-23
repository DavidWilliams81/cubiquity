#include "instance_list.h"

#include "framework/program.h"

#include <algorithm>
#include <iostream>

using namespace Cubiquity;

static const GLfloat gCubeVertices[] =
{
	-0.5f, -0.5f, -0.5f, // triangle 1 : begin
	-0.5f, -0.5f,  0.5f,
	-0.5f,  0.5f,  0.5f, // triangle 1 : end
	 0.5f,  0.5f, -0.5f, // triangle 2 : begin
	-0.5f, -0.5f, -0.5f,
	-0.5f,  0.5f, -0.5f, // triangle 2 : end
	 0.5f, -0.5f,  0.5f,
	-0.5f, -0.5f, -0.5f,
	 0.5f, -0.5f, -0.5f,
	 0.5f,  0.5f, -0.5f,
	 0.5f, -0.5f, -0.5f,
	-0.5f, -0.5f, -0.5f,
	-0.5f, -0.5f, -0.5f,
	-0.5f,  0.5f,  0.5f,
	-0.5f,  0.5f, -0.5f,
	 0.5f, -0.5f,  0.5f,
	-0.5f, -0.5f,  0.5f,
	-0.5f, -0.5f, -0.5f,
	-0.5f,  0.5f,  0.5f,
	-0.5f, -0.5f,  0.5f,
	 0.5f, -0.5f,  0.5f,
	 0.5f,  0.5f,  0.5f,
	 0.5f, -0.5f, -0.5f,
	 0.5f,  0.5f, -0.5f,
	 0.5f, -0.5f, -0.5f,
	 0.5f,  0.5f,  0.5f,
	 0.5f, -0.5f,  0.5f,
	 0.5f,  0.5f,  0.5f,
	 0.5f,  0.5f, -0.5f,
	-0.5f,  0.5f, -0.5f,
	 0.5f,  0.5f,  0.5f,
	-0.5f,  0.5f, -0.5f,
	-0.5f,  0.5f,  0.5f,
	 0.5f,  0.5f,  0.5f,
	-0.5f,  0.5f,  0.5f,
	 0.5f, -0.5f,  0.5f
};

InstanceList::InstanceList()
{
	
}

InstanceList::~InstanceList()
{
	glDeleteVertexArrays(1, &VertexArrayID);
	//delete[] mPerInstanceData;
}

void InstanceList::initialise()
{
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	//mPerInstanceData = new GLfloat[MaxParticles * 8];

	mCubeVerticesBuffer.create();
	mCubeVerticesBuffer.bind();
	mCubeVerticesBuffer.allocate(gCubeVertices, sizeof(gCubeVertices));

	// The VBO containing the positions and sizes of the particles
	mPerInstanceDataBuffer.setUsagePattern(OpenGLBuffer::StreamDraw);
	mPerInstanceDataBuffer.create();
	mPerInstanceDataBuffer.bind();

	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	mPerInstanceDataBuffer.allocate(MaxParticles * 8 * sizeof(GLfloat));

	
}

void InstanceList::render()
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	mCubeVerticesBuffer.bind();
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
	mPerInstanceDataBuffer.bind();
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
	mPerInstanceDataBuffer.bind();
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
	glDrawArraysInstanced(GL_TRIANGLES, 0, 12 * 3, ParticlesCount);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
}

void InstanceList::setPerInstanceData(const Cubiquity::Glyph* glyphs, uint32_t glyphCount)
{
	if (glyphCount > MaxParticles)
	{
		std::cerr << "Too many glyphs to upload!" << std::endl;
	}

	unsigned int glyphsToCopy = std::min(glyphCount, MaxParticles);

	if (glyphsToCopy > 0)
	{
		// Update the buffers that OpenGL uses for rendering.
		// There are much more sophisticated means to stream data from the CPU to the GPU, 
		// but this is outside the scope of this tutorial.
		// http://www.opengl.org/wiki/Buffer_Object_Streaming
		mPerInstanceDataBuffer.bind();
		mPerInstanceDataBuffer.allocate(MaxParticles * 8 * sizeof(GLfloat));
		mPerInstanceDataBuffer.write(0, &(glyphs[0]), glyphsToCopy * sizeof(GLfloat) * 8);  // Buffer orphaning, a common way to improve streaming perf. See above link for details.
	}
	else
	{
		//std::cerr << "Uploaded zero glyphs" << std::endl;
	}
	
	ParticlesCount = glyphsToCopy;
}
