#ifndef INSTANCE_LIST_H_64D74DE1
#define INSTANCE_LIST_H_64D74DE1

#include "framework/opengl_buffer.h"

#include <vector>

#include <glad/glad.h>

#include "rendering.h"

const uint32_t MaxParticles = 1000000;

class InstanceList
{
public:
	InstanceList();
	~InstanceList();

	void initialise();
	void render();

	void setPerInstanceData(const Cubiquity::Glyph* glyphs, uint32_t glyphCount);

	GLuint VertexArrayID;

	//GLfloat* mPerInstanceData;

	OpenGLBuffer mCubeVerticesBuffer;

	//GLuint mPerInstanceDataBuffer;
	OpenGLBuffer mPerInstanceDataBuffer;


	int ParticlesCount;
};

#endif // INSTANCE_LIST_H_64D74DE1
