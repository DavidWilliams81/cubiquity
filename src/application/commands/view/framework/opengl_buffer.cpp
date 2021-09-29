#include "opengl_buffer.h"

OpenGLBuffer::OpenGLBuffer()
	:mId(0)
	,mUsagePattern(StaticDraw)
{
}

OpenGLBuffer::~OpenGLBuffer()
{
	destroy();
}

void OpenGLBuffer::allocate(int count)
{
	glBufferData(GL_ARRAY_BUFFER, count, NULL, mUsagePattern);
}

void OpenGLBuffer::allocate(const void *data, int count)
{
	glBufferData(GL_ARRAY_BUFFER, count, data, mUsagePattern);
}

bool OpenGLBuffer::bind()
{
	glBindBuffer(GL_ARRAY_BUFFER, mId);
	return true; // Should add error handling
}

bool OpenGLBuffer::create()
{
	glGenBuffers(1, &mId);
	return true; // Should add error handling
}

void OpenGLBuffer::destroy()
{
	// OpenGL docs state: "glDeleteBuffers silently ignores 0's and names that do not correspond to existing
	// buffer objects.". Therefore I assume multiple deletes are safe. Not sure if glDeleteBuffers() can also
	// modify the buffer id, but I avoid resetting it to zero incase OpenGL has done something useful with it.
	glDeleteBuffers(1, &mId);
}

void OpenGLBuffer::setUsagePattern(OpenGLBuffer::UsagePattern value)
{
	mUsagePattern = value;
}

void OpenGLBuffer::write(int offset, const void *data, int count)
{
	glBufferSubData(GL_ARRAY_BUFFER, offset, count, data);
}
