#ifndef OPENGLBUFFER_H_DB940F68
#define OPENGLBUFFER_H_DB940F68

#include <glad/glad.h>

// Simple wrapper for OpenGL buffers, based on the API design in Qt5.
class OpenGLBuffer
{
public:

	enum UsagePattern // Enum rather than enum class simply to follow Qt model.
	{
		StreamDraw = GL_STREAM_DRAW,
		StaticDraw = GL_STATIC_DRAW
	};

	OpenGLBuffer();
	~OpenGLBuffer();

	void allocate(int count);
	void allocate(const void *data, int count);
	bool bind();
	bool create();
	void destroy();
	void setUsagePattern(OpenGLBuffer::UsagePattern value);
	void write(int offset, const void *data, int count);

private:

	GLuint mId;
	UsagePattern mUsagePattern;
};

#endif //OPENGLBUFFER_H_DB940F68