#ifndef CUBIQUITY_OPENGL_VIEWER_H
#define CUBIQUITY_OPENGL_VIEWER_H

#include "viewer.h"

#include <glad/glad.h>

// See https://www.geeksforgeeks.org/error-handling-in-opengl/
GLenum glCheckError_(const char* file, int line);
#define glCheckError() glCheckError_(__FILE__, __LINE__)

class OpenGLViewer : public Viewer
{
public:

	OpenGLViewer(const std::string& filename)
		: Viewer(filename, WindowType::OpenGL) {}

	void onInitialise() override;
	void onShutdown() override;

	void drawScreenAlignedQuad();

protected:
	//Screen aligned quad
	unsigned int quadVAO;

	GLuint loadShader(GLenum shaderType, const char* filePath);
	GLuint loadProgram(const char* vertex_file_path, const char* fragment_file_path);
};

#endif // CUBIQUITY_OPENGL_VIEWER_H
