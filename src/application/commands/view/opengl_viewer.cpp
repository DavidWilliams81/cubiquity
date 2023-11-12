#include "opengl_viewer.h"

#include <fstream>
#include <iostream>

float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
	// positions   // texCoords
	-1.0f,  1.0f,  0.0f, 1.0f,
	-1.0f, -1.0f,  0.0f, 0.0f,
	 1.0f, -1.0f,  1.0f, 0.0f,

	-1.0f,  1.0f,  0.0f, 1.0f,
	 1.0f, -1.0f,  1.0f, 0.0f,
	 1.0f,  1.0f,  1.0f, 1.0f
};

// See https://www.geeksforgeeks.org/error-handling-in-opengl/
GLenum glCheckError_(const char* file, int line)
{
	GLenum errorCode;
	while ((errorCode = glGetError()) != GL_NO_ERROR)
	{
		std::string error;
		switch (errorCode)
		{
		case GL_INVALID_ENUM:                  error = "INVALID_ENUM"; break;
		case GL_INVALID_VALUE:                 error = "INVALID_VALUE"; break;
		case GL_INVALID_OPERATION:             error = "INVALID_OPERATION"; break;
		case GL_STACK_OVERFLOW:                error = "STACK_OVERFLOW"; break;
		case GL_STACK_UNDERFLOW:               error = "STACK_UNDERFLOW"; break;
		case GL_OUT_OF_MEMORY:                 error = "OUT_OF_MEMORY"; break;
		case GL_INVALID_FRAMEBUFFER_OPERATION: error = "INVALID_FRAMEBUFFER_OPERATION"; break;
		}
		std::cout << error << " | " << file << "(" << line << ")" << std::endl;
	}
	return errorCode;
}

void OpenGLViewer::onInitialise()
{
	Viewer::onInitialise();

	unsigned int quadVBO;
	glGenVertexArrays(1, &quadVAO);
	glGenBuffers(1, &quadVBO);
	glBindVertexArray(quadVAO);
	glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
}

void OpenGLViewer::onShutdown()
{
	Viewer::onShutdown();
}

// Alternative approaches: https://stackoverflow.com/q/2588875
void OpenGLViewer::drawScreenAlignedQuad()
{
	glBindVertexArray(quadVAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
}

GLuint OpenGLViewer::loadShader(GLenum shaderType, const std::string& preamble, const char* filePath)
{
	// Create shader
	std::cout << "Loading " << filePath << "..." << std::endl;
	GLuint shader = glCreateShader(shaderType);

	// Load contents of file
	std::ifstream fileStream(filePath);
	if (!fileStream.is_open())
	{
		std::cerr << "Failed to open \"" << filePath << "\"" << std::endl;
	}
	std::string shaderCode(
		(std::istreambuf_iterator<char>(fileStream)),
		(std::istreambuf_iterator<char>()));

	// GLSL requires that #version (if present) is the first thing in
	// the file. Therefore any preamble has to be inserted after it.
	// Apparently this is the only way to #define preprocessor 
	// macros from C code: https://stackoverflow.com/a/19541536
	size_t preambleInsertionIndex = 0;
	auto versionIndex = shaderCode.find("#version");
	if (versionIndex != std::string::npos)
	{
		auto newlineIndex = shaderCode.find('\n', versionIndex);
		preambleInsertionIndex = newlineIndex + 1;
	}

	// Insert the preamble (and a newline, in case the user forgot).
	shaderCode.insert(preambleInsertionIndex, preamble + "\n");

	// Compile shader
	char const* shaderCodePtr = shaderCode.c_str();
	glShaderSource(shader, 1, &shaderCodePtr, NULL);
	glCompileShader(shader);

	// Check for errors
	int infoLogLength;
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
	if (infoLogLength > 0)
	{
		std::vector<char> infoLog(infoLogLength);
		glGetShaderInfoLog(shader, infoLogLength, NULL, infoLog.data());
		std::cerr << infoLog.data() << std::endl;
	}

	return shader;
}

GLuint OpenGLViewer::loadProgram(const char* vertex_file_path, const char* fragment_file_path, std::string preamble)
{
	// Create the shaders
	GLuint vertexShader = loadShader(GL_VERTEX_SHADER, preamble, vertex_file_path);
	GLuint fragmentShader = loadShader(GL_FRAGMENT_SHADER, preamble, fragment_file_path);

	// Link the program
	std::cout << "Linking program" << std::endl;
	GLuint program = glCreateProgram();
	glAttachShader(program, vertexShader);
	glAttachShader(program, fragmentShader);
	glLinkProgram(program);

	// Check the program
	int infoLogLength;
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);
	if (infoLogLength > 0)
	{
		std::vector<char> infoLog(infoLogLength);
		glGetProgramInfoLog(program, infoLogLength, NULL, infoLog.data());
		std::cerr << infoLog.data() << std::endl;
	}

	glDetachShader(program, vertexShader);
	glDetachShader(program, fragmentShader);

	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	return program;
}