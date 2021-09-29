#ifndef SHADER_HPP
#define SHADER_HPP

#include <glad/glad.h>

class Program
{
public:
	void LoadShaders(const char * vertex_file_path, const char * fragment_file_path);

	GLuint programID;
};

#endif
