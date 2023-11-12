#ifndef CUBIQUITY_INSTANCING_DEMO_H
#define CUBIQUITY_INSTANCING_DEMO_H

#include "opengl_viewer.h"

#include "visibility.h"

class InstancingDemo : public OpenGLViewer
{
public:

	InstancingDemo(const std::string& filename)
		: OpenGLViewer(filename) {}

	void onInitialise() override;
	void onUpdate(float deltaTime) override;
	void onShutdown() override;

	void onKeyDown(const SDL_KeyboardEvent& event);


	GLuint glyphProgram;

	GLuint modelMatrixID;
	GLuint viewMatrixID;
	GLuint projMatrixID;
	GLuint cameraPosID;

	GLuint materialsTexture;
	GLuint materialsTextureID;

	Cubiquity::VisibilityCalculator* mVisibilityCalculator = nullptr;

	Cubiquity::Glyph* mGlyphs;
	bool mDoGlyphUpdates = true;

	GLuint VertexArrayID;

	GLuint mCubeVerticesBuffer;

	GLuint mPerInstanceDataBuffer;


	uint32_t mGlyphCount;
};

#endif // CUBIQUITY_INSTANCING_DEMO_H
