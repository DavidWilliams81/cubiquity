#ifndef CUBIQUITY_INSTANCING_DEMO_H
#define CUBIQUITY_INSTANCING_DEMO_H

#include "opengl_viewer.h"

#include "extraction.h"

class InstancingDemo : public OpenGLViewer
{
public:

	InstancingDemo(const std::string& filename)
		: OpenGLViewer(filename) {}

	void onInitialise() override;
	void onUpdate(float deltaTime) override;
	void onShutdown() override;

	void onKeyDown(const SDL_KeyboardEvent& event);

	void onVolumeModified() override;


	GLuint glyphProgram;

	GLuint modelMatrixID;
	GLuint viewMatrixID;
	GLuint projMatrixID;
	GLuint cameraPosID;

	GLuint materialsTexture;
	GLuint materialsTextureID;

	Cubiquity::Glyph* mGlyphs;
	bool mNeedsUpdate = true;

	GLuint VertexArrayID;

	GLuint mCubeVerticesBuffer;

	GLuint mPerInstanceDataBuffer;


	uint32_t mGlyphCount;
};

#endif // CUBIQUITY_INSTANCING_DEMO_H
