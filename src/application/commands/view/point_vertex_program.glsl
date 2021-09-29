#version 330 core
layout(location = 0) in vec3 vertexPosition; // Relative to the centre of the glyph
layout(location = 1) in vec4 glyphPositionAndSize;
layout(location = 2) in vec4 glyphNormal;

// Values that stay constant for the whole mesh.
uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;

out vec4 modelSpacePosition;
out vec4 worldSpacePosition;
out vec4 worldSpaceCentre;
out vec4 voxelNormal;
 
void main()
{
	// Should be set just big enough so that the points don't clip against the edge of the glyph.
	float glyphScaleFactor = 1.5f;
	
	// Pass the model space position through to the fragment shader,
	// we use this to compute the per-fragment normal when flat-shading.
	modelSpacePosition = vec4(vertexPosition * glyphPositionAndSize.w * glyphScaleFactor + glyphPositionAndSize.xyz,1);
	
	worldSpacePosition = modelMatrix * modelSpacePosition;
	
	worldSpaceCentre = modelMatrix * vec4(glyphPositionAndSize.xyz, 1);
	
	
	
	// Pass through the per-voxel normal.
	voxelNormal = glyphNormal;
	
	mat4 modelViewProjectionMatrix = projectionMatrix * viewMatrix * modelMatrix; // Could just pass in from host
	gl_Position = modelViewProjectionMatrix * modelSpacePosition;
}