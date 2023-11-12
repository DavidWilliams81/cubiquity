#version 330 core
layout(location = 0) in vec3 vertexPositionModelSpace;
layout(location = 1) in vec4 glyphPositionAndSize;
layout(location = 2) in vec4 glyphNormalAndMaterial;

// Values that stay constant for the whole mesh.
uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;

out vec4 positionModelSpace;
out vec4 positionWorldSpace;

out float glyphSize;
out vec3  glyphNormal;
out float glyphMaterial;
out vec4  glyphCentreWorldSpace;
 
void main()
{
	glyphSize = glyphPositionAndSize.w;
	glyphNormal = glyphNormalAndMaterial.xyz;
	glyphMaterial = glyphNormalAndMaterial.w;
	glyphCentreWorldSpace = modelMatrix * vec4(glyphPositionAndSize.xyz, 1.0);

	// Pass the model space position through to the fragment shader,
	// we use this to compute the per-fragment normal when flat-shading.
	positionModelSpace = vec4(vertexPositionModelSpace * glyphSize, 1);

	positionWorldSpace = modelMatrix * vec4(positionModelSpace.xyz + glyphPositionAndSize.xyz,1);

	mat4 viewProjectionMatrix = projectionMatrix * viewMatrix;
	gl_Position = viewProjectionMatrix * positionWorldSpace;
}