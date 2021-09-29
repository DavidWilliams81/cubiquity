#version 330 core
layout(location = 0) in vec3 vertexPosition_modelspace;
layout(location = 1) in vec4 particlePosition;
layout(location = 2) in vec4 particleNormal;

// Values that stay constant for the whole mesh.
uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;

out vec4 modelPosition;
out vec4 worldPosition;
out vec4 voxelNormal;
 
void main()
{
	// Pass the model space position through to the fragment shader,
	// we use this to compute the per-fragment normal when flat-shading.
	modelPosition = vec4(vertexPosition_modelspace, 1);
	
	// Pass through the per-voxel normal.
	voxelNormal = particleNormal;
	
	worldPosition = modelMatrix * vec4(vertexPosition_modelspace * particlePosition.w + particlePosition.xyz,1);
	
	mat4 modelViewProjectionMatrix = projectionMatrix * viewMatrix * modelMatrix; // Could just pass in from host
	gl_Position = modelViewProjectionMatrix * vec4(vertexPosition_modelspace * particlePosition.w + particlePosition.xyz,1);
}