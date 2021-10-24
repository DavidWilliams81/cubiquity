#version 330 core

in vec4 modelPosition;
in vec4 worldPosition;
in vec4 voxelNormal;   // Per-voxel normal, as computed in C++ code

out vec4 color;

uniform sampler1D materials;

////////////////////////////////////////////////////////////////////////////////
// Functions for calculating normals based on model/world positions.
////////////////////////////////////////////////////////////////////////////////

// The model-space position of the fragment represent a crude approximation to the normal.
// We can perfect it by snapping it to the nearest axis as our cubes are always axis-aligned.
// The normal should never need flipping (unlike derivative version) but does suffer from
// artefacts at the edges of the cubes. These may be hidden by bevelling.
// Note: The input fragments are assumed to be from a cube covering the range -0.5 to +0.5.
vec3 flatNormalBySnapping(vec3 modelPosition)
{
	vec3 normal = modelPosition; // Initial (crude) approximation to the normal.
	vec3 signs = sign(normal); // Save these as the abs() operation will lose them
	vec3 absFaceNormal = abs(normal); // One component will be 0.5, other will be smaller.
	normal = step(0.49999f, normal); // Set single component to 1.0, other two to 0.0.
	normal = normal * signs; // Put the signs back
	return normal; // Done!
}

// This normal calculation method is generally higher quality than the model-space position
// 'snapping' method, but may have higher hardware requirements and in Unity we have seen
// issues with the normal sometimes pointing the wrong way (based on render system , edit vs.
// play mode, etc). These might be fixed by flipping the normal if it points away from the
// camera, as we never see away-facing fragments anyway?
vec3 normalFromDerivatives(vec3 position)
{
	// dFdx(position) and dFdy(position) give two vectors lying on our surface.
	// The cross product gives a vector which is orthogonal to both of them.
	vec3 normal = cross(dFdx(position), dFdy(position));
	normal = normalize(normal);	
	return normal;
}

// Raising the normals to a power distorts the distribution and pushes them
// towards the centre of the face, creating a rounded-corners effect.	
// Note: The input fragments are assumed to be from a cube covering the range -0.5 to +0.5.
vec3 roundedNormalForGlyph(vec3 modelPosition)
{	
	vec3 normal = modelPosition;
	// Do an odd number of multiplications to preserve the sign. Otherwise we
	// have to store the sign separately and then multiply by it afterwards.
	normal = normal * normal * normal * normal * normal;	
	normal = normalize(normal);	
	return normal;
}

vec3 roundedNormalForVoxel(vec3 modelPosition, vec3 worldPosition)
{
	// Take the voxel-space position as the fractional part of the world position.
	vec3 voxelPosition = fract(worldPosition);
	
	// Our position lies on the surface of a voxel, but we don't know whether to round
	// up or down to get inside the voxel. Therefore we use an approximation of the
	// normal to generate an offset, which we apply before doing the rounding operation.
	vec3 offset = modelPosition * 0.01f;	
	vec3 centre = round(voxelPosition - offset);	
	vec3 normal = voxelPosition - centre;
	
	// Similar logic to 'roundedNormalForGlyph()'.
	normal = normal * normal * normal * normal * normal;	
	normal = normalize(normal);	
	return normal;
}

float positionBasedNoise(vec4 positionAndStrength)
{
	//'floor' is more widely supported than 'round'. Offset consists of:
	//  - 0.5 to perform the rounding
	vec3 roundedPos = floor(positionAndStrength.xyz + vec3(0.5, 0.5, 0.5));
	
	// The noise function below seems to do remarkably well even for large inputs. None-the-less, it is
	// better for small inputs so we limit the input range to within a few hundred voxels of the origin.
	roundedPos += vec3(256.0, 256.0, 256.0);
	roundedPos += mod(roundedPos, vec3(512.0, 512.0, 512.0));
	
	// Based on this: http://byteblacksmith.com/improvements-to-the-canonical-one-liner-glsl-rand-for-opengl-es-2-0/
	// I've modified it to use a 3D seed with the third coefficient being a number I picked at random. Because it is 
	// using a 3D seed the magnitude of the dot product could be much larger, so I've reduced each of the coefficients
	// by a factor of 10 to limit precision problems for high seed values. We can tweak these further in the future.	
	float noise = fract(sin(mod(dot(roundedPos.xyz, vec3(1.29898,7.8233, 4.26546)), 3.14)) * 43758.5453);
	
	//Scale the noise
	float halfNoiseStrength = positionAndStrength.w * 0.5;
	noise = -halfNoiseStrength + positionAndStrength.w * noise; //http://www.opengl.org/wiki/GLSL_Optimizations#Get_MAD
	
	return noise;
}

void main()
{	
	// For visualization we derive the color of the voxel from its normal.
	//vec4 voxelColor = vec4(voxelNormal.xyz * 0.5 + vec3(0.5, 0.5, 0.5), 1.0);
	//vec4 voxelColor = vec4(decodeMaterial(voxelNormal.w), 1.0);

	int materialCount = textureSize(materials, 0);
	float u = (voxelNormal.w+0.5) / materialCount; // Map material id to range 0.0 - 1.0
	vec4 voxelColor = texture(materials, u);
	
	float noise = positionBasedNoise(vec4(worldPosition.xyz, 0.08));
	
	voxelColor.xyz += noise;

	vec3 worldNormal = normalFromDerivatives(modelPosition.xyz);
	
	//fragColor = vec4((worldNormal * 0.5f) + 0.5f, 1.0f);
	
	
	//vec3 ambientIntensity = vec3(0.3, 0.3, 0.3);
    //vec3 diffuseIntensity;
    //adModel(worldPosition, worldNormal, diffuseIntensity);
	//vec3 totalLightIntensity = ambientIntensity + diffuseIntensity;
    //fragColor = vec4( totalLightIntensity * voxelColor.xyz, 1.0 );
	

	
	// Lighting properties
	vec3 lightDir = normalize(vec3(-0.3, -0.5, -0.7));
	float ka = 0.3;
	float kd = 0.7;
	
	float lightIntensity = 0.7 * max(0.0, dot(worldNormal, -lightDir)) + 0.3;	
	
	color = voxelColor * lightIntensity;	
	color.a = 1.0f;
}