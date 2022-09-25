#version 330 core

in vec4 positionModelSpace;
in vec4 positionWorldSpace;

in float glyphSize;
in vec3  glyphNormal;
in float glyphMaterial;
in vec4  glyphCentreWorldSpace;

out vec4 color;

uniform vec3 cameraPos;
uniform sampler1D materials;
uniform unsigned int mode;

const uint CUBE = 0u;
const uint DISC = 1u;

////////////////////////////////////////////////////////////////////////////////
// Functions for calculating normals based on model/world positions.
////////////////////////////////////////////////////////////////////////////////

// The model-space position of the fragment represent a crude approximation to the normal.
// We can perfect it by snapping it to the nearest axis as our cubes are always axis-aligned.
// The normal should never need flipping (unlike derivative version) but does suffer from
// artefacts at the edges of the cubes. These may be hidden by bevelling.
// Note: The input fragments are assumed to be from a cube covering the range -0.5 to +0.5.
vec3 flatNormalBySnapping(vec3 positionModelSpace)
{
	vec3 normal = positionModelSpace; // Initial (crude) approximation to the normal.
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
vec3 roundedNormalForGlyph(vec3 positionModelSpace)
{
	vec3 normal = positionModelSpace;
	// Do an odd number of multiplications to preserve the sign. Otherwise we
	// have to store the sign separately and then multiply by it afterwards.
	normal = normal * normal * normal * normal * normal;
	normal = normalize(normal);	
	return normal;
}

vec3 roundedNormalForVoxel(vec3 positionModelSpace, vec3 positionWorldSpace)
{
	// Take the voxel-space position as the fractional part of the world position.
	vec3 voxelPosition = fract(positionWorldSpace);

	// Our position lies on the surface of a voxel, but we don't know whether to round
	// up or down to get inside the voxel. Therefore we use an approximation of the
	// normal to generate an offset, which we apply before doing the rounding operation.
	vec3 offset = positionModelSpace * 0.01f;
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
	// - 0.5 to perform the rounding
	// - Tiny offset to eliminate acne as surfaces lie *exactly* between two voxels.
	vec3 roundedPos = floor(positionAndStrength.xyz + vec3(0.501, 0.501, 0.501));

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

float min3 (vec3 v) {
  return min (min (v.x, v.y), v.z);
}

float max3 (vec3 v) {
  return max (max (v.x, v.y), v.z);
}

// See https://tavianator.com/fast-branchless-raybounding-box-intersections/
vec2 intersect(vec3 rayOrigin, vec3 rayDir, vec3 boxLower, vec3 boxUpper)
{
	// Inverse direction could be precomputed and stored in the ray
	// if we find we often intersect the same ray with multiple boxes.
	vec3 invDir = vec3(1.0, 1.0, 1.0) / rayDir;

	vec3 lower = (boxLower - rayOrigin) * invDir;
	vec3 upper = (boxUpper - rayOrigin) * invDir;

	vec3 minCorner = min(lower, upper);
	vec3 maxCorner = max(lower, upper);

	vec2 intersection;
	intersection.x = max3(minCorner);
	intersection.y = min3(maxCorner);
	return intersection;
}

bool intersectPlane(vec3 n, vec3 p0, vec3 l0, vec3 l, inout float t) 
{ 
    // assuming vectors are all normalized
    float denom = dot(n, l); 
    if (denom > 1e-6) { 
        vec3 p0l0 = p0 - l0; 
        t = dot(p0l0, n) / denom; 
        return (t >= 0); 
    } 
 
    return false; 
}

// Simple local lighting model similar to OpenGL fixed-function pipeline.
// See: https://developer.download.nvidia.com/CgTutorial/cg_tutorial_chapter05.html
vec3 light(vec3 surfaceColour, vec3 P, vec3 N, vec3 eyePosition)
{
	// Material properties
	vec3 Ke = vec3(0.0, 0.0, 0.0);
	vec3 Ka = surfaceColour;
	vec3 Kd = surfaceColour;
	vec3 Ks = vec3(0.0, 0.0, 0.0);
	float shininess = 1.0f;

	// Light properties
	vec3 lightDir = vec3(0.3, 0.5, 0.7); // Pointing towards light
	vec3 lightColor = vec3(0.7, 0.7, 0.7);
	vec3 globalAmbient = vec3(0.3, 0.3, 0.3);

	// Compute emissive contribution
	vec3 emissive = Ke;

	// Compute ambient contribution
	vec3 ambient = Ka * globalAmbient;

	// Compute diffuse contribution
	vec3 L = normalize(lightDir);
	float diffuseLight = max(dot(N, L), 0);
	vec3 diffuse = Kd * lightColor * diffuseLight;

	// Compute specular contribution
	vec3 V = normalize(eyePosition - P);
	vec3 H = normalize(L + V);
	float specularLight = pow(max(dot(N, H), 0), shininess);

	if (diffuseLight <= 0) specularLight = 0;
	vec3 specular = Ks * lightColor * specularLight;

	vec3 result = emissive + ambient + diffuse + specular;
	return result;
}

void main()
{
	if(mode == CUBE)
	{
		int materialCount = textureSize(materials, 0);
		float u = (glyphMaterial+0.5) / materialCount; // Map material id to range 0.0 - 1.0
		vec4 voxelColor = texture(materials, u);

		float noise = positionBasedNoise(vec4(positionWorldSpace.xyz, 0.08));

		voxelColor.xyz += noise;

		// The glyph normal can be zero. This might be because glyph normals weren't generated, or
		// because this one is difficult (for example, a single voxel in an otherwise empty volume).
		// As a fallback for cubic glyphs we then generate a normal from the cube face.
		vec3 worldNormal;
		if(length(glyphNormal) > 0.1)
		{
			worldNormal = glyphNormal;
		}
		else
		{
			worldNormal = normalFromDerivatives(positionModelSpace.xyz);
		}
		
		color.rgb = light(voxelColor.xyz, positionModelSpace.xyz, worldNormal, cameraPos);
		//color.rgb = glyphNormal.xyz * 0.5 + vec3(0.5, 0.5, 0.5);
		//color.rgb = normalize((worldNormal.xyz)) * 0.5 + vec3(0.5, 0.5, 0.5);
		color.a = 1.0;
	}
	
	if(mode == DISC)
	{
		vec3 rayDir = normalize(positionWorldSpace.xyz - cameraPos);
		float distToSurface = length(rayDir);
		
		float t = 0.0f;
		bool hit = intersectPlane(rayDir, glyphCentreWorldSpace.xyz, cameraPos, rayDir, t);
		
		vec3 intersection = cameraPos + (rayDir * t);
		
		vec3 temp = intersection.xyz - glyphCentreWorldSpace.xyz;
		float dist = length(temp);
		const float overlap = 1.75;
		float discRadius = glyphSize * 0.5 * overlap; // 0.5 converts diameter to radius, 1.75 factor makes discs overlap
		float weight = max(discRadius - dist, 0.0);
		weight /= discRadius; // To go to 0.0 - 1.0
		
		// The glyph normal can be zero. This might be because glyph normals weren't generated, or
		// because this one is difficult (for example, a single voxel in an otherwise empty volume).
		// As a fallback for disc glyphs we then generate a normal which is sort of spherical, by
		// interpolating between a normal pointing at the camera and one which is orthogonal to it.
		// There are probably better approaches, this is rather experimental.
		vec3 worldNormal;
		if(length(glyphNormal) > 0.1)
		{
			worldNormal = glyphNormal;
		}
		else
		{
			vec3 centreNormal = normalize(-rayDir);
			vec3 edgeNormal = normalize(temp);
			worldNormal = centreNormal * weight + edgeNormal * (1.0 - weight);
			worldNormal = normalize(worldNormal);
		}
		
		int materialCount = textureSize(materials, 0);
		float u = (glyphMaterial+0.5) / materialCount; // Map material id to range 0.0 - 1.0
		vec4 voxelColor = texture(materials, u);

		color.rgb = light(voxelColor.xyz, intersection, worldNormal, cameraPos) * weight;
		//color.rgb = (normalize((worldNormal.xyz)) * 0.5 + vec3(0.5, 0.5, 0.5)) * weight;
		color.a = weight;
	}
}