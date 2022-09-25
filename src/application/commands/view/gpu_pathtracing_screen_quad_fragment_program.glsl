#version 460 core

////////////////////////////////////////////////////////////////////////////////
// Utility code shared between instancing and pathtracing shaders
////////////////////////////////////////////////////////////////////////////////

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

vec4 getMaterial(float materialId, sampler1D materialsTexture)
{
	int materialCount = textureSize(materialsTexture, 0);
	float u = (materialId+0.5) / materialCount; // Map material id to range 0.0 - 1.0
	return texture(materialsTexture, u);
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

////////////////////////////////////////////////////////////////////////////////
// End of shared utility code
////////////////////////////////////////////////////////////////////////////////

struct RayVolumeIntersection
{
	vec3 position;
	vec3 normal;
	float distance;
	uint material;
};

struct SubDAG
{
	uint nodeIndex;
	uint padding1, padding2, padding3;
	ivec3 lowerBound;
	int padding4;
	ivec3 upperBound; // Could store size instead?
	int padding5;
};

// FIXME - Rather than hard-coding these bindings I think we should be able to retrieve them through
// glGetProgramResourceIndex() and/or glShaderStorageBlockBinding(). Unfortunately that is currently
// crashing for me inside the driver. Might be worth testing again after an upgrade/update. See
// https://www.khronos.org/opengl/wiki/Interface_Block_(GLSL)#OpenGL_binding_index_setting
// https://web.archive.org/web/20150109112923/http://www.lighthouse3d.com/tutorials/glsl-core-tutorial/3490-2/
layout(binding = 0, std430) buffer SomeName  
{
	uint dagData[][8];  
};

layout(binding = 1, std430) buffer SomeOtherName  
{
	SubDAG subDAGs[8];  
};

struct Ray3d
{
	vec3 mOrigin;
	vec3 mDir;
};

const uint MaterialCount = 255;



uint a;

uint firstNode(float tx0, float ty0, float tz0, float txm, float tym, float tzm)
{
	uint answer = 0u;   // initialize to 00000000
								// select the entry plane and set bits
	if (tx0 > ty0)
	{
		if (tx0 > tz0) // PLANE YZ
		{
			if (tym < tx0) answer |= (1 << 1);
			if (tzm < tx0) answer |= (1 << 2);
			return answer;
		}
	}
	else {
		if (ty0 > tz0) // PLANE XZ
		{
			if (txm < ty0) answer |= (1 << 0);
			if (tzm < ty0) answer |= (1 << 2);
			return answer;
		}
	}
	// PLANE XY
	if (txm < tz0) answer |= (1 << 0);
	if (tym < tz0) answer |= (1 << 1);
	return answer;
}

uint nextNode(float txm, uint x, float tym, uint y, float tzm, uint z)
{
	if (txm < tym)
	{
		if (txm < tzm) { return x; }  // YZ plane
	}
	else
	{
		if (tym < tzm) { return y; } // XZ plane
	}
	return z; // XY plane;
}

uint getNode(uint node, uint childId)
{
	return dagData[node][childId];
}

void intersectSubtreeIterative(float tx0In, float ty0In, float tz0In, float tx1In, float ty1In, float tz1In, uint nodeIndexIn, out RayVolumeIntersection intersection)
{
	const uint INVALID_CHILD = 9u; // Only eight children

	struct State
	{
		float tx0, ty0, tz0, tx1, ty1, tz1, txm, tym, tzm;
		uint nodeIndex;
		uint currNode;
	};

	int depth = 0; // Relative to subtree (not root)
	State states[33]; // FIXME - How big should this be?
	//State* pState = &(stack[0]);
	//memset(pState, 0, sizeof(State) * 33);
	State state = { tx0In, ty0In, tz0In, tx1In, ty1In, tz1In, 0.0, 0.0, 0.0, nodeIndexIn, INVALID_CHILD };
	states[depth] = state;
	
	

	while(true)
	{
		bool skip = false;
		if (depth < 0)
		{
			break;
		}

		if (states[depth].currNode == INVALID_CHILD)
		{
			//std::cout << indent(level) << state.childId << std::endl;

			if (states[depth].tx1 < 0.0 || states[depth].ty1 < 0.0 || states[depth].tz1 < 0.0)
			{
				depth--;
				skip = true;
			}

			if (skip == false)
			{
				if (states[depth].nodeIndex < MaterialCount) // is a material node
				{
					if (states[depth].nodeIndex > 0) // Occupied node
					{
						intersection.material = states[depth].nodeIndex;
						intersection.distance = max(max(states[depth].tx0, states[depth].ty0), states[depth].tz0);
						intersection.normal[0] = 0.0;
						intersection.normal[1] = 0.0;
						intersection.normal[2] = 0.0;
						if (states[depth].tx0 > states[depth].ty0 && states[depth].tx0 > states[depth].tz0)
						{
							intersection.normal[0] = -1.0;
						}
						if (states[depth].ty0 > states[depth].tx0 && states[depth].ty0 > states[depth].tz0)
						{
							intersection.normal[1] = -1.0;
						}
						if (states[depth].tz0 > states[depth].tx0 && states[depth].tz0 > states[depth].ty0)
						{
							intersection.normal[2] = -1.0;
						}

						// Flip normals if required
						if (bool(a & 1)) { intersection.normal[0] *= -1.0; }
						if (bool(a & 2)) { intersection.normal[1] *= -1.0; }
						if (bool(a & 4)) { intersection.normal[2] *= -1.0; }

						return;
					}
					depth--;
					skip = true;
				}
			}

			if (skip == false)
			{
				// FIXME - We need to handle infinite values here. Either as described in the paper, or
				// by adding a tiny offset to input vector components to make sure they are never zero.
				states[depth].txm = 0.5 * (states[depth].tx0 + states[depth].tx1);
				states[depth].tym = 0.5 * (states[depth].ty0 + states[depth].ty1);
				states[depth].tzm = 0.5 * (states[depth].tz0 + states[depth].tz1);

				states[depth].currNode = firstNode(states[depth].tx0, states[depth].ty0, states[depth].tz0, states[depth].txm, states[depth].tym, states[depth].tzm);
			}

		}

		// Note: The constants below are in the reverse order compared to the paper. Cubiquity uses the
		// LSBs in 'zyx' to index child nodes, but the paper *appears* to use 'xyz'? The paper actually
		// seems inconsistent, because Figure 1 imples 'xyz' order but Table 1 implies 'zyx'? Or they
		// are just numbering their bits differently? Maybe I am misunderstanding something.
		//State* pNextState = pState + 1;
		if (skip == false)
		{
			switch (states[depth].currNode)
			{
				case 0:
				{
					State nextState = { states[depth].tx0, states[depth].ty0, states[depth].tz0, states[depth].txm, states[depth].tym, states[depth].tzm, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,a), INVALID_CHILD };;
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].txm, 1u, states[depth].tym, 2u, states[depth].tzm, 4u);
					depth++;
					break;
				}
				case 1:
				{
					State nextState = { states[depth].txm, states[depth].ty0, states[depth].tz0, states[depth].tx1, states[depth].tym, states[depth].tzm, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,1 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].tx1, 8u, states[depth].tym, 3u, states[depth].tzm, 5u);
					depth++;
					break;
				}
				case 2:
				{
					State nextState = { states[depth].tx0, states[depth].tym, states[depth].tz0, states[depth].txm, states[depth].ty1, states[depth].tzm, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,2 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].txm, 3u, states[depth].ty1, 8u, states[depth].tzm, 6u);
					depth++;
					break;
				}
				case 3:
				{
					State nextState = { states[depth].txm, states[depth].tym, states[depth].tz0, states[depth].tx1, states[depth].ty1, states[depth].tzm, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,3 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].tx1, 8u, states[depth].ty1, 8u, states[depth].tzm, 7u);
					depth++;
					break;
				}
				case 4:
				{
					State nextState = { states[depth].tx0, states[depth].ty0, states[depth].tzm, states[depth].txm, states[depth].tym, states[depth].tz1, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,4 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].txm, 5u, states[depth].tym, 6u, states[depth].tz1, 8u);
					depth++;
					break;
				}
				case 5:
				{
					State nextState = { states[depth].txm, states[depth].ty0, states[depth].tzm, states[depth].tx1, states[depth].tym, states[depth].tz1, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,5 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].tx1, 8u, states[depth].tym, 7u, states[depth].tz1, 8u);
					depth++;
					break;
				}
				case 6:
				{
					State nextState = { states[depth].tx0, states[depth].tym, states[depth].tzm, states[depth].txm, states[depth].ty1, states[depth].tz1, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,6 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = nextNode(states[depth].txm, 7u, states[depth].ty1, 8u, states[depth].tz1, 8u);
					depth++;
					break;
				}
				case 7:
				{
					State nextState = { states[depth].txm, states[depth].tym, states[depth].tzm, states[depth].tx1, states[depth].ty1, states[depth].tz1, 0.0, 0.0, 0.0, getNode(states[depth].nodeIndex,7 ^ a), INVALID_CHILD };
					states[depth + 1] = nextState;
					states[depth].currNode = 8;
					depth++;
					break;
				}
				case 8:
				{
					depth--;
				}
			}
		}
	}
}

RayVolumeIntersection intersectNodes(Ray3d ray)
{
	RayVolumeIntersection intersection;
	intersection.material = 0;

	int minInt = -2147483648;
	ivec3 rootLowerBound = { minInt , minInt , minInt };

	// Near to far octree traversal
	uint nearestChild = 0;
	if (ray.mOrigin[0] > -0.5) nearestChild |= 0x1;
	if (ray.mOrigin[1] > -0.5) nearestChild |= 0x2;
	if (ray.mOrigin[2] > -0.5) nearestChild |= 0x4;
	
	// Used for near to far octree traversal described here:
	// https://www.flipcode.com/archives/Harmless_Algorithms-Issue_02_Scene_Traversal_Algorithms.shtml#octh
	// Note: '4' comes before '3' in this array. This is not a mistake (see link above)
	const uint nearToFar[] = { 0x00, 0x01, 0x02, 0x04, 0x03, 0x05, 0x06, 0x07 };

	for (uint i = 0; i < 8; i++) 
	{
		uint childId = nearestChild ^ nearToFar[i];
		uint childNodeIndex = subDAGs[childId].nodeIndex;

		Ray3d reflectedRay = ray;
		ivec3 reflectedLowerBound = subDAGs[childId].lowerBound;
		ivec3 reflectedUpperBound = subDAGs[childId].upperBound;

		a = 0;
		if (reflectedRay.mDir[0] < 0.0)
		{
			reflectedRay.mOrigin[0] = -reflectedRay.mOrigin[0] - 1.0;
			reflectedRay.mDir[0] = -reflectedRay.mDir[0];
			reflectedLowerBound[0] = ~reflectedLowerBound[0];
			reflectedUpperBound[0] = ~reflectedUpperBound[0];
			int temp = reflectedLowerBound[0]; // Replace below with std::swap.
			reflectedLowerBound[0] = reflectedUpperBound[0];
			reflectedUpperBound[0] = temp;
			a |= 1;
		}
		if (reflectedRay.mDir[1] < 0.0)
		{
			reflectedRay.mOrigin[1] = -reflectedRay.mOrigin[1] - 1.0;
			reflectedRay.mDir[1] = -reflectedRay.mDir[1];
			reflectedLowerBound[1] = ~reflectedLowerBound[1];
			reflectedUpperBound[1] = ~reflectedUpperBound[1];
			int temp = reflectedLowerBound[1]; // Replace below with std::swap.
			reflectedLowerBound[1] = reflectedUpperBound[1];
			reflectedUpperBound[1] = temp;
			a |= 2;
		}
		if (reflectedRay.mDir[2] < 0.0)
		{
			reflectedRay.mOrigin[2] = -reflectedRay.mOrigin[2] - 1.0;
			reflectedRay.mDir[2] = -reflectedRay.mDir[2];
			reflectedLowerBound[2] = ~reflectedLowerBound[2];
			reflectedUpperBound[2] = ~reflectedUpperBound[2];
			int temp = reflectedLowerBound[2]; // Replace below with std::swap.
			reflectedLowerBound[2] = reflectedUpperBound[2];
			reflectedUpperBound[2] = temp;
			a |= 4;
		}

		float tx0 = ((reflectedLowerBound[0] - 0.5) - reflectedRay.mOrigin[0]) / reflectedRay.mDir[0];
		float tx1 = ((reflectedUpperBound[0] + 0.5) - reflectedRay.mOrigin[0]) / reflectedRay.mDir[0];
		float ty0 = ((reflectedLowerBound[1] - 0.5) - reflectedRay.mOrigin[1]) / reflectedRay.mDir[1];
		float ty1 = ((reflectedUpperBound[1] + 0.5) - reflectedRay.mOrigin[1]) / reflectedRay.mDir[1];
		float tz0 = ((reflectedLowerBound[2] - 0.5) - reflectedRay.mOrigin[2]) / reflectedRay.mDir[2];
		float tz1 = ((reflectedUpperBound[2] + 0.5) - reflectedRay.mOrigin[2]) / reflectedRay.mDir[2];

		if (max(max(tx0, ty0), tz0) < min(min(tx1, ty1), tz1))
		{
			intersectSubtreeIterative(tx0, ty0, tz0, tx1, ty1, tz1, childNodeIndex, intersection);
		}

		if (intersection.material > 0) { break; }
	}

	intersection.position = ray.mOrigin + (ray.mDir * intersection.distance);
	return intersection;
}

// This assumes the pixel position px to be in [0,1],
// which can be done by (x+0.5)/w or (y+0.5)/h (or h-y +0.5 for screens
// with top left origin) to sample pixel centers
// See https://sibaku.github.io/computer-graphics/2017/01/10/Camera-Ray-Generation.html
vec3 createRay(vec2 px, mat4 PInv, mat4 VInv)
{
  
    // convert pixel to NDS
    // [0,1] -> [-1,1]
    vec2 pxNDS = px*2. - 1.;

    // choose an arbitrary point in the viewing volume
    // z = -1 equals a point on the near plane, i.e. the screen
    vec3 pointNDS = vec3(pxNDS, -1.);

    // as this is in homogenous space, add the last homogenous coordinate
    vec4 pointNDSH = vec4(pointNDS, 1.0);
    // transform by inverse projection to get the point in view space
    vec4 dirEye = PInv * pointNDSH;

    // since the camera is at the origin in view space by definition,
    // the current point is already the correct direction 
    // (dir(0,P) = P - 0 = P as a direction, an infinite point,
    // the homogenous component becomes 0 the scaling done by the 
    // w-division is not of interest, as the direction in xyz will 
    // stay the same and we can just normalize it later
    dirEye.w = 0.;

    // compute world ray direction by multiplying the inverse view matrix
    vec3 dirWorld = (VInv * dirEye).xyz;

    // now normalize direction
    return normalize(dirWorld); 
}

uniform mat4 VInv;
uniform mat4 PInv;
uniform vec3 cameraPos;
uniform sampler1D materials;

out vec4 FragColor;

in vec2 TexCoords;

void main()
{
	a = 0;
	vec3 rayDir = createRay(TexCoords, PInv, VInv);
	
	Ray3d ray;
	ray.mOrigin = cameraPos;
	ray.mDir = rayDir;
	
	RayVolumeIntersection rayVolumeIntersection = intersectNodes(ray);
	
	if(rayVolumeIntersection.material > 0)
	{
		vec4 voxelColor = getMaterial(rayVolumeIntersection.material, materials);
		
		FragColor.rgb = light(voxelColor.xyz, rayVolumeIntersection.position.xyz, rayVolumeIntersection.normal, cameraPos);
		
		FragColor.a = 1.0;
	}
	else
	{
		FragColor = vec4(0.0, 0.0, 0.2, 1.0);
	}
} 