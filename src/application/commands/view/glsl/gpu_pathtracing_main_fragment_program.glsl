#version 460 core

////////////////////////////////////////////////////////////////////////////////
// Utility code shared between instancing and pathtracing shaders
////////////////////////////////////////////////////////////////////////////////

// Dummy function to aid writing code that compiles as both GLSL and C/C++
void assert(bool unused) {}

float min3(vec3 vec)
{
	return min(min(vec.x, vec.y), vec.z);
}

float max3(vec3 vec)
{
	return max(max(vec.x, vec.y), vec.z);
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
	uint packedNodeAndChild;
};

struct SubDAG
{
	uint nodeIndex;
	uint depth, padding2, padding3;
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

const uint MaterialCount = 256;
const float gMaxFootprint = 0.005;



uint a;

/*uint firstNode(float tx0, float ty0, float tz0, float txm, float tym, float tzm)
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
}*/

// I seem to have arrived at a slightly simpler construction than in the ESVO paper?
// They have separate logic for first node vs next node, but I can use the logic below
// for both. Presumably I (or they) do something slightly different elsewhere.
// Also note we could use bvec type and 'lessThanEqual' GLSL function which might
// be better than the below, but then we would need to convert the bvec to/from a
// bitfield to later use it for indexing a node. Probabbly not optimal but could be
// tested. Keep in mind this: https://stackoverflow.com/a/71857863
uint lessThanOrEqualTo(vec3 vec, float threshold)
{
	uint result = 0;
	if (vec.x <= threshold) { result |= 0x1; }
	if (vec.y <= threshold) { result |= 0x2; }
	if (vec.z <= threshold) { result |= 0x4; }
	return result;
}

uint firstNode(vec3 t0, vec3 tm)
	{
		// Distance to exit cube. I believe the efficient sparse voxel octree' paper
		// tracks the cube size explicitly, but I don't think we need to?
		float tcMin = max(t0.x, max(t0.y, t0.z));
		uint result = 0;

		if (tm.x <= tcMin) { result |= 0x1; }
		if (tm.y <= tcMin) { result |= 0x2; }
		if (tm.z <= tcMin) { result |= 0x4; }

		return result;
	}

/*uint nextNode(vec3 tm, uvec3 bits)
	{
		if (tm.x < tm.y)
		{
			if (tm.x < tm.z) { return bits.x; }  // YZ plane
		}
		else
		{
			if (tm.y < tm.z) { return bits.y; } // XZ plane
		}
		return bits.z; // XY plane;
	}*/
	
uint nextNode(vec3 tm)
{
	// Distance to exit cube. I believe the efficient sparse voxel octree' paper
	// tracks the cube size explicitly, but I don't think we need to?
	float tcMax = min(tm.x, min(tm.y, tm.z));
	uint result = 0;

	if (tm.x <= tcMax) { result |= 0x1; }
	if (tm.y <= tcMax) { result |= 0x2; }
	if (tm.z <= tcMax) { result |= 0x4; }

	return result;
}

uint getNode(uint node, uint childId)
{
	return dagData[node][childId];
}

// I hope the naming here is ok, and doesn't
// conflict with built-in GLSL mix() function?
uint mix(uint h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return h;
}

float footprintSize(int depth, vec3 childTm, Ray3d ray)
{
	// Compute the child size from the parent node depth.
	float childNodeSize = float(0x01 << (32 - (depth + 1)));

	// Alternatively we can track the child node size as we move up and down
	// the tree, optionally packed into a fourth component of the child step.
	// Another option is to compute is as follows (test on newer GPU):
	// float childNodeSize = min3((childT1 - childT0) * ray.mDir);

	vec3 temp = ray.mDir * childTm; // Element-wise
	float dist = length(temp);
	float projSize = childNodeSize / dist;
	return projSize;
}

vec3 childIdToT1Offset(uint childId, vec3 tChildStep)
{
	// Map e.g. (010) to (1.0, 0.0, 1.0). Approach below seems faster
	// than using a lookup table but maybe there is a better way?
	uint negChildId = ~childId;
	vec3 negChildIdAsVec3 = {
		(negChildId & 0x1),
		(negChildId & 0x2) >> 1,
		(negChildId & 0x4) >> 2
	};
	vec3 t1Offset = negChildIdAsVec3 * tChildStep;
		
	/*uvec3 mask = {(childId & 0x1), (childId & 0x2) >> 1, (childId & 0x4) >> 2};
	mask = mask - uvec3(1, 1, 1);
	vec3 t1Offset = uintBitsToFloat(floatBitsToUint(tChildStep) & mask);*/
	
	return t1Offset;
}

// This is the core function used to intersect a ray against our DAG. The
// algorithm was originally based on 'An Efficient Parametric Algorithm for
// Octree Traversal' by J. Revelles et al, but was later reworked to take a
// ideas from 'Efficient Sparse Voxel Octrees – Analysis, Extensions, and
// Implementation' by S. Laine et al. The concepts are very similar but the
// latter is more suitable for GPU implementation.
//
// It is an iterative implementation of what is naturally a recursive
// algorithm - we recursively check the ray against each child of the node
// until an intersection is found or the ray leaves the node. A true recursive
// implementation also exists (in C++ only) which is useful for debugging.
//
// This iterative version was developed against a GTX 660 which is already
// 10 years old at the time of writing and may have different performance
// characteristics to modern GPUs. This may be worth testing at some point.
// 
// When porting the recursive algorithm to an iterative implementation I
// found that the single biggest factor affecting performance was to avoid
// pushing and popping too much data from the stack. Therefore node index and
// child id get packed together, and no other parameters get stored.
//
// There are a couple more optimisation opportunities I could potentially
// explore in the future:
//   - The child index is a bitfield which gets updated based on comparisons
//     with a threshold, but using a bvec in combination with GLSL functions
//     like 'lessThanEqual()' might allow better parallelism here. But I'd
//     probably still need to convert the bvec to/from a bitfield to store it
//     on the stack. 
//   - Instead of carefully iterating over only those children which we know
//     intersect the ray, would could iterate over all children and do an
//     intersection test against each. It would simplify the stepping forward
//     logic at the expense of some extra tests. I think something like this might
//     be described in 'The HERO Algorithm for Ray Tracing Octrees'. Rather than
//     sorting children I think we already know how to traverse them near-to-far.
//
// My old GTX 660 and drivers do not really like these nested conditional
// statements in this function and prefer to see some combined into a flatter
// hierarchy (without the 'continue's) for about a 10% speed gain. But I'm
// leaving the hierarchy in place because it makes the code clearer, reveals
// logical optimisations (e.g. only read node value if it is in front of the
// camera), and I'm hoping newer hardware handles conditionals/continue's better.
//
// I've also seen evidence that returning early (in the middle of a function/loop)
// has another 10% performance penalty, but again I'm hoping this is a result
// of old GPU/drivers and am choosing to accept this for simpler code.
void intersectSubtreeIterative(vec3 t0, vec3 t1, uint nodeIndexIn, out RayVolumeIntersection intersection, int subDagDepth, Ray3d ray, float maxFootprint)
{
	int depth = subDagDepth;

	vec3 childStep = (t1 - t0) * 0.5f;
	vec3 tm = t0 + childStep;
	float nodeEntry = max3(t0);

	// We could pack these into a single bitfield but GLSL seems to prefer that
	// we don't. Therefore we keep them separate and only pack when placing on the
	// stack, at which point saving stack accesses/memory makes a big difference.
	uint nodeIndex = nodeIndexIn;
	uint childId = lessThanOrEqualTo(tm, nodeEntry); // First child
	
	uint nextChildFlips = 0;

	uint states[32];

	do
	{
		// We compute childT1 from t1 on each iteration. An alternative approach is to only
		// compute childT1 when moving up or down the tree, and then adjust it incrementally 
		// when advancing on the same level. I did not find a measurable performance difference
		// between the two approaches and recomputing it on each iteration is simpler.
		vec3 childT1 = t1;
		if((childId & 0x1) == 0) { childT1.x -= childStep.x; }
		if((childId & 0x2) == 0) { childT1.y -= childStep.y; }
		if((childId & 0x4) == 0) { childT1.z -= childStep.z; }

		// Test if node exit is in front of ray start and node is occupied.
		// Using min3() seems faster than all(greaterThan(childT1, vec3(0.0)))
		if (min3(childT1) > 0.0) 
		{
			// Only touch memory once we are in front of the ray start.
			uint childNodeIndex = getNode(nodeIndex, childId ^ a);
			
			if(childNodeIndex > 0) // Child is occupied
			{
				vec3 childT0 = childT1 - childStep;
				vec3 childTm = (childT0 + childT1) * 0.5f; // Midpoint
				
				float tEntry = max3(childT0);

				// If we have an internal (non-leaf) node we can traverse it further.
				if(childNodeIndex >= MaterialCount)
				{
					// Apply dithering to LOD transition. Does not actually look very nice
					// so might not keep it but I'll see how it looks with sub-pixel voxels.
					// Might it instead be useful for shadow rays rather than primary rays?
					// Setting ditherRange to 1.0 means some footprints get scaled by a factor
					// of two, and hence adjacent transitions start to touch each other.
					// Setting this to zero also seems to let the optimiser remove scaling code.
					float ditherRange = 0.3; // Typically 0.0 (disabled) - 1.0
					
					// Dither based on parent node (not including child id). This means all
					// children get the same scale factor, and hence preserve the size
					// relationship. If we don't traverse into a nearer child then we won't 
					// traverse into a more distant one either.
					float footprintScale = float(mix(nodeIndex) & 0xffff); // 0 - 65535
					footprintScale *= (1.0 / 65535.0f); // 0.0 - 1.0
					footprintScale *= ditherRange; // 0.0 - ditherRange
					footprintScale += 1.0f; // 1.0 - ditherRange
					float scaledMaxFootprint = footprintScale * maxFootprint;
					
					// Traverse further is the node is large in screen space.
					if (footprintSize(depth, childTm, ray) >= scaledMaxFootprint)
					{
						// Push the current state onto the stack. Packing the variables
						// together reduces the maximum node index we can use, but it
						// is still huge and a smaller stack really helps performance.
						assert((nodeIndex >> 29) == 0); // Check upper three bits are clear.
						states[depth] = (nodeIndex << 3) | childId;
						
						nodeIndex = childNodeIndex;
						childId = lessThanOrEqualTo(childTm, tEntry); // First node

						t1 = childT1;
						childStep *= 0.5f;
						depth++;
						continue;
					}
				}

				// Node is occupied and we chose not to descend, so this is where we stop.
				// But if we didn't yet get to a leaf node then we need to peek down the
				// tree to find the correct material.
				// FIXME - We should be able to skip determining the material (and the
				// normal?) for shadow rays, but doing so currently appears to trigger a
				// bug possibly related to code being incorrectly optimised out. So I'll
				// come back to this with new GPU/drivers.
				while (childNodeIndex >= MaterialCount)
				{
					const uint nearToFar[] = { 0x00, 0x01, 0x02, 0x04, 0x03, 0x05, 0x06, 0x07 };
					for (uint i = 0; i < 8; i++)
					{
						uint bt = nearToFar[i];
						uint childId = a ^ bt;
						uint nextLevelChildIndex = getNode(childNodeIndex,childId);
						// Zero is a material used to denote empty space. But for the purpose of this function
						// we don't want to include it as it doesn't help provide a valid material for rendering.
						if (nextLevelChildIndex > 0)
						{
							childNodeIndex = nextLevelChildIndex;
							break;
						}
					}
				}

				intersection.material = childNodeIndex;
				intersection.distance = tEntry;
				
				// FIXME - Normal calculation can probably also be skipped for shadow rays.
				if (childT0.x == tEntry)
				{
					intersection.normal = vec3(-1.0, 0.0, 0.0);
				}
				if (childT0.y == tEntry)
				{
					intersection.normal = vec3(0.0, -1.0, 0.0);
				}
				if (childT0.z == tEntry)
				{
					intersection.normal = vec3(0.0, 0.0, -1.0);
				}

				// Flip normals if required
				if (bool(a & 1)) { intersection.normal[0] *= -1.0; }
				if (bool(a & 2)) { intersection.normal[1] *= -1.0; }
				if (bool(a & 4)) { intersection.normal[2] *= -1.0; }
				
				intersection.packedNodeAndChild = (nodeIndex << 3) | childId;
				
				return;

			}
		}

		// Advance forward to the next child node
		nextChildFlips = lessThanOrEqualTo(childT1, min3(childT1)); // Next node
		childId ^= nextChildFlips;

		// If we could not advance forward then move up a level in the tree and try moving
		// forward from there. Keep doing this until we do successfully move forward. It
		// might seem natural to merge this while loop with the main one to eliminate some
		// duplication but this is not as easy as it seems. Also keeping this separate tighter
		// loop lets us avoid checking against T1 and the node index as we move up the tree
		// (we know parent nodes are valid, as we came from them previously).
		while ((childId & nextChildFlips) == 0)
		{
			// Terminate with a 'miss' if we reach the start depth (and so can't go higher).
			if(depth == subDagDepth) { return; }

			depth--;
			assert(depth >= 0);
			childStep *= 2.0f;
			
			// Pop off the previous state and unpack.
			nodeIndex = states[depth] >> 3;
			childId = states[depth] & 0x7;
			
			childT1 = t1;

			// We do not store t1 on a stack (this appeared to be slow) as we move
			// down the tree so we have to recompute it as we move back up the tree.
			// It is sort of the reverse of computing childT1 from t1 on each iteration
			// and feels a little redundant, but I haven't yet found a better way.
			if((childId & 0x1) == 0) { t1.x += childStep.x; }
			if((childId & 0x2) == 0) { t1.y += childStep.y; }
			if((childId & 0x4) == 0) { t1.z += childStep.z; }

			nextChildFlips = lessThanOrEqualTo(childT1, min3(childT1)); // Next node
			childId ^= nextChildFlips;
		}

	} while (true);
}

RayVolumeIntersection intersectNodes(Ray3d ray, float maxFootprint)
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
			if (childNodeIndex > 0)
			{
				if (childNodeIndex < MaterialCount)
				{
					intersection.material = childNodeIndex;
					intersection.distance = max(max(tx0, ty0), tz0);
					intersection.normal[0] = 0.0;
					intersection.normal[1] = 0.0;
					intersection.normal[2] = 0.0;
					if (tx0 > ty0 && tx0 > tz0)
					{
						intersection.normal[0] = -1.0;
					}
					if (ty0 > tx0 && ty0 > tz0)
					{
						intersection.normal[1] = -1.0;
					}
					if (tz0 > tx0 && tz0 > ty0)
					{
						intersection.normal[2] = -1.0;
					}

					// Flip normals if required
					if (bool(a & 1)) { intersection.normal[0] *= -1.0; }
					if (bool(a & 2)) { intersection.normal[1] *= -1.0; }
					if (bool(a & 4)) { intersection.normal[2] *= -1.0; }
				}
				else
				{
					vec3 t0 = { tx0, ty0, tz0 };
					vec3 t1 = { tx1, ty1, tz1 };
					intersectSubtreeIterative(t0, t1, childNodeIndex, intersection, int(subDAGs[childId].depth), reflectedRay, maxFootprint);
				}
			}
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

bool addNoise = true;
bool includeSun = true;
bool includeSky = true;

uint nextPointInUnitSphere = 17;
const uint pointsInUnitSphereCount = 97;
vec3 pointsInUnitSphere[pointsInUnitSphereCount] = {
{-0.383666, -0.804919, 0.0944412}, {-0.443004, -0.623236, 0.093763}, {0.596212, 0.600561, -0.405941},
{0.00732541, 0.311481, 0.595857}, {0.35747, -0.202523, 0.51548}, {0.481295, 0.486265, -0.0504826},
{-0.215546, -0.155825, 0.310956}, {-0.366899, -0.446154, 0.744858}, {0.389657, -0.749635, -0.365801},
{-0.931108, 0.327211, -0.122511}, {-0.748207, -0.236883, -0.579582}, {-0.0204712, -0.0840217, -0.108828},
{-0.0248622, 0.292626, 0.58795}, {0.615062, -0.44795, 0.411548}, {0.421408, -0.674777, 0.287922},
{-0.762005, -0.0879344, -0.00327188}, {-0.319229, 0.753515, 0.170535}, {0.502534, 0.642492, -0.48981},
{0.398153, -0.174667, 0.781806}, {-0.15367, 0.918583, 0.161913}, {0.094431, -0.683885, -0.722751},
{0.62857, -0.400337, -0.51295}, {0.232089, 0.557795, -0.0534223}, {0.661657, 0.587195, 0.170528},
{0.189007, 0.0994473, 0.465597}, {0.35964, 0.5144, -0.215359}, {0.507458, 0.123115, -0.239108},
{-0.583864, 0.135643, 0.0547429}, {-0.294475, 0.0615951, 0.185648}, {0.137647, -0.210184, -0.0612187},
{-0.325755, -0.22286, -0.675635}, {-0.37757, 0.725356, 0.0570663}, {0.203964, -0.0560864, -0.474057},
{-0.319561, 0.308158, 0.0596839}, {0.378429, 0.432201, 0.496303}, {0.0109971, 0.826675, 0.116538},
{-0.695244, 0.00638008, 0.651634}, {-0.0750516, 0.0766848, 0.0931839}, {0.708902, -0.114643, 0.208463},
{-0.273627, 0.737389, 0.359039}, {-0.831128, -0.307533, -0.200435}, {0.600137, 0.320239, -0.137172},
{-0.844886, 0.159409, 0.254769}, {0.360574, 0.706062, 0.0678662}, {0.24411, -0.122666, -0.298095},
{-0.600898, 0.026499, -0.723997}, {-0.196384, -0.235334, -0.848067}, {0.00954199, -0.520095, 0.64521},
{-0.504305, -0.0214947, 0.0881122}, {0.754728, -0.220522, -0.579396}, {-0.516617, -0.454494, -0.192176},
{-0.015116, -0.807091, -0.55316}, {-0.129509, -0.53044, 0.105762}, {0.618274, 0.298231, -0.483408},
{0.463445, -0.740307, 0.295492}, {-0.0133463, -0.0981526, -0.242999}, {0.0940177, 0.43694, -0.407358},
{-0.42439, 0.489386, 0.246871}, {-0.62209, 0.623748, 0.373551}, {-0.374984, -0.632978, -0.228579},
{-0.263031, -0.312141, 0.251237}, {0.631538, 0.560455, 0.33857}, {-0.0830062, 0.551425, -0.392772},
{-0.0264167, 0.82113, -0.128283}, {-0.0227646, -0.106432, 0.51158}, {-0.387301, -0.783876, 0.0170174},
{-0.21282, 0.0215431, 0.740373}, {0.635255, -0.278279, 0.589663}, {0.834236, 0.288636, 0.087611},
{-0.242781, -0.719712, 0.623161}, {0.100313, -0.22277, 0.24495}, {0.559379, 0.174089, -0.0467238},
{-0.584515, 0.080276, -0.397507}, {0.771996, -0.610471, -0.117553}, {-0.520995, -0.544671, 0.599306},
{-0.128603, -0.0539699, -0.377795}, {-0.139585, 0.266127, -0.630367}, {0.168764, 0.809762, 0.453309},
{0.360813, -0.777762, 0.414643}, {-0.483871, -0.675343, -0.18256}, {-0.732527, 0.189792, -0.100888},
{0.594728, 0.422432, -0.664888}, {-0.350073, -0.406648, 0.31156}, {-0.362443, -0.15962, -0.151666},
{0.563818, 0.0157166, -0.773615}, {-0.022782, -0.311344, 0.15705}, {0.211856, 0.0936115, -0.546897},
{0.0422716, -0.620189, -0.536811}, {-0.196857, -0.0222045, -0.173419}, {0.24812, -0.659695, 0.358271},
{-0.401896, -0.20897, 0.272937}, {-0.344404, -0.329286, -0.741868}, {0.359456, -0.416078, -0.726894},
{0.694219, 0.442455, -0.538471}, {-0.786476, 0.00394881, 0.307515}, {0.558103, -0.426827, 0.430074},
{0.781845, -0.289303, -0.331674}
};

vec3 surfaceColour(RayVolumeIntersection intersection)
{
	vec3 colour = getMaterial(intersection.material, materials).xyz;

	if (addNoise)
	{
		// Noise is applied multiplicatively as this avoids creating overshoots
		// and undershoots for saturated or dark surface colours respectively.
		float noise = positionBasedNoise(vec4(intersection.position, 1.0));
		noise = (noise * 0.1) + 0.9; // Map 0.0 - 1.0 to range 0.9 - 1.0.
		colour *= noise;
	}

	return colour;
}

vec3 randomPointInUnitSphere()
{
	nextPointInUnitSphere = nextPointInUnitSphere++;
	vec3 result = pointsInUnitSphere[nextPointInUnitSphere % pointsInUnitSphereCount];

	return result;
}

vec3 gatherLighting(vec3 position, vec3 normal)
{
	const vec3 offset = normal * 0.001f;
	float intensity = 0.0f;

	if (includeSun)
	{
		float sunIntensity = 2.0f;
		vec3 sunDir = normalize(vec3(1.0, -2.0, 10.0 ));

		Ray3d sunShadowRay;
		sunShadowRay.mOrigin = position + offset;
		sunShadowRay.mDir = sunDir;
		RayVolumeIntersection sunShadowIntersection = intersectNodes(sunShadowRay, gMaxFootprint);
		if (sunShadowIntersection.material == 0)
		{
			intensity += max(dot(sunDir, normal), 0.0f) * sunIntensity;
		}
	}

	/*if (includeSky)
	{
		float skyIntensity = 1.0f;
		const vec3 skyDir = normalize(normal + randomPointInUnitSphere());
		Ray3d skyShadowRay;
		skyShadowRay.mOrigin = position + offset;
		skyShadowRay.mDir = skyDir;
		RayVolumeIntersection skyShadowIntersection = intersectNodes(skyShadowRay, gMaxFootprint);
		if (skyShadowIntersection.material == 0)
		{
			intensity += skyIntensity;
		}
	}*/

	return vec3(intensity);
}

vec4 traceSingleRay(const Ray3d ray, uint depth)
{
	vec4 pixelColour = { 0.0f, 0.0f, 0.0f, 0.0f };

	RayVolumeIntersection intersection0 = intersectNodes(ray, gMaxFootprint);
	if (intersection0.material > 0)
	{
		vec3 surfCol0 = surfaceColour(intersection0);
		vec3 directLighting0 = gatherLighting(intersection0.position, intersection0.normal);

		/*const vec3 reflectedDir = normalize(intersection0.normal + randomPointInUnitSphere());
		Ray3d reflectedRay;
		reflectedRay.mOrigin = intersection0.position + (intersection0.normal * 0.01);
		reflectedRay.mDir = reflectedDir;*/

		vec3 indirectLighting0 = {0.5f, 0.5f, 0.5f};
		/*RayVolumeIntersection intersection1 = intersectNodes(reflectedRay, gMaxFootprint);
		if (intersection1.material > 0)
		{
			vec3 surfCol1 = surfaceColour(intersection1);
			vec3 directLighting1 = gatherLighting(intersection1.position, intersection1.normal);
			indirectLighting0 = surfCol1 * directLighting1;
		}*/

		pixelColour.rgb = surfCol0 * (directLighting0 + indirectLighting0);

		pixelColour.a = float(intersection0.packedNodeAndChild & 0xffff);

	}

	return pixelColour;
}

// See http://byteblacksmith.com/improvements-to-the-canonical-one-liner-glsl-rand-for-opengl-es-2-0/
float rand(vec2 co)
{
    float a = 12.9898;
    float b = 78.233;
    float c = 43758.5453;
    float dt= dot(co.xy ,vec2(a,b));
    float sn= mod(dt,3.14);
    return fract(sin(sn) * c);
}

uint hashRay(Ray3d ray)
{
	uint result = 0;
	result ^= mix(floatBitsToInt(ray.mOrigin.x));
	result ^= mix(floatBitsToInt(ray.mOrigin.y));
	result ^= mix(floatBitsToInt(ray.mOrigin.z));
	result ^= mix(floatBitsToInt(ray.mDir.x));
	result ^= mix(floatBitsToInt(ray.mDir.y));
	result ^= mix(floatBitsToInt(ray.mDir.z));
	return result;
}

void main()
{
	a = 0;
	Ray3d ray;
	
	bool multisample = false;
	if(multisample)
	{
		vec2 windowSize = vec2(1600, 1200); // FIXME - Hardcoded!
		vec2 offset = vec2(1.0 / windowSize);
		vec2 lowerLeft = TexCoords - (offset * 0.5);
		
		ray.mOrigin = cameraPos;
		ray.mDir = createRay(lowerLeft + vec2(0, 0), PInv, VInv);
		nextPointInUnitSphere = hashRay(ray);
		FragColor = traceSingleRay(ray, 0);
		
		ray.mOrigin = cameraPos;
		ray.mDir = createRay(lowerLeft + vec2(offset.x, 0), PInv, VInv);
		nextPointInUnitSphere = hashRay(ray);
		FragColor += traceSingleRay(ray, 0);
		
		ray.mOrigin = cameraPos;
		ray.mDir = createRay(lowerLeft + vec2(0, offset.y), PInv, VInv);
		nextPointInUnitSphere = hashRay(ray);
		FragColor += traceSingleRay(ray, 0);
		
		ray.mOrigin = cameraPos;
		ray.mDir = createRay(lowerLeft + vec2(offset.x, offset.y), PInv, VInv);
		nextPointInUnitSphere = hashRay(ray);
		FragColor += traceSingleRay(ray, 0);
		
		FragColor /= 4.0;
	}
	else
	{
		ray.mOrigin = cameraPos;
		ray.mDir = createRay(TexCoords, PInv, VInv);
		// This variable controls the start point in the sequence of pseudorandom unit vectors.
		// It needs to be pseudorandom (at least locally) because if it is constant then all rays
		// which hit a given surface will reflect in the same way. Basing it on fragment position
		// is not sufficient, because any visible noise pattern is then static in screenspace
		// and which looks strange as the camera moves and rotates. Hence we base it on the ray
		// so that it is static when the camera is still but changes randomly as the camera is
		// transformed. It is also possible to mix in frame number which causes successive
		// frames to be different. This is necessary for progressive rendering but probably
		// undesirable otherwise (and it requires external code to provide the frame number).
		nextPointInUnitSphere = hashRay(ray);
		FragColor = traceSingleRay(ray, 0);
	}


	// For testing noise removal
	//float noise = float(hashRay(ray) & 0xFFFF) / 65535.0;
	//noise *= 2.0f; // Range 0.0 - 2.0 (so centred at 1.0)
	//FragColor.rgb *= noise;
	
	/*uvec2 windowPos = uvec2(gl_FragCoord.xy);
	
	if((windowPos.x % 100 == 0) || (windowPos.y % 100 == 0))
	{
		FragColor = vec4(1.0, 1.0, 1.0, 1.0);
	}
	else
	{
	FragColor = vec4(0.0, 0.0, 0.0, 1.0);
	}*/
} 