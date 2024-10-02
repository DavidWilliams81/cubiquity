#version 460 core

// The full DAG has 33 levels, from zero (for leaves) to 32 (for the root).
const int RootNodeHeight = 32;

const float FLT_MAX = 3.402823466e+38; // See https://stackoverflow.com/q/16069959

////////////////////////////////////////////////////////////////////////////////
// Utility code shared between instancing and pathtracing shaders
////////////////////////////////////////////////////////////////////////////////

// Definitions to aid writing code that compiles as both GLSL and C/C++
void assert(bool unused) {}
const int INT_MIN   = (-2147483647 - 1);
const int INT_MAX   =  2147483647;
const uint UINT_MAX =  0xffffffff;
#define COMPILE_AS_GLSL 1 // Not defined in C++ version

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

// This structure contains an explicit 'hit' member rather than relying on a sentinal value
// such as a negative distance or a material of zero. I have found this simplifies the code
// and we might want to skip material identification anyway (if traversal is stopped early
// due to LOD). There is some redundancy storing both the distance and the position, but we
// get the former basically for free as we traverse and the latter is very often useful.
struct RayVolumeIntersection
{
	bool hit;
	float distance;
	uint material;
	vec3 position;
	vec3 normal;
};

struct SubDAG
{
	ivec3 lowerBound;
	int nodeHeight;

	uint padding0;
	uint nodeIndex;
	uint padding1, padding2;
};

// FIXME - Rather than hard-coding these bindings I think we should be able to retrieve them through
// glGetProgramResourceIndex() and/or glShaderStorageBlockBinding(). Unfortunately that is currently
// crashing for me inside the driver. Might be worth testing again after an upgrade/update. See
// https://www.khronos.org/opengl/wiki/Interface_Block_(GLSL)#OpenGL_binding_index_setting
// https://web.archive.org/web/20150109112923/http://www.lighthouse3d.com/tutorials/glsl-core-tutorial/3490-2/
layout(binding = 0, std430) buffer SomeName  
{
	// Note: It should be possible to declare the dagData as a real 2D array as below,
	// to avoid manually doing the index calculations. However, there seems to be a
	// compiler bug on my GTX 660 which stops this work properly. It can be fixed by
	// specifying the size of the first dimension, but that is undesirable and also
	// triggers a second bug in which the shader takes a long time to compile (many
	// minutes) as the data gets large (many megabytes).
	// Therefore we stick to a single dimensional array and handle 2D lookups manually.
	//uint dagData[][8];
	//uint dagData[2000000][8];
	uint dagData[];
};

layout(binding = 1, std430) buffer SomeOtherName  
{
	SubDAG subDAGs[8];  
};

struct Ray3f
{
	vec3 mOrigin;
	vec3 mDir;
};

const uint MaterialCount = 256;
const float gMaxFootprint = 0.0035;

uint getNode(uint node, uint childId)
{
	return dagData[node * 8 + childId];
}

bool isMaterialNode(uint nodeIndex)
{
	return nodeIndex < MaterialCount;
}

ivec3 childIdToIVec3(int childId)
{
	return ivec3((childId) & 0x1, (childId >> 1) & 0x1, (childId >> 2) & 0x1);
}

SubDAG findSubDAG(uint rootNodeIndex, uint childId)
{
	// Initialised for root, but updated on first iteration of the loop.
	int childHeight = 32;
	ivec3 lowerBound = { INT_MIN , INT_MIN , INT_MIN };

	// We never return the root node as a subDAG, instead
	// we always descend into at least the first child.
	uint onlyChildId = childId;
	uint nextNodeIndex = getNode(rootNodeIndex, onlyChildId);
	uint nodeIndex = 0;
	uint childCount = 1;
	
	// Executed at least once as initialised to '1' above
	while (childCount == 1)
	{
		childHeight--;
		nodeIndex = nextNodeIndex;

		ivec3 childIdAsIVec3 = childIdToIVec3(int(onlyChildId));
		childIdAsIVec3 <<= childHeight;
		lowerBound ^= ivec3(childIdAsIVec3);

		childCount = 0;
		for (uint i = 0; i < 8; i++)
		{
			uint childNodeIndex = getNode(nodeIndex, i);
			if (childNodeIndex > 0)
			{
				// Store value we just received to avoid GPU memory access
				// cost of getting it again on next iteration of outer loop.
				nextNodeIndex = childNodeIndex;
				childCount++;
				onlyChildId = i;
			}
		}
	}

	SubDAG subDAG;
	subDAG.nodeIndex = nodeIndex;
	subDAG.lowerBound = lowerBound;
	subDAG.nodeHeight = childHeight;

	return subDAG;
}

SubDAG getSubDAG(uint rootNodeIndex, uint childId)
{
	int method = 1;
	switch (method) {
	case 1:
		// Use the real precomputed SubDAG
		return subDAGs[childId];
	case 2:
		// Compute the SubDAG on demand
		return findSubDAG(rootNodeIndex, childId);
	case 3:
	default:
		// Bypass use of the SubDAG by returning one which directly references the relevant child of the
		// root. This means the ESVO algorithm has to traverse the full tree. It's a slower option, but is
		// the simplest and can be used to validate the pecision of the ESVO traversal for large volumes.
		SubDAG bypassSubDAG;
		bypassSubDAG.lowerBound[0] = bool(childId & 0x1) ? 0 : INT_MIN;
		bypassSubDAG.lowerBound[1] = bool(childId & 0x2) ? 0 : INT_MIN;
		bypassSubDAG.lowerBound[2] = bool(childId & 0x4) ? 0 : INT_MIN;
		bypassSubDAG.nodeHeight = RootNodeHeight - 1;
		bypassSubDAG.nodeIndex = getNode(rootNodeIndex,childId);
		return bypassSubDAG;
	}
}

// Find the material of the nearest occupied child based on the direction ray is travelling (i.e.
// the direction we are viewing the node from). Note that this does not use an ESVO-type traversal
// (though perhaps it should, I haven't tested) but instead simply iterates over all child nodes
// in near-to-far order until an occupied one is found. The traversal approach is described here:
//
// https://www.flipcode.com/archives/Harmless_Algorithms-Issue_02_Scene_Traversal_Algorithms.shtml#octh
//
// It does not check for intersection and so might be faster than an ESVO-type approach (less logic,
// but maybe more memory accesses?), at the expense of some precision (the ray might actually miss
// the nearest occupied child). I'm assuming this doesn't matter too much as voxels are tiny on screen.
uint findNearestMaterial(uint nodeIndex, uint rayDirSignBits)
{
	// We use a clever packing trick for the array of child IDs to iterate over. Each of the 8 child
	// IDs (0-7) needs only three bits, and a fourth 'valid' bit is used to indicate that a value is
	// stored. We iterate over the values with right-shifting, so the valid bit is cleared for new
	// values which lets us know when we are done.
	// I'm not sure it's really any better than a simple array, but is more compact and does let us do
	// the reflection (based on ray dir sign) in one XOR prior to the loop, rather than per-iteration.
	uint nearToFarPacked = 0x76534210; // '3' and '4' swapped on purpose, see link above.
	nearToFarPacked |= 0x88888888; // Set every 4th bit high as a marker

	const uint rayDirSignBitsDup = rayDirSignBits * 0x11111111; // Duplicate lower nibble across uint
	const uint orderedChildIds = nearToFarPacked ^ rayDirSignBitsDup; // Reflect all ids in one go

	while (nodeIndex >= MaterialCount)
	{
		for (uint childIds = orderedChildIds; childIds != 0; childIds >>= 4)
		{
			uint childNodeIndex = getNode(nodeIndex, childIds & 0x7);
			if (childNodeIndex > 0) // Skip empty nodes
			{
				nodeIndex = childNodeIndex;
				break;
			}
		}
	}

	return nodeIndex;
}

uint bitMix(uint h)
{
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;

	return h;
}

ivec3 findFirstChild(float nodeEntry, vec3 rayOrigin, vec3 invRayDir, ivec3 nodeCentre)
{
	// Distances to midpoint along each axis.
	vec3 nodeTm = (vec3(nodeCentre) - rayOrigin) * invRayDir;

	// Find first child by determining the axes along which the ray has already passed the midpoint.
	// Logic reversed compared to ESVO paper, probabaly because they have a negative ray direction?
	ivec3 childId = ivec3(lessThan(nodeTm, vec3( nodeEntry, nodeEntry, nodeEntry )));

	// As an extension to the ESVO paper we skip children which are behind the ray origin. In many
	// cases there aren't any (as the ray origin is usually outside the node), but there can be
	// during the initial part of the traversal if the ray starts inside the volume. It saves us
	// having to check that the node intersection distance is greater than zero on every iteration.
	if (nodeEntry <= 0) { // Only do skipping of children if node is at least partly behind camera.
		childId |= ivec3(greaterThanEqual(rayOrigin, vec3(nodeCentre)));
	}

	return childId;
}

// This is the core function used to intersect a ray against our DAG. The algorithm is based on
// 'Efficient Sparse Voxel Octrees – Analysis, Extensions, and Implementation' by S. Laine et al
// (hereafter referred to as the ESVO paper or algorithm).
//
// This algorithm was developed on a GTX 660 which is already 10 years old at the time of writing
// and may have different performance characteristics to modern GPUs. I found that one of the
// biggest factors affecting performance was pushing and popping too much data on the stack. 
// possible cause of this is that the stack may get allocated in local memory and accessing
// (perhaps especially writing? this is much slower than accessing registers. Fortunately the ESVO
// algorithm is quite conservative in this regard. See https://stackoverflow.com/q/49818414 
//
// Note: Compared to the ESVO paper we do not need to maintain the 'active span' because that is
// used to culling child voxels against the contours (which are potentially stored for eah level
// and combined when descending the tree). ESVO paper also tracks child entry point incrementally,
// but we don't need to do that here (the exit point is enough).
RayVolumeIntersection intersectRayNodeESVO(
	uint nodeIndex, ivec3 nodePos, int nodeHeight,
	Ray3f ray, vec3 rayDirSign, uint rayDirSignBits,
	bool computeSurfaceProperties, float maxFootprint)
{
	RayVolumeIntersection intersection = { false, 0, 0, {0, 0, 0}, {0, 0, 0} }; // Miss

	uint nodeSize = 1u << nodeHeight;

	// Intersect ray wth the node. The cast to float before
	// addition prevents integer overflow for largest nodes.
	vec3 invRayDir = vec3( 1,1,1 ) / ray.mDir;
	vec3 nodeT0 = (vec3(nodePos) - ray.mOrigin) * invRayDir;
	vec3 nodeT1 = ((vec3(nodePos) + float(nodeSize)) - ray.mOrigin) * invRayDir;
	
	float nodeEntry = max3(nodeT0);
	float nodeExit = min3(nodeT1);

	// Check if the ray actually hits the node.
	if (nodeEntry < nodeExit)
	{
		const int startHeight = nodeHeight;

		// Child size fits in signed type as it is half node size.
		int childNodeSize = int(nodeSize / 2);

		// We store the child id as an vector type whereas the ESVO paper uses a simple uint bitfield.
		// Using the vector type gives shorter, cleaner code but does appear to be slightly slower.
		// Later on we have to unpack our vector type to a uint bitfield when we use it as an index.
		ivec3 childId = findFirstChild(nodeEntry, ray.mOrigin, invRayDir, nodePos + childNodeSize);

		// Child position is stored as integers to maximise precision for when dealing with very large
		// volumes. But we do quite often cast to floats, and floats can hold a wide range of integer
		// values anyway. We may want to reassess this (maybe cache the float values?) when we actually
		// have large volumes for testing.
		ivec3 childPos = nodePos + childId * childNodeSize;

		float lastExit = nodeExit; // Used to prevent unnecessary stack writes ('h' in ESVO paper)
		uint nodeStack[RootNodeHeight + 1]; // Valid range is [0, RootNodeHeight] inclusive

		do
		{
			// The exit point is computed explicitly for every child. We could instead compute it for
			// only the first child and then increment it as we move through the children, but the tiny
			// floating point error which accumulates is enough to break the stack protection (because
			// the exit point of the last child may not exactly match the exit point of the parent node).
			vec3 childT1 = (vec3(childPos + childNodeSize) - ray.mOrigin) * invRayDir;
			float tChildExit = min3(childT1);
			assert(tChildExit > 0.0); // We only process node in front of the camera

			// Only touch memory once we are in front of the ray start.
			uint childIdBits = (childId[0] & 0x1) | ((childId[1] & 0x1) << 1) | ((childId[2] & 0x1) << 2);
			uint childNodeIndex = getNode(nodeIndex,childIdBits ^ rayDirSignBits);

			if (childNodeIndex > 0) // Child is occupied
			{
				// Compute the entry point
				vec3 childT0 = (vec3(childPos) - ray.mOrigin) * invRayDir;
				float tChildEntry = max3(childT0);

				// If we have an internal (non-leaf) node we can traverse it further.
				// Traverse further if the node is large in screen space.
				const bool isInternalNode = childNodeIndex >= MaterialCount;
				const bool hasLargeFootprint = (childNodeSize / tChildExit) > maxFootprint;

				// Descend if we can and if we want to.
				if (isInternalNode && hasLargeFootprint)
				{
					// Push the current state onto the stack if not already there.
					// Note: It's not clear that the stack write protection actually helps
					// performance. This should be tested on a more recent GPU. Without it
					// we might be able to update the child exit point incrementally.
					if (tChildExit < lastExit)
					{
						nodeStack[nodeHeight] = nodeIndex;
					}
					lastExit = tChildExit;

					// As we descend the current child becomes the new current node.
					nodeHeight--;
					nodeIndex = childNodeIndex;
					childNodeSize /= 2;

					// Detemine which child of the child we will enter first.
					childId = findFirstChild(tChildEntry, ray.mOrigin, invRayDir, childPos + childNodeSize);

					// Compute position of new child and restart loop to process it.
					childPos = childPos + childId * childNodeSize;
				}
				else // Node is occupied or we chose not to descend further, so this is where we stop.
				{
					intersection.hit = true;
					intersection.distance = tChildEntry;

					if (computeSurfaceProperties)
					{
						// Find the material of the intesected node (which might not be a leaf due to LOD).
						intersection.material = findNearestMaterial(childNodeIndex, rayDirSignBits);

						// Determining a true voxel normal is tricky because we don't have easy access to neighbouring
						// voxels. Therefore we currently determine a face normal instead by looking at which side of
						// the cube the ray entered. This is fine (even desirable, depending on user requirements) and
						// we can probably smooth out voxels in screen space if we want to. But I wonder if a true 
						// voxel normal would be better when bouncing rays around inside the volume for pathtracing.
						intersection.normal = vec3(equal(vec3( tChildEntry, tChildEntry, tChildEntry ), childT0));
						intersection.normal *= -rayDirSign;
					}
				}
			}
			else
			{
				// Advance forward to the next the child node (sibling). The ESVO paper uses an equality
				// test in the text but '<=' in the sample code. I don't know why (though I can see it
				// is valid) but perhaps it adds robustness in the face of floating point inaccuracy?
				ivec3 nextChildFlipsVec = ivec3(lessThanEqual(childT1, vec3( tChildExit, tChildExit, tChildExit )));
				childId ^= nextChildFlipsVec;

				ivec3 oldChildPos = childPos;
				childPos += nextChildFlipsVec * childNodeSize;

				// Check that we sucessfully moved forward in every direction in which we tried to
				// do so. If not we have left the node and it is time to pop. Note that the sample
				// code in the ESVO paper actually checks against zero  here - I *think* that is a
				// bug and they fail to corectly handle the case where the ray leaves a child via
				// the corner and wraps around on one axis while moving forward on another.
				if (!((childId & nextChildFlipsVec) == nextChildFlipsVec))
				{
					// Find the highest node boundary crossed by the last step. This may not be the
					// direct parent (as we can ascend multiple levels in one go). Use the unsigned
					// version of GLSL findMSB() which has different behaviour to signed version.
					ivec3 differingBits = oldChildPos ^ childPos;
					int combinedDiffBits = differingBits[0] | differingBits[1] | differingBits[2];
					int msb = findMSB(uint(combinedDiffBits)); // Unsigned version

					// The next node up is the one contining both our previous and next siblings.
					nodeHeight = msb + 1;

					// We will terminate traversal if we get higher than the start of our subDAG,
					// but we should never get higher than the root of the full DAG. This also means
					// we are still within the bounds of the stack.
					assert(nodeHeight <= RootNodeHeight);

					// Retrieve the node index and compute child node properties
					nodeIndex = nodeStack[nodeHeight];
					childNodeSize = 1 << msb;
					childId = (childPos >> msb) & 1;
					childPos = (childPos >> nodeHeight) << nodeHeight;
					childPos += childId * childNodeSize;

					lastExit = 0.0f; // Prevent parent being written to stack again
				}
			}
		} while (intersection.hit == false && nodeHeight <= startHeight);
	}
	
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
uniform uint frameId;

// Static uniform branching appears to be slow on my GTX 660 so I have
// chosen to compile separate GLSL programs for preview vs. progressive mode.
// https://stackoverflow.com/questions/37827216/do-conditional-statements-slow-down-shaders
#ifdef PROGRESSIVE
const bool includeSky = true;
const uint indirectBounces = 1;
#else
const bool includeSky = false;
const uint indirectBounces = 0;
#endif

out vec4 FragColor;

in vec2 TexCoords;

bool addNoise = true;
bool includeSun = true;

uint randSeed = 0;

uint nextPointInUnitSphere = 17;

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
	vec3 result;
	do
	{
		nextPointInUnitSphere = bitMix(nextPointInUnitSphere);
		result[0] = nextPointInUnitSphere & 0x3FF;
		result[1] = (nextPointInUnitSphere >> 10) & 0x3FF;
		result[2] = (nextPointInUnitSphere >> 20) & 0x3FF;

		vec3 offset = { 511.5f, 511.5f, 511.5f };
		result = result - offset;
		result /= 511.5f;

	} while (dot(result, result) >= 1.0f);
	
	return result;
}

// Intersect a ray with a volume by performing using the ESVO on each octant.
//
// A Cubiquity volume can be very large (2^32 voxels along each side) but it practice most volume
// are much smaller than this. It is most natural for the actual voxel geometry to be in the
// centre of the volume, which is also the worst case because it means all eight octants are
// occupied. This mens a typical DAG will have a root with eight children, but each of these
// children will then contain a long sequence of nodes with no splits before the real geometry
// starts. Rather than starting traversal at the root eight times we instead precompute the first
// split point, which becaomes a sort of virtual root for each of the eight 'SubDAGs'. For each
// of these subDAGs we jump straight to the virtual root and perform an ESVO traversal from there.
//
// We could potentially compute the bounds of the actual occupied part of the volume and of the
// subdags, in order to optimise for the case when the ray mises the volume/subdag completely.
// But for now I prefer to assume that the ray hits the volume (e.g. because the camera is inside
// the volume, or because it is a secondary ray) in order to optimise for this worst case. This
// also keeps things simpler until we have determined our real requirements. We might also consider
// rendering the occupied bounds (rather than a fullscreen quad) to generate a more conservative
// set of fragments.
//
// Note: This function does not implement special handling of the case where a component of the
// ray direction is zero. This is discussed in Section 3.3 of An 'Efficient Parametric Algorithm
// for Octree Traversal'. I think that the standard behaviour of IEEE 754 handling of +/-infinity
// and NaNs might be enough but I am not certain. If it proves to be a problem (if we ever see
// NaNs?) then it can be solved by nudging tiny direction components away from zero.
RayVolumeIntersection intersectVolume(Ray3f ray, bool computeSurfaceProperties, float maxFootprint)
{
	RayVolumeIntersection intersection = { false, 0, 0, {0, 0, 0}, {0, 0, 0} }; // Miss
	
	//const Internals::NodeStore& nodes = Internals::getNodes(volume).nodes();

	const uint rootNodeIndex = 256;

	// Store the sign of the ray direction
	ivec3 rayDirSignBitsAsVec = ivec3(lessThan(ray.mDir, vec3(0,0,0))); // 0 means +ve, 1 means -ve
	uint rayDirSignBits = rayDirSignBitsAsVec[0] | (rayDirSignBitsAsVec[1] << 1) | (rayDirSignBitsAsVec[2] << 2);
	vec3 rayDirSign = vec3(rayDirSignBitsAsVec * (-2) + 1); // Binary (-1, 1) not ternary (-1, 0, 1)

	// We reflect the ray so that it points along the positive axes. This simplifies the traversal.
	// Shifting the ray origin lets us pretend voxel corners have integer coordinates (not +/-0.5).
	Ray3f reflectedRay = ray;
	reflectedRay.mOrigin += vec3( 0.5f, 0.5f, 0.5f );
	reflectedRay.mOrigin = reflectedRay.mOrigin * rayDirSign;
	reflectedRay.mDir = abs(reflectedRay.mDir);

	// Work out our starting octant by determining which of the axes we are already in front of.
	// If we are in front of a given axis, the set the corresponding distance to a large value
	// so we ignore it when later searching for the minimum.
	int childId = 0;
	vec3 distToOrigin = (-reflectedRay.mOrigin) / reflectedRay.mDir;
	if (distToOrigin[0] < 0.0) { childId += 1; distToOrigin[0] += FLT_MAX; }
	if (distToOrigin[1] < 0.0) { childId += 2; distToOrigin[1] += FLT_MAX; }
	if (distToOrigin[2] < 0.0) { childId += 4; distToOrigin[2] += FLT_MAX; }

	// Iteration over the octants is similar to the ESVO approach, except that we consider
	// the octants to be (semi) infinite. Because the root is so large it is challenging
	// (and probably unnecessary) to do proper intersection calculations against the children.
	// This also means a ray will always reach the final octant (it will never miss it)
	// unless it hits geometry first, which simplifies the loop termination condition.
	do
	{
		SubDAG subDAG = getSubDAG(rootNodeIndex, childId ^ rayDirSignBits);
		uint childNodeIndex = subDAG.nodeIndex;

		if (childNodeIndex > 0)
		{
			// FIXME - Handle fully occupied direct child of root node.
			assert(childNodeIndex >= MaterialCount);

			// We can store the lower bound as integers because we offset ray origin by 0.5.
			// On reflected axis the 'upper' corner becomes lower than the 'lower' corner.
			// Node size might be large and negative, but I *think* the logic still works.
			int nodeSize = int(1u << subDAG.nodeHeight);
			ivec3 refNodeLowerBound = ivec3(subDAG.lowerBound) * ivec3(rayDirSign);
			refNodeLowerBound -= rayDirSignBitsAsVec * ivec3( nodeSize, nodeSize, nodeSize );

			intersection = intersectRayNodeESVO(
				childNodeIndex, refNodeLowerBound, subDAG.nodeHeight,
				reflectedRay, rayDirSign, rayDirSignBits, computeSurfaceProperties, maxFootprint);

			// Return on first hit. Position computed from real ray (not reflected version).
			if (intersection.hit) {
				intersection.position = ray.mOrigin + (ray.mDir * intersection.distance);
				break;
			}
		}

		// Determine which axis is closest to work out which child we intersect next.
		float minDist = min3(distToOrigin);
		if (distToOrigin[0] <= minDist) { childId += 1; distToOrigin[0] += FLT_MAX; }
		if (distToOrigin[1] <= minDist) { childId += 2; distToOrigin[1] += FLT_MAX; }
		if (distToOrigin[2] <= minDist) { childId += 4; distToOrigin[2] += FLT_MAX; }

	} while (childId <= 7); // 8 children, number 7 is the last.

	return intersection;
}

vec3 gatherLighting(vec3 position, vec3 normal)
{
	const vec3 offset = normal * 0.001f;
	vec3 intensity = vec3(0.0f);

	if (includeSun)
	{
		const vec3 sunColour = vec3(0.1, 0.1, 0.1);
		vec3 sunDir = normalize(vec3(1.0, -2.0, 10.0 ));

		Ray3f sunShadowRay;
		sunShadowRay.mOrigin = position + offset;
		sunShadowRay.mDir = sunDir;
		RayVolumeIntersection sunShadowIntersection = intersectVolume(sunShadowRay, false, gMaxFootprint);
		if (sunShadowIntersection.hit == false)
		{
			intensity += sunColour * max(dot(sunDir, normal), 0.0f);
		}
	}

	if (includeSky)
	{
		const vec3 skyColour = vec3(1.5f, 1.5f, 1.5f);
		const vec3 skyDir = normalize(normal + randomPointInUnitSphere());
		Ray3f skyShadowRay;
		skyShadowRay.mOrigin = position + offset;
		skyShadowRay.mDir = skyDir;
		RayVolumeIntersection skyShadowIntersection = intersectVolume(skyShadowRay, false, gMaxFootprint);
		if (skyShadowIntersection.hit == false)
		{
			intensity += skyColour;
		}
	}

	return intensity;
}

vec4 traceSingleRay(const Ray3f ray, uint depth)
{
	vec4 pixelColour = { 0.8f, 0.8f, 1.0f, 1.0f }; // Light blue background

	RayVolumeIntersection intersection0 = intersectVolume(ray, true, gMaxFootprint);
	if (intersection0.material > 0)
	{
		vec3 surfCol0 = surfaceColour(intersection0);
		vec3 directLighting0 = gatherLighting(intersection0.position, intersection0.normal);

		vec3 indirectLighting0 = {0.0f, 0.0f, 0.0f};
		if(indirectBounces > 0)
		{
			const vec3 reflectedDir = normalize(intersection0.normal + randomPointInUnitSphere());
			Ray3f reflectedRay;
			reflectedRay.mOrigin = intersection0.position + (intersection0.normal * 0.01);
			reflectedRay.mDir = reflectedDir;

			RayVolumeIntersection intersection1 = intersectVolume(reflectedRay, true, gMaxFootprint);
			if (intersection1.material > 0)
			{
				vec3 surfCol1 = surfaceColour(intersection1);
				vec3 directLighting1 = gatherLighting(intersection1.position, intersection1.normal);
				indirectLighting0 = surfCol1 * directLighting1;
			}
		}

		pixelColour.rgb = surfCol0 * (directLighting0 + indirectLighting0);

		pixelColour.a = 1.0f;

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

uint hashRay(Ray3f ray)
{
	uint result = 0;
	result ^= bitMix(floatBitsToInt(ray.mOrigin.x));
	result ^= bitMix(floatBitsToInt(ray.mOrigin.y));
	result ^= bitMix(floatBitsToInt(ray.mOrigin.z));
	result ^= bitMix(floatBitsToInt(ray.mDir.x));
	result ^= bitMix(floatBitsToInt(ray.mDir.y));
	result ^= bitMix(floatBitsToInt(ray.mDir.z));
	return result;
}

void main()
{

uint tileGroup = 0;
uint currentGroup = 0;
FragColor = vec4(0);

#ifdef PROGRESSIVE
	randSeed = bitMix(frameId);
	
	uvec2 windowPos = uvec2(gl_FragCoord.xy);
	uvec2 tilePos = windowPos >> 6;
	uint tileId = (tilePos.x << 16) | (tilePos.y & 0xffff);
	tileId = bitMix(tileId);
	uint groupCount = 16;
	tileGroup = tileId % groupCount;
	
	currentGroup = frameId % groupCount;
	
#endif

if(tileGroup == currentGroup)
{
	Ray3f ray;
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
	nextPointInUnitSphere = hashRay(ray) ^ randSeed;
	FragColor = traceSingleRay(ray, 0);
}

vec4 gamma = vec4(1.0/2.2);
FragColor = pow(FragColor, gamma);

//#ifndef PROGRESSIVE
//	FragColor.g = 0.0;
//	FragColor.b = 0.0;
//#endif

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