#include "raytracing.h"

#include <cfloat>
#include <climits>

#define COMPILE_AS_CPP 1 // Not defined in GLSL version

namespace Cubiquity
{
	// The full DAG has 33 levels, from zero (for leaves) to 32 (for the root).
	const int RootNodeHeight = 32;

	ivec3 childIdToIVec3(int childId)
	{
		return ivec3({(childId) & 0x1, (childId >> 1) & 0x1, (childId >> 2) & 0x1});
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	//  ______              _                  _                                                  //
	//  | ___ \            | |                (_)                                                 //
	//  | |_/ /__ _ _   _  | |_ _ __ __ _  ___ _ _ __   __ _                                      //
	//  |    // _` | | | | | __| '__/ _` |/ __| | '_ \ / _` |                                     //
	//  | |\ \ (_| | |_| | | |_| | | (_| | (__| | | | | (_| |                                     //
	//  \_| \_\__,_|\__, |  \__|_|  \__,_|\___|_|_| |_|\__, |                                     //
	//               __/ |                              __/ |                                     //
	//              |___/                              |___/                                      //
	//                                                                                            //
	////////////////////////////////////////////////////////////////////////////////////////////////

	SubDAG findSubDAG(const Internals::NodeStore& nodes, uint rootNodeIndex, uint childId)
	{
		// Initialised for root, but updated on first iteration of the loop.
		int childHeight = 32;
		ivec3 lowerBound = { INT_MIN , INT_MIN , INT_MIN };

		// We never return the root node as a subDAG, instead
		// we always descend into at least the first child.
		uint onlyChildId = childId;
		uint nextNodeIndex = nodes[rootNodeIndex][onlyChildId];
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
				uint childNodeIndex = nodes[nodeIndex][i];
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

	SubDAGArray findSubDAGs(const Internals::NodeStore& nodes, uint32 rootNodeIndex)
	{
		SubDAGArray subDAGs;
		for (uint childId = 0; childId < 8; childId++)
		{
			subDAGs[childId] = findSubDAG(nodes, rootNodeIndex, childId);
		}
		return subDAGs;
	}

	SubDAG getSubDAG(const Internals::NodeStore& nodes, uint rootNodeIndex, const SubDAGArray& subDAGs, uint childId)
	{
		int method = 1;
		switch (method) {
		case 1:
			// Use the real precomputed SubDAG
			return subDAGs[childId];
		case 2:
			// Compute the SubDAG on demand
			return findSubDAG(nodes, rootNodeIndex, childId);
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
			bypassSubDAG.nodeIndex = nodes[rootNodeIndex][childId];
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
	uint findNearestMaterial(const Internals::NodeStore& nodes, uint nodeIndex, uint rayDirSignBits)
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
				uint childNodeIndex = nodes[nodeIndex][childIds & 0x7];
				if (childNodeIndex > 0) // Skip empty nodes
				{
					nodeIndex = childNodeIndex;
					break;
				}
			}
		}

		return nodeIndex;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	// See "An Efficient Parametric Algorithm for Octree Traversal"
	///////////////////////////////////////////////////////////////////////////////////////////////////

	float min3(vec3 vec)
	{
		return std::min(std::min(vec.x(), vec.y()), vec.z());
	}

	float max3(vec3 vec)
	{
		return std::max(std::max(vec.x(), vec.y()), vec.z());
	}

	ivec3 findFirstChild(float nodeEntry, vec3 rayOrigin, vec3 invRayDir, ivec3 nodeCentre)
	{
		// Distances to midpoint along each axis.
		vec3 nodeTm = (vec3(nodeCentre) - rayOrigin) * invRayDir;

		// Find first child by determining the axes along which the ray has already passed the midpoint.
		// Logic reversed compared to ESVO paper, probabaly because they have a negative ray direction?
		ivec3 childId = ivec3(lessThan(nodeTm, vec3({ nodeEntry, nodeEntry, nodeEntry })));

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
	RayVolumeIntersection intersectRayNodeESVO(const Internals::NodeStore& nodes,
		uint nodeIndex, ivec3 nodePos, int nodeHeight,
		Ray3f ray, vec3 rayDirSign, uint rayDirSignBits,
		bool computeSurfaceProperties, float maxFootprint)
	{
		RayVolumeIntersection intersection = { false, 0, 0, {0, 0, 0}, {0, 0, 0} }; // Miss

		uint nodeSize = 1u << nodeHeight;

		// Intersect ray wth the node. The cast to float before
		// addition prevents integer overflow for largest nodes.
		vec3 invRayDir = vec3({ 1,1,1 }) / ray.mDir;
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
				uint childNodeIndex = nodes[nodeIndex][childIdBits ^ rayDirSignBits];

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
							intersection.material = findNearestMaterial(nodes, childNodeIndex, rayDirSignBits);

							// Determining a true voxel normal is tricky because we don't have easy access to neighbouring
							// voxels. Therefore we currently determine a face normal instead by looking at which side of
							// the cube the ray entered. This is fine (even desirable, depending on user requirements) and
							// we can probably smooth out voxels in screen space if we want to. But I wonder if a true 
							// voxel normal would be better when bouncing rays around inside the volume for pathtracing.
							intersection.normal = vec3(equal(vec3({ tChildEntry, tChildEntry, tChildEntry }), childT0));
							intersection.normal *= -rayDirSign;
						}
					}
				}
				else
				{

					// Advance forward to the next the child node (sibling). The ESVO paper uses an equality
					// test in the text but '<=' in the sample code. I don't know why (though I can see it
					// is valid) but perhaps it adds robustness in the face of floating point inaccuracy?
					ivec3 nextChildFlipsVec = ivec3(lessThanEqual(childT1, vec3({ tChildExit, tChildExit, tChildExit })));
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
	RayVolumeIntersection intersectVolume(const Volume& volume, const SubDAGArray& subDAGs, Ray3f ray, bool computeSurfaceProperties, float maxFootprint)
	{
		RayVolumeIntersection intersection = { false, 0, 0, {0, 0, 0}, {0, 0, 0} }; // Miss

		const Internals::NodeStore& nodes = Internals::getNodes(volume).nodes();

		const uint rootNodeIndex = Internals::getRootNodeIndex(volume);

		// Store the sign of the ray direction
		ivec3 rayDirSignBitsAsVec = ivec3(lessThan(ray.mDir, vec3({0,0,0}))); // 0 means +ve, 1 means -ve
		uint rayDirSignBits = rayDirSignBitsAsVec[0] | (rayDirSignBitsAsVec[1] << 1) | (rayDirSignBitsAsVec[2] << 2);
		vec3 rayDirSign = vec3(rayDirSignBitsAsVec * (-2) + 1); // Binary (-1, 1) not ternary (-1, 0, 1)

		// We reflect the ray so that it points along the positive axes. This simplifies the traversal.
		// Shifting the ray origin lets us pretend voxel corners have integer coordinates (not +/-0.5).
		Ray3f reflectedRay = ray;
		reflectedRay.mOrigin += vec3({ 0.5f, 0.5f, 0.5f });
		reflectedRay.mOrigin = reflectedRay.mOrigin * rayDirSign;
		reflectedRay.mDir = abs3(reflectedRay.mDir);

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
			SubDAG subDAG = getSubDAG(nodes, rootNodeIndex, subDAGs, childId ^ rayDirSignBits);
			uint childNodeIndex = subDAG.nodeIndex;

			if (childNodeIndex > 0)
			{
				// FIXME - Handle fully occupied direct child of root node.
				assert(childNodeIndex >= MaterialCount);

				// We can store the lower bound as integers because we offset ray origin by 0.5.
				// On reflected axis the 'upper' corner becomes lower than the 'lower' corner.
				// Node size might be large and negative, but I *think* the logic still works.
				int nodeSize = int(1u << subDAG.nodeHeight);
				ivec3 refNodeLowerBound = subDAG.lowerBound * ivec3(rayDirSign);
				refNodeLowerBound -= rayDirSignBitsAsVec * ivec3({ nodeSize, nodeSize, nodeSize });

				intersection = intersectRayNodeESVO(nodes,
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
}
