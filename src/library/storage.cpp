/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#include "storage.h"
#include "utility.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <unordered_set>

#ifdef CUBIQUITY_USE_POOLSTL
#define POOLSTL_STD_SUPPLEMENT
#define POOLSTL_STD_SUPPLEMENT_FORCE
#include "poolstl.hpp"
#else
#include <execution>
#endif // CUBIQUITY_USE_POOLSTL

namespace Cubiquity
{
	using namespace Internals;

	// It is desirable for the spatial relationship between child nodes to be consistent
	// throughout the tree, as this makes it easier to e.g compute normals from children
	// and to compute intersections during ray traversal. However, the use of signed integers
	// for positions results in a different spatial relationship for children of the root
	// node due to the the sign bit being set for negative numbers.
	//
	// Therefore we map the signed ints to unsigned ints. The offset is required because 
	// otherwise negative ints get mapped to large unsigned values, inverting the spatial
	// relationship between negative and positive positions.
	void make_position_unsigned(i32 ix, i32 iy, i32 iz, u32& ux, u32& uy, u32& uz)
	{
		ux = static_cast<u32>(ix) + UINT32_C(0x80000000);
		uy = static_cast<u32>(iy) + UINT32_C(0x80000000);
		uz = static_cast<u32>(iz) + UINT32_C(0x80000000);
	}

	// Potentially optimise this by bit-interleaving these three ints into a
	// larger  96-bit int (std::bitset or 2x u64) in advance of the loop calling
	// this function. Then we would just have one shift-by-three and a mask on
	// each iteration. Shorter/simpler loop body but some setup cost. SIMD on 
	// three variables could be an alternative optimisation.
	u32 extract_next_child_id(u32& ux, u32& uy, u32& uz)
	{
		u32 childX = ux & 0x80000000;
		ux = ux << 1;
		u32 childY = uy & 0x80000000;
		uy = uy << 1;
		u32 childZ = uz & 0x80000000;
		uz = uz << 1;

		return childZ >> 29 | childY >> 30 | childX >> 31;
	}

	bool Internals::isMaterialNode(const Node& node)
	{
		// Less than MaterialCount and all matching
		return isMaterialNode(node[0]) &&
			std::all_of(node.begin(), node.end(),
				[&](int x) { return x == node[0]; });
	}

	// Not allowed to be an std::hash specialisation as Node is just a typedef.
	struct Internals::NodeHasher
	{
		u64 operator()(const Internals::Node& node) const noexcept
		{
			// Reference version baed on FNV-1a
			// return Internals::fnv1a(node.data(), sizeof(node));

			// Simple XOR-multiply hash with constants picked from MurmurHash3.
			Cubiquity::u64 hash = 0;
			hash = (hash ^ node[0]) * 0x87c37b91114253d5ULL;
			hash = (hash ^ node[1]) * 0x4cf5ad432745937fULL;
			hash = (hash ^ node[2]) * 0xff51afd7ed558ccdULL;
			hash = (hash ^ node[3]) * 0xc4ceb9fe1a85ec53ULL;
			hash = (hash ^ node[4]) * 0x87c37b91114253d5ULL;
			hash = (hash ^ node[5]) * 0x4cf5ad432745937fULL;
			hash = (hash ^ node[6]) * 0xff51afd7ed558ccdULL;
			hash = (hash ^ node[7]) * 0xc4ceb9fe1a85ec53ULL;
			hash = bit_mix(hash);
			return hash;
		}
	};

	// We populate the node vector with material nodes in which the children
	// point at themselves. This has a few uses:
	// 
	//   1. It lets us recurse to a single voxel even if the tree has been
	//      pruned (used when setting a voxel).
	//   2. New material nodes can be created by copying these (rather than
	//      needing to call 'make_node(...)') which reduces conditional
	//      logic during traversal.
	//   3. During merging they can be inserted into the map in the same way as
	//      non-material nodes (as they have a valid location in the vector).
	NodeVector Internals::make_node_vector(u32 count)
	{
		// Init with invalid nodes
		Node invalid_node = make_node(InvalidNodeIndex);
		NodeVector nodeVector(count, invalid_node);

		// Set material nodes
		assert(count >= MaterialCount);
		for (u32 i = 0; i < MaterialCount; i++) {
			nodeVector[i] = make_node(i);
		}
		return nodeVector;
	}

	NodeStore::NodeStore()
		:mNodeVector(make_node_vector())
		,mSharedNodesEnd(mNodeVector.size())
	{
		/*std::mt19937 rng;
		rng.seed(42);
		std::uniform_int_distribution<u32> u32_dist;

		NodeVector test_vec;
		int count = 100000000;
		for (int i = 0; i < count; i++)
		{
			test_vec.push_back(Node{ u32_dist(rng), u32_dist(rng),  u32_dist(rng),  u32_dist(rng),  u32_dist(rng),  u32_dist(rng),  u32_dist(rng),  u32_dist(rng) });
		}

		NodeLookupSorted node_lookup(test_vec, count);
		exit(0);*/
	}

	u32 NodeStore::getNodeChild(u32 nodeIndex, u32 childId) const
	{
		if (nodeIndex < MaterialCount) { // Can bypass memory access
			return nodeIndex;
		} else {
			return mNodeVector[nodeIndex][childId];
		}
	}

	u32 NodeStore::setNodeChild(u32 nodeIndex, u32 childId, u32 value)
	{
		// Sanity checks
		//assert(!isMaterialNode(nodeIndex) && "Can't modify material nodes");
		assert(value != nodeIndex && "New child points at parent");
		assert(value != mNodeVector[nodeIndex][childId] && "Self-assignment");

		// Copy the existing node if it is shared
		if (isShared(nodeIndex)) {
			nodeIndex = append(mNodeVector[nodeIndex]);
		}

		// Update the child
		mNodeVector[nodeIndex][childId] = value;
		return nodeIndex;
	}

	u32 NodeStore::countNodes(u32 startNodeIndex) const
	{
		std::unordered_set<u32> usedIndices;
		countNodes(startNodeIndex, usedIndices);
		return usedIndices.size();
	}

	// Counts the number of nodes at distinct locations (node indices).
	// A tree which is not fully merged will contain idetical nodes at different
	// locations in memory, and these will be counted seperately.
	void NodeStore::countNodes(u32 startNodeIndex, std::unordered_set<u32>& usedIndices) const
	{
		// It may have been more efficient to have done this test before calling
		// into this function, but the implementation is simpler this way around.
		if (isMaterialNode(startNodeIndex)) { return; }

		usedIndices.insert(startNodeIndex);
		for (const u32& childNodeIndex : mNodeVector[startNodeIndex])
		{
			countNodes(childNodeIndex, usedIndices);
		}
	}

	void NodeStore::read(std::ifstream& file)
	{
		u32 nodeCount;
		file.read(reinterpret_cast<char*>(&nodeCount), sizeof(nodeCount));
		mNodeVector.resize(MaterialCount + nodeCount);
		file.read(reinterpret_cast<char*>(&(mNodeVector[sharedNodesBegin()])), sizeof(Node) * nodeCount);
		mSharedNodesEnd = mNodeVector.size(); // One past end, no unshared nodes.
	}

	void NodeStore::write(std::ofstream& file)
	{
		u32 nodeCount = sharedNodesEnd() - sharedNodesBegin();
		file.write(reinterpret_cast<const char*>(&nodeCount), sizeof(nodeCount));
		file.write(reinterpret_cast<const char*>(&mNodeVector[sharedNodesBegin()]), sizeof(Node) * nodeCount);
	}	

	u32 NodeStore::merge(u32 index)
	{
		// Destination used as hashmap and bigger than source to reduce collisions
		NodeVector merged_node_vector = make_node_vector((mNodeVector.size() * 12) / 10);

		// Perform the recursive merge
		u32 mergedRoot = merge_node(index, merged_node_vector);

		// Compact in-place, tracking how nodes move.
		std::vector<u32> moves(merged_node_vector.size());

		size_t write = 0;
		for (size_t read = 0; read < merged_node_vector.size(); ++read) {
			if (merged_node_vector[read][0] != InvalidNodeIndex) { // Occupied
				merged_node_vector[write] = std::move(merged_node_vector[read]);
				moves[read] = write;
				write++;
			}
		}
		merged_node_vector.resize(write);

		// Update node children based on how they moved during compacting.
		for (int i = 0; i < merged_node_vector.size(); i++)	{
			for (int childId = 0; childId < 8; childId++) {
				merged_node_vector[i][childId] = moves[merged_node_vector[i][childId]];
			}
		}
		mergedRoot = moves[mergedRoot];

		mNodeVector = std::move(merged_node_vector);
		mSharedNodesEnd = mNodeVector.size(); // Everything is now potentially shared

		return mergedRoot;
	}

	u32 NodeStore::merge_node(u32 nodeIndex, NodeVector& merged_node_vector)
	{
		assert(!isMaterialNode(nodeIndex));
		Node& node = mNodeVector[nodeIndex]; // By ref, modify in-place

		// In a DAG we might hit the same node more than once due to multiple 
		// parents. Handle this case to avoid reprocessing the node.
		const int Visited = 0, NewPos = 1; // Channel identifiers, not flags.
		if (node[Visited] == InvalidNodeIndex) return node[NewPos];

		// Merge the children if needed
		for (int i = 0; i < 8; i++) {
			if (!isMaterialNode(node[i])) {
				node[i] = merge_node(node[i], merged_node_vector);
			}
		}

		u32 newNodeIndex = 0;
		if (isMaterialNode(node)) {
			newNodeIndex = node[0]; // Any element works here, they are all the same
		}
		else {
			// Determine first index (skipping material node slots)
			newNodeIndex = static_cast<u32>(NodeHasher{}(node));
			newNodeIndex %= (merged_node_vector.size() - MaterialCount);
			newNodeIndex += MaterialCount;

			while (true) {
				if (merged_node_vector[newNodeIndex] == node) break; // Already present

				if (merged_node_vector[newNodeIndex][0] == InvalidNodeIndex) { // Empty space
					merged_node_vector[newNodeIndex] = node;
					break;
				}

				if (++newNodeIndex == merged_node_vector.size()) { // Move to next slot
					newNodeIndex = MaterialCount; // Wrap around, skip material nodes.
				}
			}
		}

		// Mark node as visited and note where it went. Making these changes
		// invalidates the node but we have already copied it to the destination.
		node[Visited] = InvalidNodeIndex;
		node[NewPos] = newNodeIndex;

		return newNodeIndex;
	}

	u32 NodeStore::append(const Node& node)
	{
		mNodeVector.push_back(node);
		return mNodeVector.size() - 1; // Index of new node
	}

	u32 NodeStore::cloneRoot(u32 rootIndex)
	{
		mSharedNodesEnd = mNodeVector.size(); // Mark everything as potentially shared
		mNodeVector.push_back(mNodeVector[rootIndex]); // Append unshared copy of root
		return mNodeVector.size() - 1; // Return position of new node
	}

	void NodeStore::print(std::ostream& os)
	{
		for (u32 i = 0; i < mNodeVector.size(); i++)
		{
			if (i == mSharedNodesEnd)
			{
				os << "\n-------------------- Shared nodes end / unshared nodes begin -------------------\n\n";
			}

			os << std::setw(4) << i << ": ";
			for (int c = 0; c < 8; c++)
			{
				os << std::setw(4) << mNodeVector[i][c] << " ";
			}
			os << "\n";
		}

		os << "\nShared Nodes end = " << mSharedNodesEnd << "\n";
	}

	////////////////////////////////////////////////////////////////////////////////
	// Public member functions
	////////////////////////////////////////////////////////////////////////////////

	Volume::Volume()
	{
		mRootNodeIndices.resize(1);
		mCurrentRoot = 0;
	}

	Volume::Volume(const std::string& filename)
	{
		mRootNodeIndices.resize(1);
		mCurrentRoot = 0;

		load(filename);
	}

	void Volume::fill(MaterialId matId)
	{
		setRootNodeIndex(matId);
	}

	u32 Volume::rootNodeIndex() const
	{
		return mRootNodeIndices[mCurrentRoot];
	}

	void Volume::setRootNodeIndex(u32 newRootNodeIndex)
	{
		mRootNodeIndices[mCurrentRoot] = newRootNodeIndex;
	}

	void Volume::checkpoint()
	{
		// If we 'undo' a few times and then create a checkpoint then the 'redo'
		// checkpoints are invalidated. Discard them by cropping the vector.
		mRootNodeIndices.resize(mCurrentRoot + 1); // Clear any 'redo' items

		// Clone the current root
		u32 root = mRootNodeIndices[mCurrentRoot];
		u32 newRoot = mNodeStore.cloneRoot(root);

		// Store it
		mRootNodeIndices.push_back(newRoot);
		mCurrentRoot++;
	}

	void Volume::undo()
	{
		if (mCurrentRoot > 0) {
			mCurrentRoot--;
		}
	}

	void Volume::redo()
	{
		if (mCurrentRoot < mRootNodeIndices.size() - 1)	{
			mCurrentRoot++;
		}
	}

	void Volume::bake()
	{
		Cubiquity::Timer timer;
		u32 newRootIndex = mNodeStore.merge(rootNodeIndex());
		setRootNodeIndex(newRootIndex);
		log_debug("Baked %u nodes in %.3f seconds",
			mNodeStore.unsharedNodesEnd(), timer.elapsed_seconds());
	}

	void Volume::setVoxel(i32 x, i32 y, i32 z, MaterialId matId)
	{
		// Map to unsigned positions
		u32 ux, uy, uz;
		make_position_unsigned(x, y, z, ux, uy, uz);

		const int rootHeight = log_base_2(VolumeSideLength);
		u32 newRootNodeIndex = setVoxel(ux, uy, uz, matId, rootNodeIndex(), rootHeight);
		setRootNodeIndex(newRootNodeIndex);
	}

	// I previously had iterative versions of ths function which do perform
	// better but are a little harder to reason about. One version gives the
	// same behaviour but is faster (possibly due to easier handling of the
	// early out), and another version writes nodes on the way down the tree
	// (it assumes the voxel isn't already set) rather on the way back up - this
	// is faster still. But I don't want to maintain differnt versions at the
	// moment, and the recursive version is easier. The others can be found in
	// the commit history on the NodeStore branch (present in Jan 2026).
	u32 Volume::setVoxel(u32 ux, u32 uy, u32 uz, MaterialId matId, u32 nodeIndex, int nodeHeight)
	{
		// Get the relevant child node
		u32 childId = extract_next_child_id(ux, uy, uz); // Modifies params
		u32 childNodeIndex = mNodeStore.getNodeChild(nodeIndex, childId);

		// Check if already set to desired material
		if (childNodeIndex == matId) { return nodeIndex; }

		// Recusively process the child node until we reach a just-above-leaf node (height of 1).
		if (nodeHeight > 1)
		{
			u32 newChildNodeIndex = setVoxel(ux, uy, uz, matId, childNodeIndex, nodeHeight - 1);

			// If the child hasn't changed then we don't need to update the current node.
			if (newChildNodeIndex == childNodeIndex) { return nodeIndex; }

			return mNodeStore.setNodeChild(nodeIndex, childId, newChildNodeIndex);
		}
		else
		{
			return mNodeStore.setNodeChild(nodeIndex, childId, matId);
		}
	}

	void Volume::combine(const Volume& rhsVolume, const MaterialCombiner& combiner)
	{
		u32 newRootNodeIndex = combine(rhsVolume, rhsVolume.rootNodeIndex(), rootNodeIndex(), combiner);
		setRootNodeIndex(newRootNodeIndex);
	}

	// Fixme - This function can prbably be more efficient. Firstly through better early out when whole nodes are full/empty, and also
	// if seperate volumes shared a common memory space it would be easier to copy node directly across (maybe in compressd space?).
	u32 Volume::combine(const Volume& rhsVolume, u32 rhsNodeIndex, u32 nodeIndex, const MaterialCombiner& combiner)
	{
		for (u32 childId = 0; childId < 8; childId++)
		{
			const bool nodeIsMaterial = isMaterialNode(nodeIndex);
			const bool rhsNodeIsMaterial = isMaterialNode(rhsNodeIndex);

			// If current node is a material then just propergate it. Otherwise get the true child.
			u32 childNodeIndex = mNodeStore.getNodeChild(nodeIndex, childId);
			u32 rhsChildNodeIndex = rhsVolume.mNodeStore.getNodeChild(rhsNodeIndex, childId);

			// Process children
			u32 newChildNodeIndex = childNodeIndex;

			if (isMaterialNode(childNodeIndex) && isMaterialNode(rhsChildNodeIndex))
			{
				newChildNodeIndex = combiner(childNodeIndex, rhsChildNodeIndex);
			}
			else
			{
				newChildNodeIndex = combine(rhsVolume, rhsChildNodeIndex, childNodeIndex, combiner);
			}

			// If the child has changed then we need to update the current node.
			if (childNodeIndex != newChildNodeIndex)
			{
				nodeIndex = mNodeStore.setNodeChild(nodeIndex, childId, newChildNodeIndex);
			}
		}

		return nodeIndex;
	}

	MaterialId Volume::voxel(i32 x, i32 y, i32 z) const
	{
		// TODO: Consider starting from split point instead of root
		u32 nodeIndex = rootNodeIndex();

		// Map to unsigned positions
		u32 ux, uy, uz;
		make_position_unsigned(x, y, z, ux, uy, uz);

		// Descend tree. Memory access is apparent bottleneck. Improve node
		// ordering for better cache coherency or try memory prefetch?
		while (!isMaterialNode(nodeIndex))
		{
			u32 childId = extract_next_child_id(ux, uy, uz); // Modifies params
			nodeIndex = mNodeStore[nodeIndex][childId];
		}

		return static_cast<MaterialId>(nodeIndex);
	}

	////////////////////////////////////////////////////////////////////////////////
	// Private member functions
	////////////////////////////////////////////////////////////////////////////////

	bool Volume::load(const std::string& filename)
	{
		std::ifstream file(filename, std::ios::binary);

		// FIXME - What should we do for error handling in this function? Return codes or exceptions?
		if(!file.is_open())
		{
			return false;
		}

		u32 rootNodeIndex;
		file.read(reinterpret_cast<char*>(&rootNodeIndex), sizeof(rootNodeIndex));

		mNodeStore.read(file);

		if (rootNodeIndex >= MaterialCount)
		{
			rootNodeIndex += mNodeStore.sharedNodesBegin() - MaterialCount;
		}

		setRootNodeIndex(rootNodeIndex);

		return true;
	}

	void Volume::save(const std::string& filename)
	{
		bake(); // Possibly should be here, not part of saving?

		u32 root = rootNodeIndex();
		const Node& data = mNodeStore[mNodeStore.sharedNodesBegin()];

		std::ofstream file;
		file = std::ofstream(filename, std::ios::out | std::ios::binary);
		file.write(reinterpret_cast<const char*>(&root), sizeof(root));
		mNodeStore.write(file);
		file.close();
	}

	void Volume::print(std::ostream& os)
	{
		mNodeStore.print(os);

		os << "Root nodes: [";
		for (const auto& rootNodeIndex : mRootNodeIndices)
		{
			os << rootNodeIndex << ", ";
		}
		os << "]\n";

		os << "Current root index = " << mCurrentRoot << "\n";
		os << "Current root value = " << mRootNodeIndices[mCurrentRoot] << "\n";
	}

	namespace Internals
	{
		/// This is an advanced function which should only be used if you
		/// understand the internal memory layout of Cubiquity's volume data.
		NodeStore& getNodes(Volume& volume)
		{
			return volume.mNodeStore;
		}

		const NodeStore& getNodes(const Volume& volume)
		{
			return volume.mNodeStore;
		}

		const u32 getRootNodeIndex(const Volume& volume)
		{
			return volume.rootNodeIndex();
		}
	}
}
