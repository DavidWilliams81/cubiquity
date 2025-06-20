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

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>

namespace Cubiquity
{
	using namespace Internals;

	bool Internals::isMaterialNode(u32 nodeIndex) { return nodeIndex < MaterialCount; }

	// See https://stackoverflow.com/a/57796299
	constexpr Node makeNode(u32 value)
	{
		Node node{};
		for (auto& x : node) { x = value; }
		return node;
	}

	void NodeStore::setNode(u32 index, const Node& newNode)
	{
		assert(!isMaterialNode(index) && "Error - Cannot modify material nodes");

		for (const u32& childIndex : newNode)
		{
			assert(childIndex != index && "Error - Child points at parent");
		}

		mData[index] = newNode;
	}

	void NodeStore::setNodeChild(u32 nodeIndex, u32 childId, u32 newChildIndex)
	{
		assert(!isMaterialNode(nodeIndex));
		assert(newChildIndex != nodeIndex && "Error - Node points at self");

		mData[nodeIndex][childId] = newChildIndex;
	}

	NodeDAG::NodeDAG()
	{
		mEditNodesBegin = mNodes.size();
	}

	u32 NodeDAG::countNodes(u32 startNodeIndex) const
	{
		std::unordered_set<u32> usedIndices;
		countNodes(startNodeIndex, usedIndices);
		return usedIndices.size();
	}

	// Counts the number of nodes at distinct locations (node indices).
	// A tree which is not fully merged will contain idetical nodes at different
	// locations in memory, and these will be counted seperately.
	void NodeDAG::countNodes(u32 startNodeIndex, std::unordered_set<u32>& usedIndices) const
	{
		// It may have been more efficient to have done this test before calling
		// into this function, but the implementation is simpler this way around.
		if (isMaterialNode(startNodeIndex)) { return; }

		usedIndices.insert(startNodeIndex);
		for (const u32& childNodeIndex : mNodes[startNodeIndex])
		{
			countNodes(childNodeIndex, usedIndices);
		}
	}

	void NodeDAG::read(std::ifstream& file)
	{
		u32 nodeCount;
		file.read(reinterpret_cast<char*>(&nodeCount), sizeof(nodeCount));
		for (u32 ct = 0; ct < nodeCount; ct++)
		{
			Node node;
			file.read(reinterpret_cast<char*>(&node), sizeof(node));
			mNodes.setNode(bakedNodesBegin() + ct, node);
		}
		mBakedNodesEnd = bakedNodesBegin() + nodeCount;
	}

	void NodeDAG::write(std::ofstream& file)
	{
		u32 nodeCount = bakedNodesEnd() - bakedNodesBegin();
		file.write(reinterpret_cast<const char*>(&nodeCount), sizeof(nodeCount));
		for (u32 ct = 0; ct < nodeCount; ct++)
		{
			file.write(reinterpret_cast<const char*>(&mNodes[bakedNodesBegin() + ct]), sizeof(Node));
		}
		
	}

	bool NodeDAG::isPrunable(const Node& node) const
	{
		if (!isMaterialNode(node[0])) { return false; }
		for (u32 i = 1; i < 8; i++)
		{
			if (node[i] != node[0]) { return false; }
		}

		// All children represent the same solid material, so this node can be pruned.
		return true;
	}

	void NodeDAG::merge(u32 index)
	{
		std::unordered_map<Node, u32> map;
		u32 mergedEnd = mEditNodesBegin;
		u32 nextSpace = mergedEnd - 1;

		if (isMaterialNode(index))
		{
			mBakedNodesEnd = bakedNodesBegin();
		}
		else
		{
			u32 mergedRoot = mergeNode(index, map, nextSpace);

			u32 actualNodeCount = mergedEnd - mergedRoot;

			mBakedNodesEnd = bakedNodesBegin() + actualNodeCount;

			// FIXME - This offset value seems to get large. Is the logic backwards
			// but we ar wrapping around the array so it happens to work?
			u32 offset = bakedNodesBegin() - mergedRoot;
			for (u32 nodeIndex = bakedNodesBegin(); nodeIndex < bakedNodesEnd(); nodeIndex++)
			{
				Node node = mNodes[nodeIndex - offset];
				for (u32& childIndex : node)
				{
					if (childIndex > MaxMaterial)
					{
						childIndex += offset;
					}
				}
				mNodes.setNode(nodeIndex, node);
			}
		}
	}

	u32 NodeDAG::mergeNode(u32 nodeIndex, std::unordered_map<Node, u32>& map, u32& nextSpace)
	{
		assert(!isMaterialNode(nodeIndex));
		const Node& oldNode = mNodes[nodeIndex];
		Node newNode;
		for (int i = 0; i < 8; i++)
		{
			u32 oldChildIndex = oldNode[i];
			if (!isMaterialNode(oldChildIndex))
			{
				newNode[i] = mergeNode(oldChildIndex, map, nextSpace);
			}
			else
			{
				newNode[i] = oldNode[i];
			}
		}

		u32 newNodeIndex = 0;
		auto iter = map.find(newNode);
		if (iter == map.end())
		{
			//newNodeIndex = insertNode(newNode)
			//newNodeIndex = mDAG.mInitialNodes.insert(newNode);

			mNodes.setNode(nextSpace, newNode);
			newNodeIndex = nextSpace;
			nextSpace--;

			map.insert({ newNode, newNodeIndex });
		}
		else
		{
			newNodeIndex = iter->second;
		}
		return newNodeIndex;
		//return nodes.insert(newNode);
	}

	u32 NodeDAG::insert(const Node& node)
	{
		if (mEditNodesBegin > bakedNodesEnd())
		{
			mEditNodesBegin--;
			u32 index = mEditNodesBegin;
			mNodes.setNode(index, node);
			return index;
		}

		assert(false && "Out of space for unshared edits!");
		log_warning("Out of space for unshared edits!");
		exit(1);
		return 0; // Indicates error (we don't use this function to allocate the zeroth node).
	}

	// Update the child of a node. If a copy needs to be made then this is returned,
	// otherwise the return value is empty to indicate that the update was done in-place.
	u32 NodeDAG::updateNodeChild(u32 nodeIndex, u32 childId, u32 newChildNodeIndex, bool forceCopy)
	{
		assert(newChildNodeIndex != mNodes[nodeIndex][childId]); // Watch for self-assignment (wasteful)
		assert(newChildNodeIndex != nodeIndex); // Don't let child point to parent.

		// Edit nodes can be modified in-place as they are unshared, unless the users
		// requests they are copied (e.g. for the purpose of maintaining an undo history).
		const bool modifyInPlace = isEditNode(nodeIndex) && (!forceCopy);

		if(modifyInPlace)
		{
			// Modify the existing node in-place by updating only the relevant child.
			mNodes.setNodeChild(nodeIndex, childId, newChildNodeIndex);

			// The modification may have made the node prunable. If so it must have taken on the value
			// of the new child so return that, otherwise return the updated node that we were given.
			return isPrunable(mNodes[nodeIndex]) ? newChildNodeIndex : nodeIndex;
		}
		else
		{
			const bool nodeIsMaterial = isMaterialNode(nodeIndex);

			// Make a copy of the existing node and then update the child.
			Node newNode = nodeIsMaterial ? makeNode(nodeIndex) : mNodes[nodeIndex];
			newNode[childId] = newChildNodeIndex;

			// If the copy becomes prunable as a result of the modification
			// then we can skip inserting it, which saves time and space.
			return isPrunable(newNode) ? newChildNodeIndex : insert(newNode);
		}
		
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
		if (mTrackEdits)
		{
			// When tracking edits a root node should not be set to itself. However, this
			// might happen when *not* tracking edits because we are then modifying in-place.
			assert(newRootNodeIndex != mRootNodeIndices[mCurrentRoot]);

			mCurrentRoot++;
			if (mCurrentRoot >= mRootNodeIndices.size())
			{
				mRootNodeIndices.resize(mCurrentRoot + 1);
			}
		}

		mRootNodeIndices[mCurrentRoot] = newRootNodeIndex;
	}

	void Volume::setTrackEdits(bool trackEdits)
	{
		mTrackEdits = trackEdits;
		if (!mTrackEdits)
		{
			mRootNodeIndices[0] = mRootNodeIndices[mCurrentRoot];
			mRootNodeIndices.resize(1);
			mCurrentRoot = 0;
		}
	}

	bool Volume::undo()
	{
		if (mCurrentRoot > 0)
		{
			mCurrentRoot--;
			return true;
		}

		return false; // Nothing to undo
	}

	bool Volume::redo()
	{
		if (mCurrentRoot < mRootNodeIndices.size() - 1)
		{
			mCurrentRoot++;
			return true;
		}

		return false; // Nothing to redo
	}

	void Volume::bake()
	{
		mDAG.merge(rootNodeIndex());
		setRootNodeIndex(mDAG.bakedNodesBegin());
	}

	void Volume::setVoxelRecursive(i32 x, i32 y, i32 z, MaterialId matId)
	{
		// Do we need this mapping to unsiged space? Or could we eliminate it in both voxel()/
		//setVoxel() and get the same behaviour? Or is that confusing, e.g. with raycasting?
		u32 ux = static_cast<u32>(x) ^ (1UL << 31);
		u32 uy = static_cast<u32>(y) ^ (1UL << 31);
		u32 uz = static_cast<u32>(z) ^ (1UL << 31);

		const int rootHeight = logBase2(VolumeSideLength);
		u32 newRootNodeIndex = setVoxelRecursive(ux, uy, uz, matId, rootNodeIndex(), rootHeight);
		setRootNodeIndex(newRootNodeIndex);
	}

	u32 Volume::setVoxelRecursive(u32 ux, u32 uy, u32 uz, MaterialId matId, u32 nodeIndex, int nodeHeight)
	{
		assert(nodeHeight > 0);
		u32 childHeight = nodeHeight - 1;
		
		// We could possibly remove these variable bitshifts, but I'm not sure it's worth the effort.
		u32 childX = (ux >> childHeight) & 0x01;
		u32 childY = (uy >> childHeight) & 0x01;
		u32 childZ = (uz >> childHeight) & 0x01;
		u32 childId = childZ << 2 | childY << 1 | childX;

		const bool nodeIsMaterial = isMaterialNode(nodeIndex);

		// If current node is a material then just propergate it. Otherwise get the true child.
		u32 childNodeIndex = nodeIsMaterial ? nodeIndex : mDAG[nodeIndex][childId];

		// If the child node is set to the desired material then the voxel is 
		// already set. Return invalid node indicating nothing to update.
		if (childNodeIndex == matId) { return nodeIndex; }

		// Recusively process the child node until we reach a just-above-leaf node (height of 1).
		if (nodeHeight > 1)
		{
			u32 newChildNodeIndex = setVoxelRecursive(ux, uy, uz, matId, childNodeIndex, nodeHeight - 1);

			// If the child hasn't changed then we don't need to update the current node.
			if (newChildNodeIndex == childNodeIndex) { return nodeIndex; }

			return mDAG.updateNodeChild(nodeIndex, childId, newChildNodeIndex, mTrackEdits);
		}
		else
		{
			return mDAG.updateNodeChild(nodeIndex, childId, matId, mTrackEdits);
		}
	}

	void Volume::setVoxel(i32 x, i32 y, i32 z, MaterialId matId)
	{
		//return setVoxelRecursive(x, y, z, matId);

		struct NodeState
		{
			u32 mIndex;
			bool mProcessedNode;
		};

		// Note that the first two elements of this stack never actually get used.
		// Leaf and almost-leaf nodes(heights 0 and 1) never get put on the stack.
		// We accept this wasted space, rather than subtracting two on every access.
		const int maxStackDepth = 33;
		NodeState nodeStateStack[maxStackDepth];

		const int rootHeight = logBase2(VolumeSideLength);
		int nodeHeight = rootHeight;
		nodeStateStack[nodeHeight].mIndex = rootNodeIndex();
		nodeStateStack[nodeHeight].mProcessedNode = false;

		while (true)
		{
			// This loop does not go right down to leaf nodes, it stops one level above. But for any
			// given node it manipulates it's children, which means leaf nodes can get modified.
			assert(nodeHeight >= 1);

			NodeState& nodeState = nodeStateStack[nodeHeight];
			//const Node* node = &(mDAG[nodeState.mIndex]);

			// Find which subtree we are in.
			u32 childHeight = nodeHeight - 1;
			int tx = (x ^ (1UL << 31)); // Could precalculte these.
			int ty = (y ^ (1UL << 31));
			int tz = (z ^ (1UL << 31));
			u32 childX = (tx >> childHeight) & 0x01;
			u32 childY = (ty >> childHeight) & 0x01;
			u32 childZ = (tz >> childHeight) & 0x01;
			u32 childId = childZ << 2 | childY << 1 | childX;

			if (nodeState.mProcessedNode == false) // Executed the first time we see a node - i.e. as we move *down* the tree.
			{
				const bool nodeIsMaterial = isMaterialNode(nodeState.mIndex);

				// If current node is a material then just propergate it. Otherwise get the true child.
				u32 childNodeIndex = nodeIsMaterial ? nodeState.mIndex : mDAG[nodeState.mIndex][childId];

				// If the child node is set to the desired material then the voxel is already set. 
				if (childNodeIndex == matId) { return; }

				if (nodeHeight >= 2)
				{
					NodeState& childNodeState = nodeStateStack[nodeHeight - 1];
					childNodeState.mIndex = childNodeIndex;
					childNodeState.mProcessedNode = false;

					nodeHeight -= 1;
				}

				nodeState.mProcessedNode = true;
			}
			else // Executed the second time we see a node - i.e. as we move *up* the tree.
			{
				if (nodeHeight > 1)
				{
					// If the child has changed then we need to update the current node.
					const NodeState& childNodeState = nodeStateStack[nodeHeight - 1];
					if (mDAG[nodeState.mIndex][childId] != childNodeState.mIndex)
					{
						nodeState.mIndex = mDAG.updateNodeChild(nodeState.mIndex, childId, childNodeState.mIndex, mTrackEdits);
					}
				}
				else
				{
					nodeState.mIndex = mDAG.updateNodeChild(nodeState.mIndex, childId, matId, mTrackEdits);
				}

				// Move up the tree to process parent node next, until we reach the root.
				nodeHeight += 1;
				if (nodeHeight > rootHeight)
				{
					break;
				}
			}
		}

		setRootNodeIndex(nodeStateStack[rootHeight].mIndex);
	}

	void Volume::fillBrush(const Brush& brush, MaterialId matId)
	{
		const int rootHeight = logBase2(VolumeSideLength);
		int nodeHeight = rootHeight;
		u32 newIndex = matId;

		constexpr i32 rootLowerBound = std::numeric_limits<i32>::min();

		u32 newRootNodeIndex = fillBrush(brush, matId, rootNodeIndex(), nodeHeight, rootLowerBound, rootLowerBound, rootLowerBound);
		setRootNodeIndex(newRootNodeIndex);
	}

	u32 Volume::fillBrush(const Brush& brush, MaterialId matId, u32 nodeIndex, int nodeHeight, i32 nodeLowerX, i32 nodeLowerY, i32 nodeLowerZ)
	{
		u32 childHeight = nodeHeight - 1;
		//int tx = (x ^ (1UL << 31)); // Could precalculte these.
		//int ty = (y ^ (1UL << 31));
		//int tz = (z ^ (1UL << 31));

		for (u32 childZ = 0; childZ <= 1; childZ++)
		{
			for (u32 childY = 0; childY <= 1; childY++)
			{
				for (u32 childX = 0; childX <= 1; childX++)
				{
					//u32 childX = (tx >> childHeight) & 0x01;
					//u32 childY = (ty >> childHeight) & 0x01;
					//u32 childZ = (tz >> childHeight) & 0x01;
					u32 childId = childZ << 2 | childY << 1 | childX;

					u32 childSideLength = 1 << (childHeight);
					i32 childLowerX = nodeLowerX + (childSideLength * childX);
					i32 childLowerY = nodeLowerY + (childSideLength * childY);
					i32 childLowerZ = nodeLowerZ + (childSideLength * childZ);

					i32 childUpperX = childLowerX + (childSideLength - 1);
					i32 childUpperY = childLowerY + (childSideLength - 1);
					i32 childUpperZ = childLowerZ + (childSideLength - 1);

					Box3f childBounds(vec3f({ static_cast<float>(childLowerX), static_cast<float>(childLowerY), static_cast<float>(childLowerZ) }), vec3f({ static_cast<float>(childUpperX),static_cast<float>(childUpperY), static_cast<float>(childUpperZ) }));

					if (!overlaps(brush.bounds(), childBounds))
					{
						continue;
					}

					bool allCornersInsideBrush = true;
					if (!brush.contains(vec3f({ static_cast<float>(childLowerX), static_cast<float>(childLowerY), static_cast<float>(childLowerZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childLowerX), static_cast<float>(childLowerY), static_cast<float>(childUpperZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childLowerX), static_cast<float>(childUpperY), static_cast<float>(childLowerZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childLowerX), static_cast<float>(childUpperY), static_cast<float>(childUpperZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childUpperX), static_cast<float>(childLowerY), static_cast<float>(childLowerZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childUpperX), static_cast<float>(childLowerY), static_cast<float>(childUpperZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childUpperX), static_cast<float>(childUpperY), static_cast<float>(childLowerZ) }))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f({ static_cast<float>(childUpperX), static_cast<float>(childUpperY), static_cast<float>(childUpperZ) }))) { allCornersInsideBrush = false; }

					const bool nodeIsMaterial = isMaterialNode(nodeIndex);

					// If current node is a material then just propergate it. Otherwise get the true child.
					u32 childNodeIndex = nodeIsMaterial ? nodeIndex : mDAG[nodeIndex][childId];

					// If the child node is set to the desired material then the voxel is already set.
					if (childNodeIndex == matId) { continue; }

					// Process children
					u32 newChildNodeIndex = childNodeIndex;
					//if (nodeHeight >= 2)
					if(nodeHeight >= 2 && !allCornersInsideBrush)
					{
						newChildNodeIndex = fillBrush(brush, matId, childNodeIndex, nodeHeight - 1, childLowerX, childLowerY, childLowerZ);
					}
					else
					{
						if (allCornersInsideBrush)
						{
							newChildNodeIndex = matId;
						}
					}

					// If the child has changed then we need to update the current node.
					if (childNodeIndex != newChildNodeIndex)
					{
						nodeIndex = mDAG.updateNodeChild(nodeIndex, childId, newChildNodeIndex, mTrackEdits);
					}
				}
			}
		}

		return nodeIndex;
	}

	void Volume::addVolume(const Volume& rhsVolume)
	{
		const int rootHeight = logBase2(VolumeSideLength);
		int nodeHeight = rootHeight;

		constexpr i32 rootLowerBound = std::numeric_limits<i32>::min();

		u32 newRootNodeIndex = addVolume(rhsVolume, rhsVolume.rootNodeIndex(), rootNodeIndex(), nodeHeight, rootLowerBound, rootLowerBound, rootLowerBound);
		setRootNodeIndex(newRootNodeIndex);
	}

	// Fixme - This function can prbably be more efficient. Firstly through better early out when whole nodes are full/empty, and also
	// if seperate volumes shared a common memory space it would be easier to copy node directly across (maybe in compressd space?).
	u32 Volume::addVolume(const Volume& rhsVolume, u32 rhsNodeIndex, u32 nodeIndex, int nodeHeight, i32 nodeLowerX, i32 nodeLowerY, i32 nodeLowerZ)
	{
		u32 childHeight = nodeHeight - 1;
		//int tx = (x ^ (1UL << 31)); // Could precalculte these.
		//int ty = (y ^ (1UL << 31));
		//int tz = (z ^ (1UL << 31));

		for (u32 childZ = 0; childZ <= 1; childZ++)
		{
			for (u32 childY = 0; childY <= 1; childY++)
			{
				for (u32 childX = 0; childX <= 1; childX++)
				{
					//u32 childX = (tx >> childHeight) & 0x01;
					//u32 childY = (ty >> childHeight) & 0x01;
					//u32 childZ = (tz >> childHeight) & 0x01;
					u32 childId = childZ << 2 | childY << 1 | childX;

					u32 childSideLength = 1 << (childHeight);
					i32 childLowerX = nodeLowerX + (childSideLength * childX);
					i32 childLowerY = nodeLowerY + (childSideLength * childY);
					i32 childLowerZ = nodeLowerZ + (childSideLength * childZ);

					i32 childUpperX = childLowerX + (childSideLength - 1);
					i32 childUpperY = childLowerY + (childSideLength - 1);
					i32 childUpperZ = childLowerZ + (childSideLength - 1);

					Box3f childBounds(vec3f({ static_cast<float>(childLowerX), static_cast<float>(childLowerY), static_cast<float>(childLowerZ) } ), vec3f({ static_cast<float>(childUpperX), static_cast<float>(childUpperY), static_cast<float>(childUpperZ) }));

					/*if (!overlaps(brush.bounds(), childBounds))
					{
						continue;
					}*/

					/*bool allCornersInsideBrush = true;
					if (!brush.contains(vec3f(childLowerX, childLowerY, childLowerZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childLowerX, childLowerY, childUpperZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childLowerX, childUpperY, childLowerZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childLowerX, childUpperY, childUpperZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childUpperX, childLowerY, childLowerZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childUpperX, childLowerY, childUpperZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childUpperX, childUpperY, childLowerZ))) { allCornersInsideBrush = false; }
					if (!brush.contains(vec3f(childUpperX, childUpperY, childUpperZ))) { allCornersInsideBrush = false; }*/

					const bool nodeIsMaterial = isMaterialNode(nodeIndex);
					const bool rhsNodeIsMaterial = isMaterialNode(rhsNodeIndex);

					// If current node is a material then just propergate it. Otherwise get the true child.
					u32 childNodeIndex = nodeIsMaterial ? nodeIndex : mDAG[nodeIndex][childId];
					u32 rhsChildNodeIndex = rhsNodeIsMaterial ? rhsNodeIndex : rhsVolume.mDAG[rhsNodeIndex][childId];

					// If the child node is set to the desired material then the voxel is already set.
					if (isMaterialNode(childNodeIndex) && isMaterialNode(rhsChildNodeIndex) && (childNodeIndex == rhsChildNodeIndex)) { continue; }

					// Process children
					u32 newChildNodeIndex = childNodeIndex;

					// We have to treat at least one material as empty space, otherwise we are just copying 
					// every voxel from source to destination, which just result in a copy of the source. 
					// Hard-code it to zero for now, but should think about how else it might be controlled.
					MaterialId emptyMaterial = 0;
					if (rhsChildNodeIndex != emptyMaterial)
					{
						if (isMaterialNode(rhsChildNodeIndex))
						{
							newChildNodeIndex = rhsChildNodeIndex;
						}
						else
						{
							newChildNodeIndex = addVolume(rhsVolume, rhsChildNodeIndex, childNodeIndex, nodeHeight - 1, childLowerX, childLowerY, childLowerZ);;
						}
					}

					// If the child has changed then we need to update the current node.
					if (childNodeIndex != newChildNodeIndex)
					{
						nodeIndex = mDAG.updateNodeChild(nodeIndex, childId, newChildNodeIndex, mTrackEdits);
					}
				}
			}
		}

		return nodeIndex;
	}

	MaterialId Volume::voxel(i32 x, i32 y, i32 z) const
	{
		u32 nodeIndex = rootNodeIndex();
		u32 height = logBase2(VolumeSideLength);

		// FIXME - think whether we need the line below - I think we do for empty/solid volumes?
		//if (mDAG.isMaterialNode(mRootNodeIndex)) { return static_cast<MaterialId>(mRootNodeIndex); }

		while (height >= 1)
		{			
			// If we reach a full node then the requested voxel is occupied.
			if (isMaterialNode(nodeIndex)) { return static_cast<MaterialId>(nodeIndex); }

			// Otherwise find which subtree we are in.
			// Optimization - Note that the code below requires shifting by a variable amount which can be slow.
			// Alternatively I think we can simply shift x, y, and z by one bit per iteration, but this requires us
			// to reverse the order of the bits at the start of this function. It would be a one-time cost for a
			// faster loop, and testing is needed (on a real, large volume) to determine whether it is beneficial.
			u32 childHeight = height - 1;
			int tx = (x ^ (1UL << 31)); // Could precalculte these.
			int ty = (y ^ (1UL << 31));
			int tz = (z ^ (1UL << 31));
			u32 childX = (tx >> childHeight) & 0x01;
			u32 childY = (ty >> childHeight) & 0x01;
			u32 childZ = (tz >> childHeight) & 0x01;
			u32 childId = childZ << 2 | childY << 1 | childX;

			// Prepare for next iteration.
			nodeIndex = mDAG[nodeIndex][childId];
			height--;
		}

		// We have reached a height of zero so the node must be a material node.
		assert(height == 0 && isMaterialNode(nodeIndex));
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
		//assert(rootNodeIndex == mDAG.arrayBegin());
		//setRootNodeIndex(RefCountedNodeIndex(rootNodeIndex, &mDAG));

		//mDAG.read(file);

		// This fill is not required but can be useful for debugging
		//const Node InvalidNode = makeNode(0xffffffff);
		//std::fill(mDAG.mNodes.begin(), mDAG.mNodes.end(), InvalidNode);

		mDAG.read(file);

		if (rootNodeIndex >= MaterialCount)
		{
			rootNodeIndex += mDAG.bakedNodesBegin() - MaterialCount;
		}

		setRootNodeIndex(rootNodeIndex);

		return true;
	}

	void Volume::save(const std::string& filename)
	{
		bake();

		u32 root = rootNodeIndex();
		const Node& data = mDAG[mDAG.bakedNodesBegin()];

		std::ofstream file;
		file = std::ofstream(filename, std::ios::out | std::ios::binary);
		file.write(reinterpret_cast<const char*>(&root), sizeof(root));
		mDAG.write(file);
		file.close();
	}

	namespace Internals
	{
		/// This is an advanced function which should only be used if you
		/// understand the internal memory layout of Cubiquity's volume data.
		NodeDAG& getNodes(Volume& volume)
		{
			return volume.mDAG;
		}

		const NodeDAG& getNodes(const Volume& volume)
		{
			return volume.mDAG;
		}

		const u32 getRootNodeIndex(const Volume& volume)
		{
			return volume.rootNodeIndex();
		}
	}
}
