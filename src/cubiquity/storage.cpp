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
#include <iostream>

namespace Cubiquity
{
    const int ArraySizePower = 24;
    const int ArraySize = 1 << ArraySizePower;

	using namespace Internals;

	// Internal utility functions
	namespace Internals
	{
		// Node is initially empty
		Node::Node(uint32_t initialChildValue)
		{
			mRefCount = 0;

			for (int i = 0; i < 8; i++)
			{
				mChildren[i] = initialChildValue;
			}
		}

		// Note: It is perhaps misleading to use this as the equality operator as it ignores
		// the reference count... think if that matters? Should this function be renamed?
		bool Node::operator== (const Node& rhs) const
		{
			for (int i = 0; i < 8; i++)
			{
				if (mChildren[i] != rhs.mChildren[i])
				{
					return false;
				}
			}
			return true;
		}

		bool Node::operator!= (const Node& rhs) const
		{
			return !(*this == rhs);
		}

		uint32_t Node::refCount() const
		{
			return mRefCount;
		}

		uint32_t Node::child(uint32_t id) const
		{
			assert(id < 8);
			return mChildren[id];
		}

		void Node::setChild(uint32_t id, uint32_t value)
		{
			assert(id < 8);
			mChildren[id] = value;
		}

		bool Node::allChildrenAre(uint32_t value) const
		{
			for (int i = 0; i < 8; i++)
			{
				if (mChildren[i] != value)
				{
					return false;
				}
			}
			return true;
		}

		bool Node::allChildrenAreLessThan(uint32_t value) const
		{
			for (int i = 0; i < 8; i++)
			{
				if (mChildren[i] >= value)
				{
					return false;
				}
			}
			return true;
		}
	}

	NodeArray::NodeArray(uint32_t size, MaterialId initialValue)
		: mSize(0), mData(nullptr)
	{
		// Validate the minimum size (though such
		// a tiny array is useless in practice).
		const uint32_t minSize = RootNodeIndex + 1;
		assert(size >= minSize);
		size = std::max(size, minSize);

		// Allocate the data
		mSize = size;
		mData = new Node[mSize];

		// The first nodes are dummy data and should never be touched.
		// Index values of '0' and '1' represent empty and full nodes
		// respectively, and should not be used for array lookups.
		memset(mData, 0xFF, sizeof(Node) * MaterialCount);

		// Set up the root node with the initial values.
		mData[RootNodeIndex].mRefCount = 0;
		for (auto& child : mData[RootNodeIndex].mChildren)
		{
			child = initialValue;
		}

		// The rest of the array is zero-initiialized.
		const uint32_t remainingNodeStart = RootNodeIndex + 1;
		memset(&(mData[remainingNodeStart]), 0, sizeof(Node) * (mSize - remainingNodeStart));

		// We consider the root node to have already been allocated.
		// Subsequent allocations carry on from there.
		mLastAllocation = RootNodeIndex;
	}

	NodeArray::NodeArray(const NodeArray& other)
	{
		mSize = other.mSize;
		mData = new Node[mSize];
		memcpy(mData, other.mData, sizeof(Node) * mSize);
		mLastAllocation = other.mLastAllocation;
	}

	NodeArray::~NodeArray()
	{
		// Free the data
		delete[] mData;
	}

	NodeArray& NodeArray::operator=(const NodeArray& other)
	{
		mSize = other.mSize;
		mData = new Node[mSize];
		memcpy(mData, other.mData, sizeof(Node) * mSize);
		mLastAllocation = other.mLastAllocation;
		return *this;
	}

	uint32_t NodeArray::allocateNode(uint32_t initialChildValue)
	{
		Node defaultNode(initialChildValue);

		// Search from the last allocation to the end of the array
		for (uint32_t i = mLastAllocation + 1; i < mSize; i++)
		{
			if (mData[i].mRefCount == 0)
			{
				assert(i > 2); // Nodes 0, 1, and 2 are reserved
				mLastAllocation = i;
				mData[i] = defaultNode;
				return i;
			}
		}

		// No luck so far, so try searching the first part of the array
		// Note that we skip the reserved elements
		for (uint32_t i = RootNodeIndex + 1; i <= mLastAllocation; i++) 
		{
			if (mData[i].mRefCount == 0)
			{
				assert(i > 2); // Nodes 0, 1, and 2 are reserved
				mLastAllocation = i;
				mData[i] = defaultNode;
				return i;
			}
		}

		assert(false);
		std::cout << "Out of memory!" << std::endl;

		// Needs a better way of indicating this?
		return 0; // Indicates error (we don't use this function to allocate the zeroth node).
	}

	void NodeArray::resize(uint32_t size)
	{
		Node* newData = new Node[size];
		memcpy(newData, mData, size * sizeof(Node));
		delete[] mData;

		mData = newData;
		mSize = size;
	}

	uint32_t NodeArray::size() const
	{
		return mSize;
	}

	void NodeArray::setChild(uint32_t nodeIndex, uint32_t childId, uint32_t newChildIndex)
	{
		//  Can't set the reserved nodes.
		assert(nodeIndex >= MaterialCount);

		// A node cannot point at itself.
		assert(newChildIndex != nodeIndex);

		uint32_t oldChildIndex = mData[nodeIndex].child(childId);

		// FIXME - Need to think carefully what happens if we set a node to it's current value.
		// We can make a mess if we decrease the ref count and delete it!
		//assert(newChildIndex != oldChildIndex);

		if (newChildIndex == oldChildIndex) // Should we let this happen?
		{
			// Nothing to do?
			return;
		}

		if (oldChildIndex >= MaterialCount)
		{
			assert(mData[oldChildIndex].mRefCount > 0);
			mData[oldChildIndex].mRefCount--;

			if (mData[oldChildIndex].mRefCount == 0)
			{
				for (int i = 0; i < 8; i++)
				{
					setChild(oldChildIndex, i, 0); // FIXME - What should we set to here?
				}
			}
		}

		mData[nodeIndex].setChild(childId, newChildIndex);

		if (newChildIndex >= MaterialCount)
		{
			mData[newChildIndex].mRefCount++;
		}
	}

	const Internals::Node& NodeArray::operator[](uint32_t index) const
	{
		assert(index < mSize);

		return mData[index];
	}

	/*void NodeArray::mergeOctree()
	{
	const int hashTableLength = 0xFFFF;
	uint32_t hashTable[hashTableLength];

	//FIXME - How should this be initialised? Need some way of marking unfilled entries?
	//memset(hashTable, 0xff, hashTableLength * sizeof(uint32_t));
	for(int h = 0; h < hashTableLength; h++)
	{
	hashTable[h] = 0x12345678;
	}

	Internals::NodeHash nodeHasher;

	//Should this loop start at 0?
	for (uint32_t i = 0; i < mSize; i++)
	{
	Node& node = mData[i];
	if (node.mRefCount > 0)
	{
	size_t hash = nodeHasher(node);

	// If we've already got a node for this hash then don't overwrite
	// it. This shold let us keep element '1' in the table?
	if (hashTable[hash % hashTableLength] == 0x12345678)
	{
	hashTable[hash % hashTableLength] = i;
	}
	}
	}

	// Should this loop start at 0?
	for (uint32_t index = 2; index < mSize; index++)
	{
	Node& node = mData[index];
	if (node.mRefCount > 0)
	{
	for (int i = 0; i < 8; i++)
	{
	uint32_t childIndex = node.child(i);

	const Node& childNode = mData[childIndex];

	size_t childHash = nodeHasher(childNode);

	uint32_t potentialMatchIndex = hashTable[childHash % hashTableLength];

	if (childNode == mData[potentialMatchIndex])
	{
	setChild(index, i, potentialMatchIndex);
	}
	//else
	//{
	//	hashTable[childHash % hashTableLength] = childIndex;
	//}
	}
	}
	}
	}*/

	void NodeArray::mergeOctree()
	{
		const int hashTableLength = 0xFFFFFF;
		uint32_t* hashTable = new uint32_t[hashTableLength];

		Internals::NodeHash nodeHasher;

		int passIndex = 0;

		bool somethingChanged = false;
		do
		{
			passIndex++;

			somethingChanged = false;

			uint32_t lastValidNodeIndex = RootNodeIndex;
			for (uint32_t i = RootNodeIndex; i < mSize; i++)
			{
				Node& node = mData[i];
				if (node.mRefCount > 0)
				{
					lastValidNodeIndex = i;
				}
			}

			//FIXME - How should this be initialised? Need some way of marking unfilled entries?
			//memset(hashTable, 0xff, hashTableLength * sizeof(uint32_t));
			for (int h = 0; h < hashTableLength; h++)
			{
				hashTable[h] = 0x12345678;
			}

			// We skip the root node because no other node ever points at it.
			/*for (uint32_t i = RootNodeIndex + 1; i <= lastValidNodeIndex; i++)
			{
			Node& node = mData[i];
			if (node.mRefCount > 0)
			{
			size_t hash = nodeHasher(node);
			assert(hash != 0x12345678); //This is our uninitialised value... how should we handle this?

			uint32_t hashTableIndex = hash % hashTableLength;

			if (hashTable[hashTableIndex] == 0x12345678)
			{
			hashTable[hashTableIndex] = i;
			}
			}
			}*/

			for (uint32_t innerIndex = RootNodeIndex; innerIndex <= lastValidNodeIndex; innerIndex++)
			{
				Node& innerNode = mData[innerIndex];
				// We do potentially merge the children of the root node,
				// even though the root node itself has a ref count of zero.
				if (innerNode.mRefCount == 0 && innerIndex != RootNodeIndex) continue;

				for (int i = 0; i < 8; i++)
				{
					uint32_t innerChildIndex = innerNode.child(i);

					if (innerChildIndex < RootNodeIndex) continue;

					const Node& innerChildNode = mData[innerChildIndex];

					size_t innerChildNodeHash = nodeHasher(innerChildNode, passIndex);
					assert(innerChildNodeHash != 0x12345678); //This is our uninitialised value... how should we handle this?

					uint32_t potentialMatchIndex = hashTable[innerChildNodeHash % hashTableLength];

					if (potentialMatchIndex == 0x12345678)
					{
						hashTable[innerChildNodeHash % hashTableLength] = innerChildIndex;
					}
					else
					{
						const Node& potentialMatchNode = mData[potentialMatchIndex];

						if (innerChildIndex != potentialMatchIndex)
						{
							if (innerChildNode == potentialMatchNode)
							{
								setChild(innerIndex, i, potentialMatchIndex);
								somethingChanged = true;
							}
						}
					}
				}
			}

		} while (somethingChanged);

		delete[] hashTable;
	}

	void NodeArray::mergeOctreeBruteForce()
	{
		bool somethingChanged = false;
		do
		{
			somethingChanged = false;

			uint32_t lastValidNodeIndex = RootNodeIndex;
			for (uint32_t i = RootNodeIndex; i < mSize; i++)
			{
				Node& node = mData[i];
				if (node.mRefCount > 0)
				{
					lastValidNodeIndex = i;
				}
			}

			// We skip the root node because no other node ever points at it.
			for (uint32_t outerIndex = RootNodeIndex + 1; outerIndex <= lastValidNodeIndex; outerIndex++)
			{
				Node& outerNode = mData[outerIndex];
				if (outerNode.mRefCount == 0) continue;

				for (uint32_t innerIndex = RootNodeIndex; innerIndex <= lastValidNodeIndex; innerIndex++)
				{
					Node& innerNode = mData[innerIndex];
					// We do potentially merge the children of the root node,
					// even though the root node itself has a ref count of zero.
					if (innerNode.mRefCount == 0 && innerIndex != RootNodeIndex) continue;

					for (int i = 0; i < 8; i++)
					{
						uint32_t innerChildIndex = innerNode.child(i);

						// The following line should have made it faster but it doesn't?
						if (innerChildIndex < RootNodeIndex) continue;

						const Node& innerChildNode = mData[innerChildIndex];

						if (innerChildIndex != outerIndex)
						{
							if (innerChildNode == outerNode)
							{
								setChild(innerIndex, i, outerIndex);
								somethingChanged = true;
							}
						}
					}
				}
			}
		} while (somethingChanged);
	}

	// This function could be significantly improved. Currently it copies one node at a time,
	// but memmove (not memcpy) could be used to shift large groups of nodes down at a time.
	// The map array could also store offsets instead of absolute positions, and it could then
	// be efficeiently stored with run-length encoding which would probably decrease memory usage.
	void NodeArray::defragment()
	{
		uint32_t* map = new uint32_t[size()];
		memset(map, 0xFF, sizeof(uint32_t) * size());

		uint32_t newPos = RootNodeIndex + 1;
		for (uint32 i = RootNodeIndex + 1; i < size(); i++)
		{
			Node node = mData[i];
			if (node.mRefCount > 0)
			{
				if (newPos != i)
				{
					mData[newPos] = node;
					memset((&mData[i]), 0, sizeof(mData[i]));
				}
				map[i] = newPos;
				newPos++;
			}
		}

		// Now fix up the child indices to point to the new locations in the defragmented array.
		// FIXME - Isn't it enough to iterate up to the number of copied nodes here, rather than the full length of the array?
		for (uint32 n = RootNodeIndex; n < size(); n++)
		{
			Node& node = mData[n];

			for (int i = 0; i < 8; i++)
			{
				//if (compactedNode.refCount() > 0)
				{
					if (node.child(i) >= MaterialCount)
					{
						uint32_t index = node.child(i);
						uint32_t mappedIndex = map[index];
						node.setChild(i, mappedIndex);
					}
				}
			}
		}

		delete[] map;
	}

	void NodeArray::read(std::ifstream& file)
	{
		uint32_t nodeCount;
		file.read(reinterpret_cast<char*>(&nodeCount), sizeof(nodeCount));

		// Rather than resizing, we should be using a fixed size array here and we should memset() the remainder to zero.
		//resize(RootNodeIndex + nodeCount);
		// The rest of the array is zero-initiialized.
		const uint32_t remainingNodeStart = RootNodeIndex + 1;
		memset(&(mData[remainingNodeStart]), 0, sizeof(Node) * (mSize - remainingNodeStart));
		file.read(reinterpret_cast<char*>(&mData[RootNodeIndex]), nodeCount * sizeof(mData[0]));
	}

	void NodeArray::write(std::ofstream& file)
	{
		defragment();

		// Count the number of nodes to write out. After defragmenting the valid nodes are all consecutive in the.
		// array. So the end of the valid data is marked by the first node with a ref count of zero. Note that
		// we start counting one past the root node, as the root node actually has a ref count of zero as well.
		uint32_t nodeCount = 1;
		for (uint32_t i = RootNodeIndex + 1; i < mSize && mData[i].mRefCount > 0; i++)
		{
			nodeCount++;
		}

		file.write(reinterpret_cast<const char*>(&nodeCount), sizeof(nodeCount));
		file.write(reinterpret_cast<const char*>(&mData[RootNodeIndex]), nodeCount * sizeof(mData[0]));
	}

	////////////////////////////////////////////////////////////////////////////////
	// Public member functions
	////////////////////////////////////////////////////////////////////////////////

	Volume::Volume(MaterialId initialValue)
		:mNodeArray(ArraySize, initialValue)
	{
	}

	Volume::Volume(const std::string& filename)
		:mNodeArray(ArraySize, 0)
	{
		load(filename);
	}

	void Volume::compact()
	{
		mNodeArray.mergeOctree();
		mNodeArray.defragment();
	}

	void Volume::setVoxel(int32_t x, int32_t y, int32_t z, MaterialId matId, bool checkIfAlreadySet)
	{
		struct NodeState
		{
			NodeState()
				: mIndex(0)
				, mProcessedNode(false) {}

			uint32_t mIndex;
			bool mProcessedNode;
		};

		// Setting a voxel can be relatively expensive (compared to checking one) because
		// we may have to spend time unsharing nodes. The implementation is also simpler
		// and more efficient if we know that a value is actually changing, otherwise we
		// may spend time unsharing as we work down the tree, only to find that the value
		// doesn't need to change and so all the unsharing was unnecessary. Therefore we
		// perform the check here.
		//
		// It should be noted that both this check and the following set involve traversing
		// the tree. It would be nice to only do this traversal once, but I haven't yet
		// come up with an elegant way to do it.
		//
		// It would feel nice if we could traverse the tree to the required voxel, check
		// if it needs changing, and if so do the unsharing on the way back up. But I
		// think the unsharing needs to be done in a top-down order, otherwise you are
		// changing nodes (children) which haven't yet been unshared themselves. Besides,
		// is a top-down followed by bottom up traversal really any better than two top
		// down traversals? Perhaps not, so maybe the current approach is ok.
		//
		// Overall I don't know whether or not we should perform this check. It feels ugly,
		// in the sense that it should be perfectly safe to set a voxel to it's current
		// value and so I don't like checking for it. However, because it can cause unsharing
		// to occur, setting a voxel to it's current value causes no observable change to the
		// volume (except to it's memory usage) but does degrade its internal state.

		if (checkIfAlreadySet)
		{
			if (voxel(x, y, z) == matId)
			{
				return; // Already set, nothing to do.
			}
		}

		// Note that the first two elements of this stack never actually get used.
		// Leaf and almost-leaf nodes(heights 0 and 1) never get put on the stack.
		// We accept this wasted space, rather than subtracting two on every access.
		const int maxStackDepth = 33;
		NodeState nodeStateStack[maxStackDepth];

		const int rootHeight = logBase2(VolumeSideLength);
		int nodeHeight = rootHeight;
		nodeStateStack[nodeHeight].mIndex = RootNodeIndex;

		while (true)
		{
			// This loop does not go right down to leaf ndes, it stops one level above. But for any
			// given node it manipulates it's children, which means leaf nodes can get modified.
			assert(nodeHeight >= 1);

			NodeState& nodeState = nodeStateStack[nodeHeight];
			const Node* node = &(mNodeArray[nodeState.mIndex]);

			// Find which subtree we are in.
			uint32_t childHeight = nodeHeight - 1;
			int tx = (x ^ (1UL << 31)); // Could precalculte these.
			int ty = (y ^ (1UL << 31));
			int tz = (z ^ (1UL << 31));
			uint32_t childX = (tx >> childHeight) & 0x01;
			uint32_t childY = (ty >> childHeight) & 0x01;
			uint32_t childZ = (tz >> childHeight) & 0x01;
			uint32_t childId = childZ << 2 | childY << 1 | childX;

			if (nodeState.mProcessedNode == false) // Executed the first time we see a node - i.e. as we move *down* the tree.
			{
				nodeState.mProcessedNode = true;


				if (nodeHeight >= 2)
				{
					if (node->child(childId) == matId)
					{
						return;
					}

					auto current = node->child(childId);
					if (current < MaterialCount)
					{
						uint32_t i = mNodeArray.allocateNode(current);
						mNodeArray.setChild(nodeState.mIndex, childId, i);
					}

					const Node* childNode = &(mNodeArray[node->child(childId)]);

					if (childNode->refCount() > 1)
					{
						//mNodeArray.unshareChildOfNode(nodeState.mIndex, childIndex);
						splitChild(nodeState.mIndex, childId);
						childNode = &(mNodeArray[node->child(childId)]);
						assert(childNode->refCount() == 1);
					}

					NodeState& childNodeState = nodeStateStack[nodeHeight - 1];
					childNodeState.mIndex = node->child(childId);
					childNodeState.mProcessedNode = false;
					nodeHeight -= 1;
				}
				else // nodeHeight == 1
				{
					// We've reached a just-above-leaf node, so set the approriate child (leaf) and we are done.

					assert(nodeHeight == 1);

					// Note: Do we want a check here for whether it is already set?
					mNodeArray.setChild(nodeState.mIndex, childId, matId);
				}

			}
			else // Executed the second time we see a node - i.e. as we move *up* the tree.
			{
				// We have deleted a voxel by setting one of the children of a node to be invalid.
				// If the other children were already invalid, then all the children are now
				// invalid and it looks like the node is full. This is not true (it is actually
				// empty) and so we must delete the 'full' node.
				//
				// Note that we cannot delete the current node as we would need to know who its
				// parent is and we don't have this information. So we implement the logic by
				// looking at the child of the current node and deleting that if it is full.

				// Do we need this test?
				if (node->child(childId) >= MaterialCount)
				{
					const Node& childNode = mNodeArray[node->child(childId)];
					if (childNode.allChildrenAre(matId))
					{
						mNodeArray.setChild(nodeState.mIndex, childId, matId);
					}
				}

				// If we get this far then the current node has no more children to process, so
				// we are done with it. Pop it's parent off the stack and carry on processing that.
				nodeHeight += 1;
				if (nodeHeight > rootHeight)
				{
					break;
				}
			}
		}
	}

	MaterialId Volume::voxel(int32_t x, int32_t y, int32_t z)
	{
		uint32_t nodeIndex = RootNodeIndex;
		uint32_t height = logBase2(VolumeSideLength);

		while (height >= 1)
		{			
			// If we reach a full node then the requested voxel is occupied.
			//if (mNodeArray[nodeIndex].isFull()) return true;

			// Otherwise find which subtree we are in.
			// Optimization - Note that the code below requires shifting by a variable amount which can be slow.
			// Alternatively I think we can simply shift x, y, and z by one bit per iteration, but this requires us
			// to reverse the order of the bits at the start of this function. It would be a one-time cost for a
			// faster loop, and testing is needed (on a real, large volume) to determine whether it is beneficial.
			uint32_t childHeight = height - 1;
			int tx = (x ^ (1UL << 31)); // Could precalculte these.
			int ty = (y ^ (1UL << 31));
			int tz = (z ^ (1UL << 31));
			uint32_t childX = (tx >> childHeight) & 0x01;
			uint32_t childY = (ty >> childHeight) & 0x01;
			uint32_t childZ = (tz >> childHeight) & 0x01;
			uint32_t childId = childZ << 2 | childY << 1 | childX;

			uint32_t childIndex = mNodeArray[nodeIndex].child(childId);

			if(childIndex < MaterialCount) return static_cast<MaterialId>(childIndex);

			// Prepare for next iteration. Can we replace 'nodeIndex' and 'childIndex' with a single variable?
			nodeIndex = childIndex;
			height--;
		}

		assert(false);
		return 0; // Should never get here!
	}

	////////////////////////////////////////////////////////////////////////////////
	// Private member functions
	////////////////////////////////////////////////////////////////////////////////

	// This function ensures that a child node has only a single parent (that it is not shared).
	// It does this by making a copy of the child and setting that as the child of the original parent.
	void Volume::splitChild(uint32_t nodeIndex, uint32_t childId)
	{
		// Get the index of the child node.
		uint32_t childIndex = mNodeArray[nodeIndex].child(childId);

		// We should avoid calling this if the child is not actually
		// shared. I don't think it will matter, but it is wasteful.
		assert(mNodeArray[childIndex].refCount() > 1);

		// Allocate a new node to be a copy of the child.
		uint32_t newChildIndex = mNodeArray.allocateNode();

		// Copy the old child's data to the new child. That is, iterate
		// over each of the child's children and copy them across.
		for(int i = 0; i < 8; i++)
		{
			mNodeArray.setChild(newChildIndex, i, mNodeArray[childIndex].child(i));
		}

		// Now replace the original child with the new one.
		mNodeArray.setChild(nodeIndex, childId, newChildIndex);

		// Same check we had at the start, but now checking
		// it is *not* shared rather than that it *is* shared.
		childIndex = mNodeArray[nodeIndex].child(childId);
		assert(mNodeArray[childIndex].refCount() == 1);
	}

	bool Volume::load(const std::string& filename)
	{
		std::ifstream file(filename, std::ios::binary);

		// FIXME - What should we do for error handling in this function? Return codes or exceptions?
		if(!file.is_open())
		{
			return false;
		}

		mNodeArray.read(file);
		return true;
	}

	void Volume::save(const std::string& filename)
	{
		// Should we compact here? This results in a merge and a defragmentation, both of which (and especially the
		// former) might be unexpected from a save operation. But they don't affect the observable behaviour of the
		// volume except that it will now use less memory and might be faster. I think the defrag is probably always
		// desirable, but certain internal operation in Cubiquity operate on an unmerged node array (e.g. the
		// voxelization), so I will need to ensure I don't save if I need this property.
		//
		// I might also later find that one of the operations is too slow to be performed as part of the save, but I'll
		// address that if I come to it.
		compact();

		std::ofstream file;
		file = std::ofstream(filename, std::ios::out | std::ios::binary);

		mNodeArray.write(file);

		file.close();
	}

	namespace Internals
	{
		/// This is an advanced function which should only be used if you
		/// understand the internal memory layout of Cubiquity's volume data.
		NodeArray& getNodeArray(Volume& volume)
		{
			return volume.mNodeArray;
		}

		const NodeArray& getNodeArray(const Volume& volume)
		{
			return volume.mNodeArray;
		}
	}
}
