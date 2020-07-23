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
#ifndef CUBIQUITY_VOLUME_H
#define CUBIQUITY_VOLUME_H

#include "base.h"

#include <string>
#include <vector>

namespace Cubiquity
{
	static const MaterialId MinMaterial = 0;
	static const MaterialId MaxMaterial = 65535;
	static const uint64 VolumeSideLength = UINT64_C(1) << 32;

	namespace Internals
	{
		static const uint32_t EmptyNodeIndex = 0;
		static const uint32_t MaterialCount = static_cast<uint32_t>(MaxMaterial) + 1;
		static const uint32_t RootNodeIndex = MaterialCount;

		class Node
		{
		public:

			Node(uint32_t initialChildValue = 0);

			bool operator== (const Node& rhs) const;
			bool operator!= (const Node& rhs) const;

			uint32_t refCount() const;

			uint32_t child(uint32_t id) const;
			void setChild(uint32_t id, uint32_t value);

			bool allChildrenAre(uint32_t value) const;
			bool allChildrenAreLessThan(uint32_t value) const;

		public:
			friend class NodeArray;

			uint32_t mRefCount;
			uint32_t mChildren[8];
		};

		// custom hash can be a standalone function object:
		struct NodeHash
		{
			uint32_t operator()(Node const& node, uint32_t seed = 0) const noexcept
			{
				/*size_t hashResult = 0;

				for (int i = 0; i < 8; i++)
				{
				hashResult = hashResult ^ std::hash<uint32_t>{}(node.child(i));
				}

				return hashResult;*/

				return murmurHash3(&(node.mChildren[0]), sizeof(uint32_t) * 8, seed);
			}

			std::vector<Internals::Node>* mNodeData = nullptr;
		};
	}

	namespace Internals
	{
		class NodeArray
		{
		public:
			NodeArray(uint32_t size, MaterialId initialValue);
			NodeArray(const NodeArray& other);
			~NodeArray();

			NodeArray& operator=(const NodeArray& other);

			uint32_t allocateNode(uint32_t initialChildValue = 0);
			void resize(uint32_t size);
			uint32_t size() const;

			void setChild(uint32_t nodeIndex, uint32_t childId, uint32_t newChildIndex);

			// This function chenges the structure of the tree by merging nodes, but it does not
			// move the nodes around in memory (so the array may be fragmented after calling.
			void mergeOctree();
			void mergeOctreeBruteForce();

			// This function does not change the structure of the tree, but
			// it moves the nodes around in memory so that they are contiguous.
			// We should also check if it lays the nodes out in an optimal
			// order (whatever that might be).
			void defragment();

			const Internals::Node& operator[](uint32_t index) const;

			void read(std::ifstream& file);
			void write(std::ofstream& file);

		private:
			uint32_t mSize;
			Internals::Node* mData;

			uint32_t mLastAllocation;
		};
	}

	class Volume;
	namespace Internals
	{
		// The Volume class provides a high-level interface to the voxel data, but for some purposes (e.g. rendering) it's
		// faster to work on the raw data. This should be done with care as it requires understanding the internal format.
		// Such code can access the internals only through this 'Internals' namespace, which serves as a warning that they
		// probably shouldn't be doing it.

		/// Provides access to the raw node data.
		NodeArray& getNodeArray(Volume& volume);
		const NodeArray& getNodeArray(const Volume& volume);
	}

	class Volume
	{
	public:

		Volume(MaterialId initialValue = 0);
		Volume(const std::string& filename);

		template <typename ArrayType>
		void setVoxel(const ArrayType& position, MaterialId matId, bool checkIfAlreadySet = true);
		void setVoxel(int32_t x, int32_t y, int32_t z, MaterialId matId, bool checkIfAlreadySet = true);

		template <typename ArrayType>
		MaterialId voxel(const ArrayType& position);
		MaterialId voxel(int32_t x, int32_t y, int32_t z);

		void compact();

		bool load(const std::string& filename);
		void save(const std::string& filename);


	private:

		void splitChild(uint32_t nodeIndex, uint32_t childIndex);

		friend Internals::NodeArray& Internals::getNodeArray(Volume& volume);
		friend const Internals::NodeArray& Internals::getNodeArray(const Volume& volume);

		Internals::NodeArray mNodeArray;
	};

	// Implementation of templatised accessors
	template <typename ArrayType> MaterialId Volume::voxel(const ArrayType& position)
	{
		return voxel(position[0], position[1], position[2]);
	}

	template <typename ArrayType> void Volume::setVoxel(const ArrayType& position, MaterialId matId, bool checkIfAlreadySet)
	{
		setVoxel(position[0], position[1], position[2], matId, checkIfAlreadySet);
	}
}

#endif //CUBIQUITY_VOLUME_H
