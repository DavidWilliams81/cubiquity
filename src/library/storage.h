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
#include "geometry.h"

#include <array>
#include <climits>
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Cubiquity
{
	//! An integer type used to represent a material.
	/*!
		An environment in Cubiquity is stored as a 3D grid (see Volume) of material identifiers, with each one specifying the material
		at a given point in space. That is, MaterialId is the type used for each voxel. The material identifiers are simple integer
		values which the user is expected to map so a more complete material description via some application-specific method. For example,
		an application might use [1='rock', 2='wood', 3='glass', ...], or might map striaght to colours [1='red', 2='yellow', 3='blue', ...].

		The MaterialId is a 16-bit value but it is *not* expected that you use the entire 16-bit range. Most applications will only have a
		small number of materials, perhaps up to a few hundred. The compression method used in Cubiquity works by identifying identical
		regions in the scene and storing them only once, and using more materials will generally result in larger file sizes and/or higher
		memory usage.

		Given the above, it is reasonable to ask why an 8-bit data type is not used instead? There are a few potential resons why a 16 bit type is still useful:

		* 256 is too few materials
		* Encoding colours or other sparse values
		* Perhaps spatially large volumes aren't needed.

		Note typedef should not be changed to differnt integer type.
	*/
	typedef u8 MaterialId;

	namespace Internals
	{
		constexpr MaterialId MinMaterial = std::numeric_limits< MaterialId>::min();
		constexpr MaterialId MaxMaterial = std::numeric_limits< MaterialId>::max();
		constexpr u32      MaterialCount = static_cast<u32>(MaxMaterial) + 1;
		constexpr u64   VolumeSideLength = UINT64_C(1) << 32; // FIXME - Doesn't belong here?

		 // Reserving a node index costs almost nothing and it sometimes useful
		const u32 InvalidNodeIndex = 0xffffffff;

		// A node is just an array of child indices
		typedef std::array<u32, 8> Node;
		inline Node make_node(u32 i) { return Node{ i, i, i, i, i, i, i, i }; }
		inline bool isMaterialNode(u32 nodeIndex) { return nodeIndex < MaterialCount; }
		bool isMaterialNode(const Node& node);

		// Container for our nodes
		typedef std::vector<Node> NodeVector;
		NodeVector make_node_vector(u32 count = MaterialCount);

		class NodeStore
		{

		public:
			NodeStore();

			// Non-copyable (deleted)
			NodeStore(const NodeStore&) = delete;
			NodeStore& operator=(const NodeStore&) = delete;

			// Element access
			const Node& operator[](u32 index) const { return mNodeVector[index]; }
			u32 getNodeChild(u32 nodeIndex, u32 childId) const;
			u32 setNodeChild(u32 nodeIndex, u32 childId, u32 newChildNodeIndex);

			u32 sharedNodesBegin() const { return MaterialCount; }
			u32 sharedNodesEnd() const { return mSharedNodesEnd; }
			u32 unsharedNodesBegin() const { return mSharedNodesEnd; } // As above
			u32 unsharedNodesEnd() const { return mNodeVector.size(); }			

			u32 countNodes(u32 startNodeIndex) const;
			void countNodes(u32 startNodeIndex, std::unordered_set<u32>& usedIndices) const;

			void read(std::ifstream& file);
			void write(std::ofstream& file);

			u32 merge(u32 index);
			u32 merge_node(u32 nodeIndex, NodeVector& merged_node_vector);

			const void* rawBytesPtr() const { return mNodeVector.data(); }
			u32 rawBytesCount() const { return mNodeVector.size() * sizeof(mNodeVector[0]); }

			u32 cloneRoot(u32 rootIndex);

			void print(std::ostream& os);

		private:
			u32 append(const Node& node);

			bool isShared(u32 index) const { return index < sharedNodesEnd(); } // Includes material nodes

			NodeVector mNodeVector;
			u32 mSharedNodesEnd; // Also begining of unshared nodes
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
		NodeStore& getNodes(Volume& volume);
		const NodeStore& getNodes(const Volume& volume);
		//u32& getRootNodeIndex(Volume& volume);
		const u32 getRootNodeIndex(const Volume& volume);
	}

	typedef std::function<MaterialId(MaterialId, MaterialId)> MaterialCombiner;

	class Volume
	{
	public:

		Volume();
		Volume(const std::string& filename);

		// Non-copyable (deleted)
		Volume(const Volume&) = delete;
		Volume& operator=(const Volume&) = delete;

		// FIXME - Default move constructor doesn't work, I think because there
		// is no move constructor on NodeStore (which manages a raw pointer) and
		// so it falls back on the copy constructor. 
		// But moveable (defaulted)
		//Volume(Volume&& vol) = default;
		//Volume& operator=(Volume&& vol) = default;

		void fill(MaterialId matId);

		u32 rootNodeIndex() const;
		void setRootNodeIndex(u32 newRootNodeIndex);

		void checkpoint();
		void undo();
		void redo();

		void setVoxel(i32 x, i32 y, i32 z, MaterialId matId);
		u32 setVoxel(u32 ux, u32 uy, u32 uz, MaterialId matId, u32 nodeIndex, int nodeHeight);

		template <typename ArrayType>
		void setVoxel(const ArrayType& position, MaterialId matId);

		void combine(const Volume& rhsVolume, const MaterialCombiner& combiner);
		u32 combine(const Volume& rhsVolume, u32 rhsNodeIndex, u32 nodeIndex, const MaterialCombiner& combiner);

		template <typename ArrayType>
		MaterialId voxel(const ArrayType& position) const;
		MaterialId voxel(i32 x, i32 y, i32 z) const;

		void bake();

		u32 countNodes() const { return mNodeStore.countNodes(rootNodeIndex()); };

		bool load(const std::string& filename);
		void save(const std::string& filename);

		void print(std::ostream& os);

	private:

		friend Internals::NodeStore& Internals::getNodes(Volume& volume);
		friend const Internals::NodeStore& Internals::getNodes(const Volume& volume);
		//friend u32& Internals::getRootNodeIndex(Volume& volume);
		friend const u32 Internals::getRootNodeIndex(const Volume& volume);
		
	public:
		Internals::NodeStore mNodeStore;

	private:
		std::vector<u32> mRootNodeIndices;
		int mCurrentRoot = 0;
	};

	// Implementation of templatised accessors
	template <typename ArrayType> MaterialId Volume::voxel(const ArrayType& position) const
	{
		return voxel(position[0], position[1], position[2]);
	}

	template <typename ArrayType> void Volume::setVoxel(const ArrayType& position, MaterialId matId)
	{
		setVoxel(position[0], position[1], position[2], matId);
	}
}

namespace std
{
	static Cubiquity::u64 fnv1a_64(const Cubiquity::u8* data, size_t size)
	{
		Cubiquity::u64 hash = UINT64_C(0xcbf29ce484222325);
		for (size_t i = 0; i < size; i++) {
			hash ^= data[i];
			hash *= UINT64_C(0x00000100000001B3);
		}
		return hash;
	}

	template<>
	struct hash<Cubiquity::Internals::Node>;
}

#endif //CUBIQUITY_VOLUME_H
