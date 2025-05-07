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
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Cubiquity
{
	namespace Internals
	{
		constexpr MaterialId MinMaterial = std::numeric_limits< MaterialId>::min();
		constexpr MaterialId MaxMaterial = std::numeric_limits< MaterialId>::max();
		constexpr u32     MaterialCount = static_cast<u32>(MaxMaterial) + 1;
		constexpr u64     VolumeSideLength = UINT64_C(1) << 32;

		bool isMaterialNode(u32 nodeIndex);

		typedef std::array<u32, 8> Node;

		class NodeStore
		{
		public:
			NodeStore() { mData = new Node[size()]; }
			~NodeStore() { delete[] mData; }

			const Node& operator[](u32 index) const { return mData[index]; }

			void setNode(u32 index, const Internals::Node& node);
			void setNodeChild(u32 nodeIndex, u32 childId, u32 newChildIndex);
			Node* data() const { return mData; }
			u32 size() const { return 0x3FFFFFF; }
			
		private:
			Node* mData = nullptr;
		};

		class NodeDAG
		{
		public:
			NodeDAG();

			const Node& operator[](u32 index) const { return mNodes[index]; }

			u32 bakedNodesBegin() const { return MaterialCount; }
			u32 bakedNodesEnd() const { return mBakedNodesEnd; }
			u32 editNodesBegin() const { return mEditNodesBegin; }
			u32 editNodesEnd() const { return mNodes.size(); }

			bool isBakedNode(u32 index) const { return index >= bakedNodesBegin() && index < bakedNodesEnd(); }
			bool isEditNode(u32 index) const { return index >= editNodesBegin() && index < editNodesEnd(); }

			NodeStore& nodes() { return mNodes; }
			const NodeStore& nodes() const { return mNodes; }

			u32 countNodes(u32 startNodeIndex) const;
			void countNodes(u32 startNodeIndex, std::unordered_set<u32>& usedIndices) const;

			void read(std::ifstream& file);
			void write(std::ofstream& file);


			bool isPrunable(const Node& node) const;

			u32 insert(const Node& node);
			u32 updateNodeChild(u32 nodeIndex, u32 childId, u32 newChildNodeIndex, bool forceCopy);

			void merge(u32 index);
			u32 mergeNode(u32 nodeIndex, std::unordered_map<Internals::Node, u32>& map, u32& nextSpace);

		private:
			NodeStore mNodes;
			u32 mBakedNodesEnd = MaterialCount;
			u32 mEditNodesBegin = 0;
		};
	}

	class Brush
	{
	public:
		virtual bool contains(const vec3f& point) const = 0;
		virtual Box3f bounds() const = 0;

		vec3f mCentre;
	};

	class SphereBrush : public Brush
	{
	public:
		SphereBrush(float x, float y, float z, float radius)
			:SphereBrush(vec3f(x, y, z), radius) {}
		SphereBrush(const vec3f& centre, float radius)
			:mRadiusSquared(radius * radius)
		{
			mCentre = centre;
			vec3f radiusAsVec = vec3f(radius);
			mBounds = Box3f(centre - radiusAsVec, centre + radiusAsVec);
		}

		bool contains(const vec3f& point) const
		{
			// WARNING - Dubious precision here - these values can be huge!
			float distX = point.x - mCentre.x;
			float distY = point.y - mCentre.y;
			float distZ = point.z - mCentre.z;

			i64 distSq = (distX * distX + distY * distY + distZ * distZ);

			return distSq < mRadiusSquared;
		}

		Box3f bounds() const
		{
			return mBounds;
		}

	public:
		
		float mRadiusSquared;
		Box3f mBounds;
	};

	class Volume;
	namespace Internals
	{
		// The Volume class provides a high-level interface to the voxel data, but for some purposes (e.g. rendering) it's
		// faster to work on the raw data. This should be done with care as it requires understanding the internal format.
		// Such code can access the internals only through this 'Internals' namespace, which serves as a warning that they
		// probably shouldn't be doing it.

		/// Provides access to the raw node data.
		NodeDAG& getNodes(Volume& volume);
		const NodeDAG& getNodes(const Volume& volume);
		//u32& getRootNodeIndex(Volume& volume);
		const u32 getRootNodeIndex(const Volume& volume);
	}

	class Volume
	{
	public:

		Volume();
		Volume(const std::string& filename);

		void fill(MaterialId matId);

		u32 rootNodeIndex() const;
		void setRootNodeIndex(u32 newRootNodeIndex);

		void setTrackEdits(bool trackEdits);
		bool undo();
		bool redo();

		void setVoxelRecursive(i32 x, i32 y, i32 z, MaterialId matId);
		u32 setVoxelRecursive(u32 ux, u32 uy, u32 uz, MaterialId matId, u32 nodeIndex, int nodeHeight);

		template <typename ArrayType>
		void setVoxel(const ArrayType& position, MaterialId matId);
		void setVoxel(i32 x, i32 y, i32 z, MaterialId matId);

		void fillBrush(const Brush& brush, MaterialId matId);
		u32 fillBrush(const Brush& brush, MaterialId matId, u32 nodeIndex, int nodeHeight, i32 nodeLowerX, i32 nodeLowerY, i32 nodeLowerZ);

		void addVolume(const Volume& rhsVolume);
		u32 addVolume(const Volume& rhsVolume, u32 rhsNodeIndex, u32 nodeIndex, int nodeHeight, i32 nodeLowerX, i32 nodeLowerY, i32 nodeLowerZ);

		template <typename ArrayType>
		MaterialId voxel(const ArrayType& position) const;
		MaterialId voxel(i32 x, i32 y, i32 z) const;

		void bake();

		u32 countNodes() const { return mDAG.countNodes(rootNodeIndex()); };

		bool load(const std::string& filename);
		void save(const std::string& filename);

	private:

		friend Internals::NodeDAG& Internals::getNodes(Volume& volume);
		friend const Internals::NodeDAG& Internals::getNodes(const Volume& volume);
		//friend u32& Internals::getRootNodeIndex(Volume& volume);
		friend const u32 Internals::getRootNodeIndex(const Volume& volume);
		
		Internals::NodeDAG mDAG;

		// Note: The undo system is linear. If we apply operation A, undo it, and then apply operation B
		// then A is lost. But it we wanted to we could still track it, because the appropriate edit still
		// exists (it is jst unreferenced). But having such an undo history tree will probably be confusing
		// for the user, moving forwards and backwards is probably enough.
		bool mTrackEdits = false;
		std::vector<u32> mRootNodeIndices;
		u32 mCurrentRoot = 0;
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
	template<>
	struct hash<Cubiquity::Internals::Node>
	{
		Cubiquity::u32 seed = 0;

		std::size_t operator()(const Cubiquity::Internals::Node& node) const noexcept
		{
			return Cubiquity::Internals::murmurHash3(&(node[0]), sizeof(uint32_t) * 8, seed);
		}
	};
}

#endif //CUBIQUITY_VOLUME_H
