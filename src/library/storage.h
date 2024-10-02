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
		constexpr uint32     MaterialCount = static_cast<uint32>(MaxMaterial) + 1;
		constexpr uint64     VolumeSideLength = UINT64_C(1) << 32;

		bool isMaterialNode(uint32 nodeIndex);

		typedef std::array<uint32_t, 8> Node;

		class NodeStore
		{
		public:
			NodeStore() { mData = new Node[size()]; }
			~NodeStore() { delete[] mData; }

			const Node& operator[](uint32_t index) const { return mData[index]; }

			void setNode(uint32 index, const Internals::Node& node);
			void setNodeChild(uint32 nodeIndex, uint32 childId, uint32 newChildIndex);
			Node* data() const { return mData; }
			uint32 size() const { return 0x3FFFFFF; }
			
		private:
			Node* mData = nullptr;
		};

		class NodeDAG
		{
		public:
			NodeDAG();

			const Node& operator[](uint32_t index) const { return mNodes[index]; }

			uint32 bakedNodesBegin() const { return MaterialCount; }
			uint32 bakedNodesEnd() const { return mBakedNodesEnd; }
			uint32 editNodesBegin() const { return mEditNodesBegin; }
			uint32 editNodesEnd() const { return mNodes.size(); }

			bool isBakedNode(uint32 index) const { return index >= bakedNodesBegin() && index < bakedNodesEnd(); }
			bool isEditNode(uint32 index) const { return index >= editNodesBegin() && index < editNodesEnd(); }

			NodeStore& nodes() { return mNodes; }
			const NodeStore& nodes() const { return mNodes; }

			uint32 countNodes(uint32 startNodeIndex) const;
			void countNodes(uint32 startNodeIndex, std::unordered_set<uint32>& usedIndices) const;

			void read(std::ifstream& file);
			void write(std::ofstream& file);


			bool isPrunable(const Node& node) const;

			uint32 insert(const Node& node);
			uint32 updateNodeChild(uint32 nodeIndex, uint32 childId, uint32 newChildNodeIndex, bool forceCopy);

			void merge(uint32 index);
			uint32 mergeNode(uint32 nodeIndex, std::unordered_map<Internals::Node, uint32>& map, uint32& nextSpace);

		private:
			NodeStore mNodes;
			uint32 mBakedNodesEnd = MaterialCount;
			uint32 mEditNodesBegin = 0;
		};
	}

	class Brush
	{
	public:
		virtual bool contains(const Vector3f& point) const = 0;
		virtual Box3f bounds() const = 0;

		Vector3f mCentre;
	};

	class SphereBrush : public Brush
	{
	public:
		SphereBrush(const Vector3f& centre, float radius)
			:mRadiusSquared(radius * radius)
		{
			mCentre = centre;
			Vector3f radiusAsVec = Vector3f::filled(radius);
			mBounds = Box3f(centre - radiusAsVec, centre + radiusAsVec);
		}

		bool contains(const Vector3f& point) const
		{
			// WARNING - Dubious precision here - these values can be huge!
			float distX = point.x() - mCentre.x();
			float distY = point.y() - mCentre.y();
			float distZ = point.z() - mCentre.z();

			int64 distSq = (distX * distX + distY * distY + distZ * distZ);

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
		//uint32& getRootNodeIndex(Volume& volume);
		const uint32 getRootNodeIndex(const Volume& volume);
	}

	class Volume
	{
	public:

		Volume();
		Volume(const std::string& filename);

		void fill(MaterialId matId);

		uint32 rootNodeIndex() const;
		void setRootNodeIndex(uint32 newRootNodeIndex);

		void setTrackEdits(bool trackEdits);
		bool undo();
		bool redo();

		void setVoxelRecursive(int32_t x, int32_t y, int32_t z, MaterialId matId);
		uint32 setVoxelRecursive(uint32 ux, uint32 uy, uint32 uz, MaterialId matId, uint32 nodeIndex, int nodeHeight);

		template <typename ArrayType>
		void setVoxel(const ArrayType& position, MaterialId matId);
		void setVoxel(int32_t x, int32_t y, int32_t z, MaterialId matId);

		void fillBrush(const Brush& brush, MaterialId matId);
		uint32 fillBrush(const Brush& brush, MaterialId matId, uint32 nodeIndex, int nodeHeight, int32 nodeLowerX, int32 nodeLowerY, int32 nodeLowerZ);

		void addVolume(const Volume& rhsVolume);
		uint32 addVolume(const Volume& rhsVolume, uint32 rhsNodeIndex, uint32 nodeIndex, int nodeHeight, int32 nodeLowerX, int32 nodeLowerY, int32 nodeLowerZ);

		template <typename ArrayType>
		MaterialId voxel(const ArrayType& position) const;
		MaterialId voxel(int32_t x, int32_t y, int32_t z) const;

		void bake();

		uint32 countNodes() const { return mDAG.countNodes(rootNodeIndex()); };

		bool load(const std::string& filename);
		void save(const std::string& filename);

	private:

		friend Internals::NodeDAG& Internals::getNodes(Volume& volume);
		friend const Internals::NodeDAG& Internals::getNodes(const Volume& volume);
		//friend uint32& Internals::getRootNodeIndex(Volume& volume);
		friend const uint32 Internals::getRootNodeIndex(const Volume& volume);
		
		Internals::NodeDAG mDAG;

		// Note: The undo system is linear. If we apply operation A, undo it, and then apply operation B
		// then A is lost. But it we wanted to we could still track it, because the appropriate edit still
		// exists (it is jst unreferenced). But having such an undo history tree will probably be confusing
		// for the user, moving forwards and backwards is probably enough.
		bool mTrackEdits = false;
		std::vector<uint32> mRootNodeIndices;
		uint32 mCurrentRoot = 0;
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
		Cubiquity::uint32 seed = 0;

		std::size_t operator()(const Cubiquity::Internals::Node& node) const noexcept
		{
			return Cubiquity::Internals::murmurHash3(&(node[0]), sizeof(uint32_t) * 8, seed);
		}
	};
}

#endif //CUBIQUITY_VOLUME_H
