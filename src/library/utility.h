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
#ifndef CUBIQUITY_ALGORITHMS_H
#define CUBIQUITY_ALGORITHMS_H

#include "base.h"
#include "geometry.h"
#include "storage.h"

#include <cinttypes>
#include <map>
#include <memory>

namespace Cubiquity
{
	// For debugging and profiling (usually instantiated as a global)
	class Counter
	{
	public:
		// 'name' is not copied, so it must outlive the Counter (e.g. a string literal).
		Counter(const char* name) : mName(name) {}
		~Counter() { Internals::log_debug("%s = %" PRIu64, mName, mValue); }
		void inc() { mValue++; }
	private:
		const char* mName;
		u64 mValue = 0;
	};

	// Simple timer class for profiling
	class Timer
	{
	public:
		Timer(bool bAutoStart = true)
		{
			if (bAutoStart)	{
				start();
			}
		}

		void start(void) { m_start = seconds_since_epoch(); };
		double elapsed_seconds(void) { return seconds_since_epoch() - m_start; };
		double elapsed_milliseconds(void) { return elapsed_seconds() * 1000.0; };

	private:
		double m_start = 0; // Plain double so we don't need <chrono>
		double seconds_since_epoch();
	};

	// FIXME - Can we make this take a const volume reference?
	template<typename Functor>
	void visitVolumeNodes(Cubiquity::Volume& volume, Functor&& callback)
	{
		Internals::NodeStore& mNodeStore = Internals::getNodes(volume);
		const u32 rootNodeIndex = Internals::getRootNodeIndex(volume);

		const u32 rootHeight = 32; // FIXME - Make this a constant somewhere?
		const Box3i rootBounds = Box3i::max();

		// Call the handler on the root.
		const bool processChildren = callback(mNodeStore, rootNodeIndex, rootBounds);

		// Process the root's children if requested and possible.
		const bool hasChildren = !Internals::isMaterialNode(rootNodeIndex);
		if (hasChildren && processChildren)
		{
			visitChildNodes(mNodeStore, rootNodeIndex, rootBounds, rootHeight, callback);
		}
	}

	// Call the callback on each child of the specified node.
	template<typename Functor>
	void visitChildNodes(Internals::NodeStore& mNodeStore, u32 nodeIndex, const Box3i& bounds, u32 height, Functor&& callback)
	{
		// Determine which bit may need to be flipped
		// to derive child bounds from parent bounds.
		const u32 childHeight = height - 1;
		const u32 bitToFlip = 0x01 << childHeight;

		for(u32 childId = 0; childId < 8; childId++)
		{
			const u32 childNodeIndex = mNodeStore[nodeIndex][childId];

			// Set the child bounds to be the same as the parent bounds and then collapse 
			// three of the faces (selected via the child id). Note that we could actually
			// skip the copy and modify the parent bounds in place as long as we flipped
			// the bits again after proessing. We could also iterate over the three child
			// components directly (rather than iterating over the child id and then
			// extracting the components), which also lets us avoid flipping all three
			// axes on every iteration. However, these changes complicate the code and
			// gave only a small speed improvement.
			Box3i childBounds = bounds;
			childBounds.mExtents[((~childId) >> 0) & 0x01][0] ^= bitToFlip;
			childBounds.mExtents[((~childId) >> 1) & 0x01][1] ^= bitToFlip;
			childBounds.mExtents[((~childId) >> 2) & 0x01][2] ^= bitToFlip;

			// Call the handler on this child.
			const bool processChildren = callback(mNodeStore, childNodeIndex, childBounds);

			// Process this child's children if requested and possible.
			const bool hasChildren = !Internals::isMaterialNode(childNodeIndex);
			if (hasChildren && processChildren)
			{
				visitChildNodes(mNodeStore, childNodeIndex, childBounds, childHeight, callback);
			}
		}
	}
	Cubiquity::Box3i computeBounds(Cubiquity::Volume& volume, MaterialId externalMaterial);
	std::pair<u16, Cubiquity::Box3i> estimateBounds(Cubiquity::Volume& volume);

	struct HistogramEntry
	{
		bool overflow = false;
		u64 count = 0;
	};

	typedef std::map<MaterialId, HistogramEntry> Histogram;
	Histogram computeHistogram(Volume& volume);
	void printHistogram(const Histogram& histogram);
}

#endif // CUBIQUITY_ALGORITHMS_H
