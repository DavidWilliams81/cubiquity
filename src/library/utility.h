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

#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace Cubiquity
{
	// A little utility class useful for debugging and profiling.
	class Counter
	{
	public:
		Counter(const std::string& name) : mName(name) {}
		~Counter()
		{
			std::cout << mName << " = " << mValue << std::endl;
		}
		void inc() { mValue++; }
	private:
		std::string mName;
		uint64 mValue = 0;
	};

	class Timer
	{
	public:
		Timer(bool bAutoStart = true)
		{
			if (bAutoStart)
			{
				start();
			}
		}

		void start(void)
		{
			m_start = clock::now();
		}

		float elapsedTimeInSeconds(void)
		{
			std::chrono::duration<float> elapsed_seconds = clock::now() - m_start;
			return elapsed_seconds.count();
		}

		float elapsedTimeInMilliSeconds(void)
		{
			std::chrono::duration<float, std::milli> elapsed_milliseconds = clock::now() - m_start;
			return elapsed_milliseconds.count();
		}

		float elapsedTimeInMicroSeconds(void)
		{
			std::chrono::duration<float, std::micro> elapsed_microseconds = clock::now() - m_start;
			return elapsed_microseconds.count();
		}

	private:
		typedef std::chrono::system_clock clock;
		std::chrono::time_point<clock> m_start;
	};

	// FIXME - Can we make this take a const volume reference?
	template<typename Functor>
	void visitVolumeNodes(Cubiquity::Volume& volume, Functor&& callback)
	{
		Internals::NodeDAG& mDAG = Internals::getNodes(volume);
		const uint32_t rootNodeIndex = Internals::getRootNodeIndex(volume);

		const uint32 rootHeight = 32; // FIXME - Make this a constant somewhere?
		const Box3i rootBounds = Box3i::max();

		// Call the handler on the root.
		const bool processChildren = callback(mDAG, rootNodeIndex, rootBounds);

		// Process the root's children if requested and possible.
		const bool hasChildren = !Internals::isMaterialNode(rootNodeIndex);
		if (hasChildren && processChildren)
		{
			visitChildNodes(mDAG, rootNodeIndex, rootBounds, rootHeight, callback);
		}
	}

	// Call the callback on each child of the specified node.
	template<typename Functor>
	void visitChildNodes(Internals::NodeDAG& mDAG, uint32_t nodeIndex, const Box3i& bounds, uint32 height, Functor&& callback)
	{
		// Determine which bit may need to be flipped
		// to derive child bounds from parent bounds.
		const uint32 childHeight = height - 1;
		const uint32 bitToFlip = 0x01 << childHeight;

		for(uint32 childId = 0; childId < 8; childId++)
		{
			const uint32 childNodeIndex = mDAG[nodeIndex][childId];

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
			const bool processChildren = callback(mDAG, childNodeIndex, childBounds);

			// Process this child's children if requested and possible.
			const bool hasChildren = !Internals::isMaterialNode(childNodeIndex);
			if (hasChildren && processChildren)
			{
				visitChildNodes(mDAG, childNodeIndex, childBounds, childHeight, callback);
			}
		}
	}
	Cubiquity::Box3i computeBounds(Cubiquity::Volume& volume, MaterialId externalMaterial);
	std::pair<uint16_t, Cubiquity::Box3i> estimateBounds(Cubiquity::Volume& volume);

	struct HistogramEntry
	{
		bool overflow = false;
		uint64 count = 0;
	};

	typedef std::map<MaterialId, HistogramEntry> Histogram;
	Histogram computeHistogram(Volume& volume);
	void printHistogram(const Histogram& histogram);
}

#endif // CUBIQUITY_ALGORITHMS_H
