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

	Box3i childBounds(Box3i nodeBounds, uint childId);

	template<typename Functor>
	void traverseNodesRecursive(Cubiquity::Volume& volume, Functor&& callback)
	{
		Internals::NodeDAG& mDAG = Internals::getNodes(volume);
		uint32_t rootNodeIndex = Internals::getRootNodeIndex(volume);

		Vector3i rootLowerBound = Vector3i::filled(std::numeric_limits<int32>::min());
		Vector3i rootUpperBound = Vector3i::filled(std::numeric_limits<int32>::max());

		traverseNodesRecursive(mDAG, rootNodeIndex, Box3i(rootLowerBound, rootUpperBound), callback);
	}

	template<typename Functor>
	void traverseNodesRecursive(Internals::NodeDAG& mDAG, uint32_t nodeIndex, Box3i nodeBounds, Functor&& callback)
	{
		bool processChildren = callback(mDAG, nodeIndex, nodeBounds);

		if (processChildren && (!Internals::isMaterialNode(nodeIndex)))
		{
			for(uint32 childId = 0; childId < 8; childId++)
			{
				const uint32 childNodeIndex = mDAG[nodeIndex][childId];
				const Box3i childNodeBounds = childBounds(nodeBounds, childId);
				traverseNodesRecursive(mDAG, childNodeIndex, childNodeBounds, callback);
			}
		}
	}

	typedef std::map<MaterialId, uint64> Histogram;
	Cubiquity::Box3i computeBounds(Cubiquity::Volume& volume, MaterialId externalMaterial);
	std::pair<uint16_t, Cubiquity::Box3i> estimateBounds(Cubiquity::Volume& volume);
	Histogram computeHistogram(Volume& volume, const Box3i& bounds);
	void printHistogram(const Histogram& histogram);

	// Implementation of a Linear Feedback Shift Register
	//
	// A list of masks for maximul-length LSFRs at various bit sizes is available here:
	//
	//     https://users.ece.cmu.edu/~koopman/lfsr/index.html
	//
	// Note that the start state need to be one of those visited otherwise we'll never get back to it.
	// For a maximul-length LSFR this might mean the bit size of the start state has to be no bigger
	// than the bit size of the mask? I didn't check this...
	class GaloisLFSR
	{
	public:
		GaloisLFSR(uint32_t mask, uint32_t startState = 1);

		void next();
		uint32_t state();

		uint64_t computePeriod();

	private:
		uint32_t mMask;
		uint32_t mState;
	};

	class ShuffledSequence
	{
	public:
		ShuffledSequence(uint32 sequenceLength);

		void next();
		uint32_t state();

		static uint32_t maximulLengthMask(int sizeInBits);

	private:
		uint32 mSequenceLength;
		GaloisLFSR mGaloisLFSR;

		static const uint32_t MaximulLengthMasks[33];
	};
}

#endif // CUBIQUITY_ALGORITHMS_H
