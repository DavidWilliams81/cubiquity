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
#include <map>
#include <memory>
#include <string>
#include <iostream>

namespace Cubiquity
{
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

	struct Colour
	{
	public:
		Colour() {}
		Colour(uint8 l, uint8 a = 0xff) : red(l), green(l), blue(l), alpha(a) {}
		Colour(uint8 r, uint8 g, uint8 b, uint8 a = 0xff) : red(r), green(g), blue(b), alpha(a) {}

		uint8 red, green, blue, alpha;
	};

	class Image
	{
	public:
		Image(uint width, uint height) : mWidth(width), mHeight(height)
		{
			//mData = std::make_unique<Colour[]>(width*height);
			mData = std::unique_ptr<Colour[]>{ new Colour[width*height] };
		}

		const Colour& pixel(uint x, uint y) { return mData[y * mWidth + x]; }
		void setPixel(uint x, uint y, const Colour& col) { mData[y * mWidth + x] = col; }

		uint32 hash();
		void save(const std::string& filename);
	private:
		uint mWidth;
		uint mHeight;
		std::unique_ptr<Colour[]> mData;
	};

	// When implementing and/or using a ProgressMonitor, keep in mind that a task with x steps
	// will have x+1 states. One for no steps completed, and then one for each of the steps.
	class ProgressMonitor
	{
	public:
		virtual void startTask(const std::string& name) = 0;
		virtual void setProgress(int minimum, int current, int maximum) = 0;
		virtual void finishTask() = 0;
	};

	class TerminalProgressBar : public ProgressMonitor
	{
		const uint NameLength = 30;
		const uint BarLength = 60;

	public:
		void startTask(const std::string& name)
		{
			mTaskName = name;

			// Fixed length simplifies layout
			if (mTaskName.size() > NameLength) { mTaskName.replace(NameLength - 3, 3, "..."); }
			mTaskName.resize(NameLength);

			mTimer.start();
		}

		void setProgress(int minimum, int current, int maximum)
		{
			assert(minimum < maximum && current >= minimum && current <= maximum);

			// Work out progress
			float progress = static_cast<float>(current - minimum) /
				static_cast<float>(maximum - minimum);
			long markerCount = lroundf(progress * BarLength);

			// Draw progress bar ().			
			std::cout << mTaskName << "[";
			for (long i = 0; i < BarLength; i++) { std::cout << (i < markerCount ? '=' : ' '); }
			std::cout << "] ";

			// Write the time with fixed precision to make sure it overwrites the previous value.
			std::streamsize oldPrec = std::cout.precision(3);
			std::ios_base::fmtflags oldFlags = std::cout.setf(std::ios::fixed);
			std::cout << mTimer.elapsedTimeInSeconds() << "s";
			std::cout.precision(oldPrec);
			std::cout.flags(oldFlags);

			std::cout << "\r" << std::flush; // '\r' without '\n' goes back to start of the line
		}

		void finishTask()
		{
			setProgress(0, 1, 1); // Max progress
			std::cout << std::endl;
		}

	private:
		std::string mTaskName;
		Timer mTimer;
	};

	// FIXME - Thse volumes should be const referenes.
	// WARNING - I think the use of startPoint has not actually been tested.
	template<typename Functor>
	void traverseNodes(Cubiquity::Volume& volume, Functor&& callback, Vector3d startPoint = Vector3d(0.0))
	{
		Internals::NodeArray& mNodeArray = Internals::getNodeArray(volume);

		struct NodeState
		{
			NodeState()
				: mIndex(0)
				, mProcessedChildCount(0)
				, mLowerCorner(std::numeric_limits<int32_t>::min())
				, mCentre(0.0f, 0.0f, 0.0f)
				, mNearestChild(0) {}

			uint32_t mIndex;
			uint32_t mProcessedChildCount;
			Vector3i mLowerCorner;
			Vector3f mCentre;
			uint8_t mNearestChild; // Index (0 - 7) of nearest child.
		};

		Box3i bounds;

		const int rootHeight = Internals::logBase2(VolumeSideLength);
		// Note that the first two elements of this stack never actually get used.
		// Leaf and almost-leaf nodes(heights 0 and 1) never get put on the stack.
		// We accept this wasted space, rather than subtracting two on every access.
		std::vector<NodeState> nodeStateStack(rootHeight + 1);

		int nodeHeight = rootHeight;
		NodeState& rootNodeState = nodeStateStack[rootHeight];
		rootNodeState.mIndex = Internals::RootNodeIndex; // Root node is always at index 0.
		rootNodeState.mCentre = Vector3f(-0.5f);

		Vector3i rootLowerBound(std::numeric_limits<int32>::min());
		Vector3i rootUpperBound(std::numeric_limits<int32>::max());

		if (startPoint.x() > rootNodeState.mCentre.x()) rootNodeState.mNearestChild |= 0x01;
		if (startPoint.y() > rootNodeState.mCentre.y()) rootNodeState.mNearestChild |= 0x02;
		if (startPoint.z() > rootNodeState.mCentre.z()) rootNodeState.mNearestChild |= 0x04;

		callback(mNodeArray, Internals::RootNodeIndex, Box3i(rootLowerBound, rootUpperBound));

		while (true)
		{
			NodeState& nodeState = nodeStateStack[nodeHeight];

			if (nodeState.mProcessedChildCount < 8)
			{
				// Based on Octree traversal method here: https://www.flipcode.com/archives/Harmless_Algorithms-Issue_02_Scene_Traversal_Algorithms.shtml#octh
				const uint8_t bitToggles[] = { 0x00, 0x01, 0x02, 0x04, 0x03, 0x05, 0x06, 0x07 };
				uint32_t childId = nodeState.mNearestChild ^ bitToggles[nodeState.mProcessedChildCount];

				nodeState.mProcessedChildCount++;

				Internals::Node node = mNodeArray[nodeState.mIndex];

				NodeState& childNodeState = nodeStateStack[nodeHeight - 1];

				uint32_t childHeight = nodeHeight - 1;
				uint32_t childSideLength = 1 << (childHeight);
				uint32_t childX = (childId >> 0) & 0x01;
				uint32_t childY = (childId >> 1) & 0x01;
				uint32_t childZ = (childId >> 2) & 0x01;

				// Compute the centre of the child based on the parent and child index. Could stick the
				// '(float(childIndex >> 0 & 0x01) - 0.5f)' part in an eight-element LUT if we find that it is faster?
				childNodeState.mCentre.data[0] = nodeState.mCentre.x() + (float(childId >> 0 & 0x01) - 0.5f) * childSideLength;
				childNodeState.mCentre.data[1] = nodeState.mCentre.y() + (float(childId >> 1 & 0x01) - 0.5f) * childSideLength;
				childNodeState.mCentre.data[2] = nodeState.mCentre.z() + (float(childId >> 2 & 0x01) - 0.5f) * childSideLength;

				// childX/Y/Z are all zero or one.
				childNodeState.mLowerCorner[0] = nodeState.mLowerCorner[0] + (childSideLength * childX);
				childNodeState.mLowerCorner[1] = nodeState.mLowerCorner[1] + (childSideLength * childY);
				childNodeState.mLowerCorner[2] = nodeState.mLowerCorner[2] + (childSideLength * childZ);

				if (startPoint.x() > childNodeState.mCentre.x()) childNodeState.mNearestChild |= 0x01;
				if (startPoint.y() > childNodeState.mCentre.y()) childNodeState.mNearestChild |= 0x02;
				if (startPoint.z() > childNodeState.mCentre.z()) childNodeState.mNearestChild |= 0x04;

				childNodeState.mIndex = node.child(childId);

				Vector3i childLowerBound = childNodeState.mLowerCorner;
				Vector3i childUpperBound = childNodeState.mLowerCorner + Vector3i(childSideLength - 1, childSideLength - 1, childSideLength - 1);

				callback(mNodeArray, childNodeState.mIndex, Box3i(childLowerBound, childUpperBound));

				// Descend to any further children
				if (childNodeState.mIndex >= Internals::MaterialCount)
				{
					childNodeState.mProcessedChildCount = 0;

					nodeHeight -= 1;
				}
			}
			else
			{
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

	typedef std::map<MaterialId, uint64> Histogram;
	uint64_t countNodes(Cubiquity::Volume& volume);
	Cubiquity::Box3i computeBounds(Cubiquity::Volume& volume, bool(*include)(Cubiquity::MaterialId));
	std::pair<uint16_t, Cubiquity::Box3i> estimateBounds(Cubiquity::Volume& volume);
	Histogram computeHistogram(Volume& volume, const Box3i& bounds);
	void printHistogram(const Histogram& histogram);
	void saveVolumeAsImages(Cubiquity::Volume& volume, const std::string& filename, ProgressMonitor* progMon = nullptr);

	Geometry loadObjFile(const std::string& path, const std::string& filename);

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
