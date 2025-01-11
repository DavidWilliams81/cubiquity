#include "test_volume.h"

#include "position_enumerator.h"

#include "base/logging.h"

#include "cubiquity.h"
#include "utility.h"
#include "storage.h"

#include "fractal_noise.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <thread>
#include <mutex>
#include <set>

using namespace Cubiquity;
using namespace Cubiquity::Internals;

bool checkIntegrity(Volume& volume)
{
	const NodeStore& nodes = getNodes(volume).nodes();

	// Reserved nodes always contain dummy data.
	//for (uint32_t i = 0; i < MaterialCount; i++)
	{
		/*if (!nodes[i].refCount() == 0xFFFFFFFF)
			return false;*/
		/*if (!nodes[i].allChildrenAre(0xFFFFFFFF))
			return false;*/
	}

	// The root node should always have a ref count of zero.
	// The children can be any mix of zeros, ones, and actual indices.
	/*if (nodes[RootNodeIndex].refCount() != 0)
		return false;*/

	// Now check the rest of the nodes.
	//for (uint32_t nodeIndex = nodes.hashMapBegin() + 1; nodeIndex < nodes.hashMapEnd(); nodeIndex++)
	{
		// If a node is not referenced then it should not have any children. I'm not
		// certain that we need to enforce this but it might get complex if we don't.

		//FIXME - Used to use ref count, wht to do here?
		/*if (nodes[nodeIndex].refCount() == 0)
		{
			if (!nodes[nodeIndex].allChildrenAreLessThan(MaterialCount))
				return false;
		}*/
		// Otherwise if it is referenced then it should not have all children
		// pointing to the same material (in that case the referring
		// node should have been pointed directly at the relevent child).
		// Howeer, note that a node can have eight identical children if they
		// are pointing at a real (non-material) node. For example, a 3D
		// checkerboard will be stored as a single chain of such nodes.
		/*else
		{
			for (uint32_t material = 0; material < MaterialCount; material++)
			{
				if (nodes[nodeIndex].allChildrenAre(material))
					return false;
			}
		}*/
	}

	// Also check that no node points at itself.
	/*for (uint32_t nodeIndex = nodes.hashMapBegin(); nodeIndex < nodes.hashMapEnd(); nodeIndex++)
	{
		for (uint32_t childId = 0; childId < 8; childId++)
		{
			uint32_t childIndex = nodes[nodeIndex][childId];
			if (childIndex == nodeIndex)
				return false;
		}
	}*/

	// Fixme - should also check all nodes are pruned.

	// Lastly we re-compute ref counts for all nodes
	// and make sure they match the existing values.

	//FIXME - Used to use ref count, wht to do here?
	/*bool allRefCountsCorrect = true;
	uint32_t* correctRefCounts = new uint32_t[nodes.size()];
	memset(correctRefCounts, 0, sizeof(uint32_t) * nodes.size());

	// Compute ref counts
	for (uint32_t nodeIndex = RootNodeIndex; nodeIndex < nodes.size(); nodeIndex++)
	{
		for(uint32_t childId = 0; childId < 8; childId++)
		{
			uint32_t childIndex = nodes[nodeIndex].child(childId);
			correctRefCounts[childIndex]++;
		}
	}

	// Check ref counts
	for(uint32_t nodeIndex = RootNodeIndex; nodeIndex < nodes.size(); nodeIndex++)
	{
		if (nodes[nodeIndex].mRefCount != correctRefCounts[nodeIndex])
		{
			allRefCountsCorrect = false;
			break;
		}
	}

	delete[] correctRefCounts;
	return allRefCountsCorrect;*/

	return true;
}

std::set< std::pair<uint32_t, uint32_t> > mergeOpportunities(Volume& volume)
{
	NodeDAG& nodes = getNodes(volume);

	std::set< std::pair<uint32_t, uint32_t> > result;

	// Early out until we fix use of ref count.
	return result;

	/*std::vector< std::vector<uint32_t> > hashTable;

	hashTable.resize(100000);


	// Skip empty and full nodes, as they contain identical (but dummy) data giving a false merge opportunity.
	for (uint32_t index = nodes.mSharedEdits.mBegin; index <= nodes.mSharedEdits.mEnd; index++)
	{
		const Node& node = nodes.mNodes[index];

		//FIXME - Used to use ref count, wht to do here?
		//if (node.mRefCount == 0) continue;

		uint32_t hash = std::hash<Node>{}(node);
		hashTable[hash % hashTable.size()].push_back(index);
	}

	// Now find the matches in each bucket
	for (uint32_t bucketIndex = 0; bucketIndex < hashTable.size(); bucketIndex++)
	{
		const std::vector<uint32_t>& bucket = hashTable[bucketIndex];

		for (uint32_t outerIndex = 0; outerIndex < bucket.size(); outerIndex++)
		{
			uint32_t outerNodeIndex = bucket[outerIndex];
			const Node& outerNode = nodes.mNodes[outerNodeIndex];
			//assert(outerNode.mRefCount > 0);

			for (uint32_t innerIndex = outerIndex + 1; innerIndex < bucket.size(); innerIndex++)
			{
				uint32_t innerNodeIndex = bucket[innerIndex];
				const Node& innerNode = nodes.mNodes[innerNodeIndex];
				//assert(innerNode.mRefCount > 0);

				// '==' Ignores ref count, so two node are identical even if thier ref counts differ.
				if(innerNode == outerNode)
				{
					auto mergeOpportunity = std::make_pair(outerNodeIndex, innerNodeIndex);
					result.insert(mergeOpportunity);
				}
			}
		}
	}

	return result;*/
}

template <typename PositionEnumeratorType, typename Function>
void applyFunction(Volume* volume, const Box3i& bounds, Function function)
{
	PositionEnumeratorType pe(bounds);

	do
	{
		auto functionResult = function(pe.x(), pe.y(), pe.z());
		volume->setVoxel(pe.x(), pe.y(), pe.z(), functionResult);
	}
	while(pe.next());
}

template <typename PositionEnumeratorType, typename Function>
std::pair<uint32_t, uint32_t> validateFunction(Volume* volume, const Box3i& bounds, Function function, uint64_t maxTests = std::numeric_limits<uint64_t>::max())
{
	uint32_t matches = 0;
	uint32_t mismatches = 0;
    uint64_t testCount = 0;

	PositionEnumeratorType pe(bounds);

	do
	{
		auto functionResult = function(pe.x(), pe.y(), pe.z());
		volume->voxel(pe.x(), pe.y(), pe.z()) == functionResult ? matches++ : mismatches++;
        testCount++;
	}
	while(pe.next() && (testCount < maxTests));

	return std::make_pair(matches, mismatches);
}

uint8_t randomMaterial(uint32_t /*x*/, uint32_t /*y*/, uint32_t /*z*/)
{
	return (rand() % 2 == 1) ? 2 : 7;
}

// FIXME - This should be a 3D checkerboard instead of a 2D one.
uint8_t checkerboard(uint32_t x, uint32_t y, uint32_t /*z*/)
{
	return (x % 2 == y % 2) ? 8 : 3;
}

bool testBounds()
{
	int sideLength = 256;
	std::unique_ptr<Volume> volume(new Volume);

	Box3i refResult;
	Box3i fullVolumeBox3i(Vector3i::filled(0), Vector3i::filled(sideLength - 1));

	/*Box3iSampler sampler(fullVolumeBox3i);

	for(int ct = 0; ct < 20; ct++)
	{
		Vector3i pos = sampler.next();
		volume->setVoxel(pos.x(), pos.y(), pos.z(), 1);
		refResult.accumulate(pos);
	}*/

	for(auto pos : Box3iSampler2(20, fullVolumeBox3i))
	{
		volume->setVoxel(pos.x(), pos.y(), pos.z(), 1);
		refResult.accumulate(pos);
	}

	const MaterialId externalMaterial = 0;
	Box3i bounds = computeBounds(*volume, externalMaterial);

	log_info("Lower     = ({},{},{})", bounds.lower().x(), bounds.lower().y(), bounds.lower().z());
	log_info("Ref Lower = ({},{},{})", refResult.lower().x(), refResult.lower().y(), refResult.lower().z());

	log_info("Upper     = ({},{},{})", bounds.upper().x(), bounds.upper().y(), bounds.upper().z());
	log_info("Ref Upper = ({},{},{})", refResult.upper().x(), refResult.upper().y(), refResult.upper().z());

	return true;
}

bool testBasics()
{
	std::pair<uint32_t, uint32_t> result;

	log_info("");
	log_info("Basic tests:");
	log_info("------------");


	// Create a volume for some simple tests
	int sideLength = 64;
	const Box3i bounds(Vector3i::filled(0), Vector3i::filled(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	// Volume should start empty
	result = validateFunction<RandomPositionEnumerator>(volume.get(), bounds, [](uint32_t, uint32_t, uint32_t) { return 0; });
	log_info("Empty volume node count = {}", volume->countNodes());
	log_info("Empty volume has {} matches and {} mismatches", result.first, result.second);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	// Fill it
	applyFunction<RandomPositionEnumerator>(volume.get(), bounds,  [](uint32_t, uint32_t, uint32_t) { return 5; });
	result = validateFunction<RandomPositionEnumerator>(volume.get(), bounds, [](uint32_t, uint32_t, uint32_t) { return 5; });
	log_info("\nFull volume node count = {}", volume->countNodes());
	log_info("Full volume has {} matches and {} mismatches", result.first, result.second);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	// Empty it again
	applyFunction<RandomPositionEnumerator>(volume.get(), bounds, [](uint32_t, uint32_t, uint32_t) { return 0; });
	result = validateFunction<RandomPositionEnumerator>(volume.get(), bounds, [](uint32_t, uint32_t, uint32_t) { return 0; });
	log_info("\nEmpty volume node count = {}", volume->countNodes());
	log_info("Empty volume has {} matches and {} mismatches", result.first, result.second);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	return true;
}

// A checkerboard pattern is interesting because it
// cannot be pruned, but it can be aggressively merged.
bool testCheckerboard()
{
	std::pair<uint32_t, uint32_t> result;

	log_info("");
	log_info("Checkerboard tests:");
	log_info("-------------------");

	// Create a volume for some simple tests
	int sideLength = 64;
	const Box3i bounds(Vector3i::filled(0), Vector3i::filled(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	// Fill it
	applyFunction<RandomPositionEnumerator>(volume.get(), bounds,  checkerboard);
	result = validateFunction<RandomPositionEnumerator>(volume.get(), bounds, checkerboard);

	log_info("Node count before bake = {}", volume->countNodes());
	volume->bake();
	log_info("Node count after bake = {}", volume->countNodes());

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	return true;
}

bool testRandomAccess()
{
	log_info("");
	log_info("Random access tests:");
	log_info("--------------------");

	// Create a volume for some simple tests
	int sideLength = 64;
	const Box3i bounds(Vector3i::filled(0), Vector3i::filled(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	Timer timer;
	for (int i = 0; i < 10; i++)
	{
		applyFunction<RandomPositionEnumerator>(volume.get(), bounds, randomMaterial);
	}

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	log_info("Node count = {}", volume->countNodes());

	log_info("Completed test in {} seconds", timer.elapsedTimeInSeconds());

	return true;
}

bool testSerialization()
{
	log_info("");
	log_info("Serialization test:");
	log_info("-------------------");

	// Create a volume for some simple tests
	int sideLength = 128;
	const Box3i bounds(Vector3i::filled(0), Vector3i::filled(sideLength - 1));
	Volume* volume = new Volume;

	// Write in simplex noise
	FractalNoise fractalNoise(7, 0, 0, 0);
	applyFunction<RandomPositionEnumerator>(volume, bounds, fractalNoise);

	volume->save("testSerialization.dag");
	delete volume;
	volume = new Volume("testSerialization.dag");

	auto validationResult = validateFunction<RandomPositionEnumerator>(volume, bounds, fractalNoise);

	// Test the result
	log_info("Serialization test gave {} matches and {} mismatches", validationResult.first, validationResult.second);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	delete volume;

	return true;
}

bool testFractalNoise()
{
	log_info("");
	log_info("Simplex noise tests:");
	log_info("--------------------");

	// Create a volume for some simple tests
	int sideLength = 128;
	const Box3i bounds(Vector3i::filled(0), Vector3i::filled(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	Timer timer;
	for (int i = 0; i < 10; i++)
	{
		// Each iteration samples from a differnt region of the fractal noise field to give different data each time.
		// We could have used 4D noise to achive this (sampling from different temporal slices) but evaluating 4D
		// noise is presumably slower, and we don't need the coherence which it can bring.
		int offset = i * 1000;

		// Write in simplex noise
		FractalNoise fractalNoise(7, offset, offset, offset);
		applyFunction<RandomPositionEnumerator>(volume.get(), bounds, fractalNoise);

		// Sometimes bake the octree
		if (i % 2 == 0)
		{
			volume->bake();
		}

		auto validationResult = validateFunction<RandomPositionEnumerator>(volume.get(), bounds, fractalNoise);

		// Test the result
		log_info("Simplex noise test gave {} matches and {} mismatches, node count = {}",
				 validationResult.first, validationResult.second, volume->countNodes());

		// Just another check that the node count is as expected
		if (i == 0) { assert(volume->countNodes() == 22817); }
	}

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	log_info("Completed test in {} seconds", timer.elapsedTimeInSeconds());

	return true;
}

bool testMerging()
{
	log_info("");
	log_info("Merge test:");
	log_info("-----------");

	// Create a volume for some simple tests
	int sideLength = 256;
	const Box3i bounds(Vector3i::filled(0), Vector3i::filled(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	Timer timer;

	// Write in simplex noise
	FractalNoise fractalNoise(9);

	MortonPositionEnumerator pe(bounds);

	do
	{
		auto result = fractalNoise(pe.x(), pe.y(), pe.z());
		if (result > 0)
		{
			volume->setVoxel(pe.x(), pe.y(), pe.z(), result);
		}
	} while (pe.next());

	//applyFunction<MortonPositionEnumerator>(volume.get(), bounds, fractalNoise);

	log_info("{} : Wrote data", timer.elapsedTimeInSeconds());

	volume->bake();
	log_info("{} : Baked", timer.elapsedTimeInSeconds());

	auto validationResult = validateFunction<RandomPositionEnumerator>(volume.get(), bounds, fractalNoise, 1000000);

	// Test the result
	log_info("Merge test gave {} matches and {} mismatches, node count = {}",
			 validationResult.first, validationResult.second, volume->countNodes());

	// Just another check that the node count is as expected
	assert(volume->countNodes() == 67065);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	log_info("{} : Completed", timer.elapsedTimeInSeconds());

	return true;
}

bool testSphere()
{
	return true;
}

bool testCSG()
{
	Volume building;
	building.load("../data/voxelized.dag");
	Volume shapes;
	shapes.load("../data/shapes.dag");

	building.addVolume(shapes);

	building.save("../data/csg.dag");

	uint8 outside_material;
	int32 lower_x, lower_y, lower_z, upper_x, upper_y, upper_z;
	cubiquity_estimate_bounds(&building, &outside_material, &lower_x, &lower_y, &lower_z, &upper_x, &upper_y, &upper_z);
	log_info("({},{},{}) ({},{},{})", lower_x, lower_y, lower_z, upper_x, upper_y, upper_z);

	int64_t histogram[256];
	cubiquity_compute_histogram(&building, histogram);
	for (int i = 0; i < 256; i++)
	{
		if (histogram[i] != 0) // Note that -1 can occur to indicate overflow
		{
			log_info("Material {}: {} voxels", static_cast<uint16_t>(i), histogram[i]);
		}
	}

	return true;
}

bool testVolume()
{
	srand(12345);

	// Slow test, so not usually run.
	/*if(!PositionEnumerator::test())
	{
		log_error("PositionEnumerator::test() failed!");(
	}*/

	testBounds();
	testBasics();
	//testCSG();
	testCheckerboard();
	testRandomAccess();
	testFractalNoise();
	testMerging();
	testSerialization();

	return true;
}
