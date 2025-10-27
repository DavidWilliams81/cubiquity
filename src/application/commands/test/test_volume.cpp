#include "test_volume.h"

#include "position_generator.h"

#include "base/bounds.h"
#include "base/noise.h"
#include "base/logging.h"
#include "base/random3d.h"


#include "cubiquity.h"
#include "utility.h"
#include "storage.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <thread>
#include <mutex>
#include <set>

using Cubiquity::Volume;

bool checkIntegrity(Volume& volume)
{
	const Cubiquity::Internals::NodeStore& nodes = Cubiquity::Internals::getNodes(volume).nodes();

	// Reserved nodes always contain dummy data.
	//for (u32 i = 0; i < MaterialCount; i++)
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
	//for (u32 nodeIndex = nodes.hashMapBegin() + 1; nodeIndex < nodes.hashMapEnd(); nodeIndex++)
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
			for (u32 material = 0; material < MaterialCount; material++)
			{
				if (nodes[nodeIndex].allChildrenAre(material))
					return false;
			}
		}*/
	}

	// Also check that no node points at itself.
	/*for (u32 nodeIndex = nodes.hashMapBegin(); nodeIndex < nodes.hashMapEnd(); nodeIndex++)
	{
		for (u32 childId = 0; childId < 8; childId++)
		{
			u32 childIndex = nodes[nodeIndex][childId];
			if (childIndex == nodeIndex)
				return false;
		}
	}*/

	// Fixme - should also check all nodes are pruned.

	// Lastly we re-compute ref counts for all nodes
	// and make sure they match the existing values.

	//FIXME - Used to use ref count, wht to do here?
	/*bool allRefCountsCorrect = true;
	u32* correctRefCounts = new u32[nodes.size()];
	memset(correctRefCounts, 0, sizeof(u32) * nodes.size());

	// Compute ref counts
	for (u32 nodeIndex = RootNodeIndex; nodeIndex < nodes.size(); nodeIndex++)
	{
		for(u32 childId = 0; childId < 8; childId++)
		{
			u32 childIndex = nodes[nodeIndex].child(childId);
			correctRefCounts[childIndex]++;
		}
	}

	// Check ref counts
	for(u32 nodeIndex = RootNodeIndex; nodeIndex < nodes.size(); nodeIndex++)
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

std::set< std::pair<u32, u32> > mergeOpportunities(Volume& volume)
{
	Cubiquity::Internals::NodeDAG& nodes = Cubiquity::Internals::getNodes(volume);

	std::set< std::pair<u32, u32> > result;

	// Early out until we fix use of ref count.
	return result;

	/*std::vector< std::vector<u32> > hashTable;

	hashTable.resize(100000);


	// Skip empty and full nodes, as they contain identical (but dummy) data giving a false merge opportunity.
	for (u32 index = nodes.mSharedEdits.mBegin; index <= nodes.mSharedEdits.mEnd; index++)
	{
		const Node& node = nodes.mNodes[index];

		//FIXME - Used to use ref count, wht to do here?
		//if (node.mRefCount == 0) continue;

		u32 hash = std::hash<Node>{}(node);
		hashTable[hash % hashTable.size()].push_back(index);
	}

	// Now find the matches in each bucket
	for (u32 bucketIndex = 0; bucketIndex < hashTable.size(); bucketIndex++)
	{
		const std::vector<u32>& bucket = hashTable[bucketIndex];

		for (u32 outerIndex = 0; outerIndex < bucket.size(); outerIndex++)
		{
			u32 outerNodeIndex = bucket[outerIndex];
			const Node& outerNode = nodes.mNodes[outerNodeIndex];
			//assert(outerNode.mRefCount > 0);

			for (u32 innerIndex = outerIndex + 1; innerIndex < bucket.size(); innerIndex++)
			{
				u32 innerNodeIndex = bucket[innerIndex];
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

// Fractal noise geometry with materials assigned based on Voronoi cells.
// Useful for producing detailed test volumes of potentially unlimited size. 
class FractalWorleyNoise
{
public:
	FractalWorleyNoise(int octaves,
		               int offset_x = 0, int offset_y = 0, int offset_z = 0)
		: m_octaves(octaves),
		  m_offset_x(offset_x), m_offset_y(offset_y), m_offset_z(offset_z) {
	}

	u8 operator()(int x, int y, int z)
	{
		// Fractal noise determines whether voxel is occupied
		if (fractal_noise(x + m_offset_x,
			              y + m_offset_y,
			              z + m_offset_z,
			              m_octaves) >= 0.0f) {

			// If occupied then assign material via Worley noise
			return worley_noise(x, y, z, 32, 2);
		}
		else {
			return 0;
		}
	}

private:
	int m_octaves  = 1;
	int m_offset_x = 0;
	int m_offset_y = 0;
	int m_offset_z = 0;
};

template <typename Function>
void applyFunction(Volume* volume,
	               const ivec3& lower, const ivec3& upper,
	               Function function,
	               bool shuffle = false, u16 seed = 0)
{
	PositionGenerator cycle(lower, upper, shuffle, seed);
	std::for_each(cycle.begin(), cycle.end(), [&](ivec3 pos) {
		auto functionResult = function(pos.x, pos.y, pos.z);
		volume->setVoxel(pos.x, pos.y, pos.z, functionResult);
	});
}

template <typename Function>
std::pair<u32, u32> validateFunction(Volume* volume,
	                                 const ivec3& lower, const ivec3& upper,
	                                 Function function,
	                                 bool shuffle = false, u16 seed = 0)
{
	u32 matches = 0;
	u32 mismatches = 0;

	PositionGenerator cycle(lower, upper, shuffle, seed);
	std::for_each(cycle.begin(), cycle.end(), [&](ivec3 pos) {
		auto functionResult = function(pos.x, pos.y, pos.z);
		volume->voxel(pos.x, pos.y, pos.z) == functionResult ?
			matches++ : mismatches++;
	});

	return std::make_pair(matches, mismatches);
}

u8 randomMaterial(u32 /*x*/, u32 /*y*/, u32 /*z*/)
{
	return (rand() % 2 == 1) ? 2 : 7;
}

// FIXME - This should be a 3D checkerboard instead of a 2D one.
u8 checkerboard(u32 x, u32 y, u32 /*z*/)
{
	return (x % 2 == y % 2) ? 8 : 3;
}

void print_positions(PositionGenerator pos_gen)
{
	std::vector<ivec3> buf;
	std::for_each(pos_gen.begin(), pos_gen.end(), [&](ivec3 pos) {
		buf.push_back(pos);
	});

	log_info("\t{}", buf[0]);
	log_info("\t{}", buf[1]);
	log_info("\t{}", buf[2]);
	log_info("\t{}", buf[3]);
	log_info("\t{}", buf[4]);
	log_info("\t...");
	log_info("\t{}", buf[buf.size() - 5]);
	log_info("\t{}", buf[buf.size() - 4]);
	log_info("\t{}", buf[buf.size() - 3]);
	log_info("\t{}", buf[buf.size() - 2]);
	log_info("\t{}", buf[buf.size() - 1]);
	log_info("\tTotal positions = {}", buf.size());
}

bool testPositionGenerator()
{
	log_info("Testing PositionGenerator");
	log_info("-------------------------");
	ivec3 lower = { -101, -103, -107 };
	ivec3 upper = { 109, 113, 127 };

	log_info("Sequential generator:");
	print_positions(PositionGenerator(lower, upper, false, 0));

	for (u16 seed = 0; seed <= 4; seed++) {
		log_info("\nShuffled generator (seed {}):", seed);
		print_positions(PositionGenerator(lower, upper, true, seed));
	}

	log_info("");
	return true;
}


bool testBounds()
{
	int sideLength = 256;
	std::unique_ptr<Volume> volume(new Volume);

	ivec3 ref_lower(1000000000);
	ivec3 ref_upper(-1000000000);

	uniform_vec3i_distribution random_vec3i(ivec3(0), ivec3(sideLength - 1));

	for(int ct = 0; ct < 20; ct++)
	{
		ivec3 pos = random_vec3i();
		volume->setVoxel(pos.x, pos.y, pos.z, 1);
		ref_lower = min(pos, ref_lower);
		ref_upper = max(pos, ref_upper);
	}

	auto [lower, upper] = find_bounds(*volume);

	log_info("Lower     = ({})", lower);
	log_info("Ref Lower = ({},{},{})", ref_lower.x, ref_lower.y, ref_lower.z);

	log_info("Upper     = ({})", upper);
	log_info("Ref Upper = ({},{},{})", ref_upper.x, ref_upper.y, ref_upper.z);

	return true;
}

bool testBasics()
{
	std::pair<u32, u32> result;

	log_info("");
	log_info("Basic tests:");
	log_info("------------");


	// Create a volume for some simple tests
	int sideLength = 64;
	const ivec3 lower(ivec3(0));
	const ivec3 upper(ivec3(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	// Volume should start empty
	result = validateFunction(volume.get(), lower, upper, [](u32, u32, u32) { return 0; });
	log_info("Empty volume node count = {}", volume->countNodes());
	log_info("Empty volume has {} matches and {} mismatches", result.first, result.second);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	// Fill it
	applyFunction(volume.get(), lower, upper,  [](u32, u32, u32) { return 5; });
	result = validateFunction(volume.get(), lower, upper, [](u32, u32, u32) { return 5; });
	log_info("\nFull volume node count = {}", volume->countNodes());
	log_info("Full volume has {} matches and {} mismatches", result.first, result.second);

	if(!checkIntegrity(*volume))
	{
		log_error("Integrity check failed!!!");
	}

	// Empty it again
	applyFunction(volume.get(), lower, upper, [](u32, u32, u32) { return 0; });
	result = validateFunction(volume.get(), lower, upper, [](u32, u32, u32) { return 0; });
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
	std::pair<u32, u32> result;

	log_info("");
	log_info("Checkerboard tests:");
	log_info("-------------------");

	// Create a volume for some simple tests
	int sideLength = 64;
	const ivec3 lower(ivec3(0));
	const ivec3 upper(ivec3(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	// Fill it
	applyFunction(volume.get(), lower, upper,  checkerboard);
	result = validateFunction(volume.get(), lower, upper, checkerboard);

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
	const ivec3 lower(ivec3(0));
	const ivec3 upper(ivec3(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	Cubiquity::Timer timer;
	for (int i = 0; i < 10; i++)
	{
		applyFunction(volume.get(), lower, upper, randomMaterial);
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
	const ivec3 lower(ivec3(0));
	const ivec3 upper(ivec3(sideLength - 1));
	Volume* volume = new Volume;

	// Write in simplex noise
	FractalWorleyNoise noise_func(7, 0, 0, 0);
	applyFunction(volume, lower, upper, noise_func);

	volume->save("testSerialization.dag");
	delete volume;
	volume = new Volume("testSerialization.dag");

	auto validationResult = validateFunction(volume, lower, upper, noise_func);

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
	const ivec3 lower(ivec3(0));
	const ivec3 upper(ivec3(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	Cubiquity::Timer timer;
	for (int i = 0; i < 10; i++)
	{
		// Each iteration samples from a differnt region of the fractal noise field to give different data each time.
		// We could have used 4D noise to achive this (sampling from different temporal slices) but evaluating 4D
		// noise is presumably slower, and we don't need the coherence which it can bring.
		int offset = i * 1000;

		// Write in simplex noise
		FractalWorleyNoise noise_func(6, offset, offset, offset);
		applyFunction(volume.get(), lower, upper, noise_func, true, 1);

		// Sometimes bake the octree
		if (i % 2 == 0)
		{
			volume->bake();
		}

		auto validationResult =
			validateFunction(volume.get(), lower, upper, noise_func, true, 2);

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

	//test_position_generator();

	// Create a volume for some simple tests
	int sideLength = 256;
	const ivec3 lower(ivec3(0));
	const ivec3 upper(ivec3(sideLength - 1));
	std::unique_ptr<Volume> volume(new Volume);

	Cubiquity::Timer timer;

	// Write in simplex noise
	FractalWorleyNoise noise_func(7);

	applyFunction(volume.get(), lower, upper, noise_func, true, 3);

	log_info("{} : Wrote data", timer.elapsedTimeInSeconds());

	volume->bake();
	log_info("{} : Baked", timer.elapsedTimeInSeconds());

	volume->save("fractalNoise.dag");

	auto validationResult =
		validateFunction(volume.get(), lower, upper, noise_func, true, 4);

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

	// Bounds may have changed as a result of CSG, so find and print them.
	auto [lower, upper] = find_bounds(building);
	log_info("({}) ({})", lower, upper);

	i64 histogram[256];
	cubiquity_compute_histogram(&building, histogram);
	for (int i = 0; i < 256; i++)
	{
		if (histogram[i] != 0) // Note that -1 can occur to indicate overflow
		{
			log_info("Material {}: {} voxels", static_cast<u16>(i), histogram[i]);
		}
	}

	return true;
}

bool testVolume()
{
	srand(12345);

	testPositionGenerator();
	testBounds();
	testBasics();
	testCSG();
	testCheckerboard();
	testRandomAccess();
	testFractalNoise();
	testMerging();
	testSerialization();

	return true;
}
