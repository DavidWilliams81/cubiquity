#include "test_base.h"

#include "storage.h"

#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>

using namespace Cubiquity;
using namespace std;

bool testBase()
{
	// Test utility function 
	assert(Internals::roundUpToPowerOf2(1023) == 1024);
	assert(Internals::roundUpToPowerOf2(1024) == 1024);

	// Test integer hashing
	assert(Internals::mix(42) == 0x087fcd5c);

	// Test data hashing
	std::vector<int> data(1000000);
	std::iota (std::begin(data), std::end(data), 0);
	uint32_t dataHash = Internals::murmurHash3(data.data(), sizeof(data[0]) * data.size(), 123);
	assert(dataHash == 0x6e894d1d);

	return true;
}
