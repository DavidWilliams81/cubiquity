#include "test_base.h"

#include "framework.h"

#include "storage.h"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using namespace Cubiquity;
using namespace std;

uint64 mixPeriod()
{
	uint32 initialState = 42;
	uint32 state = initialState;
	uint64 steps = 0;
	do
	{
		state = mixBits(state);
		steps++;
	} while (state != initialState);

	return steps;
}

bool testBase()
{
	// Test utility function 
	assert(Internals::roundUpToPowerOf2(1023) == 1024);
	assert(Internals::roundUpToPowerOf2(1024) == 1024);

	// Test integer hashing
	assert(Internals::mixBits(42) == 0x087fcd5c);


	// Check the period of the mix function. In the future
	// I hope to find one with only one cycle (max period).
	check(mixPeriod(), UINT64_C(722985151));
	//check(mixPeriod(), UINT64_C(0x100000000)); // This would be better

	// Test data hashing
	std::vector<int> data(1000000);
	std::iota (std::begin(data), std::end(data), 0);
	uint32_t dataHash = Internals::murmurHash3(data.data(), sizeof(data[0]) * data.size(), 123);
	assert(dataHash == 0x6e894d1d);

	return true;
}
