#include "test_utility.h"

#include "utility.h"

#include "framework.h"

#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>

using namespace Cubiquity;
using namespace std;

void testShuffledSequence(uint sequenceLength)
{
	ShuffledSequence sequence(sequenceLength);

	// Although the sequence is shuffled it always starts at zero. Knowing where
	// we start can be useful for testing when the sequence has wrapped around.
	// We could add a user-specified start point in the future if desired.
	uint32 firstValue = sequence.state();
	check(firstValue, 0);

	std::vector<bool> generated(sequenceLength);
	for (uint ct = 0; ct < sequenceLength; ct++)
	{
		uint32 value = sequence.state();
		sequence.next();

		check(generated[value], false); // Ensure we don't hit a number twice
		generated[value] = true;
	}

	// Ensure every number gets hit
	check(std::all_of(generated.begin(), generated.end(), [](bool v) { return v; }), true);
}

void testShuffledSequences()
{
	testShuffledSequence(14);
	testShuffledSequence(15);
	testShuffledSequence(16);
	testShuffledSequence(17);

	testShuffledSequence(65534);
	testShuffledSequence(65535);
	testShuffledSequence(65536);
	testShuffledSequence(65537);
}

void testUtility()
{
	// Test ShuffledSequence function
	testShuffledSequences();
}
