#include "test_base.h"

#include "framework.h"

#include <cassert>
#include <cstdint>

using namespace std;
using namespace Cubiquity::Internals;

bool testBase()
{
	// Test utility function
	check("is_power_of_2(1023)", is_power_of_2(1023), false);
	check("is_power_of_2(1024)", is_power_of_2(1024), true);
	check("log_base_2(1024)", log_base_2(1024), 10);
	check("next_power_of_2(1023)", next_power_of_2(1023), 1024);
	check("next_power_of_2(1024)", next_power_of_2(1024), 1024);

	// Test integer hashing
	check("bit_mix(0x0123456789abcdef)", bit_mix(0x0123456789abcdef), 0x960cbea3c15f985a);

	// Test data hashing
	check("fnv1a(\"hello\")", fnv1a("hello", 5    ), 0xa430d84680aabd0b);
	check("fnv1a(\"hello\", seed 42)", fnv1a("hello", 5, 42), 0x2dc456e562e2b4a2);

	// Test hash combining
    std::size_t running_hash = 0;
	running_hash = hash_combine(running_hash, hash_value(1.0f));
	running_hash = hash_combine(running_hash, hash_value(2.0f));
	running_hash = hash_combine(running_hash, hash_value(3.0f));
    check("Combined hash", running_hash, UINT64_C(9174456649812731690));

	// Test that +0.0f and -0.0f hash the same
	check("Zero-canonicalization", hash_value(0.0f), hash_value(-0.0f));

	return true;
}
