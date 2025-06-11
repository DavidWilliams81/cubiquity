#ifndef CUBIQUITY_TEST_H
#define CUBIQUITY_TEST_H

#include <string>

enum class Test {
	all,
	base,
	volume,
	voxelization,
	raytracing
};

//inline constexpr Tests operator&(Tests x, Tests y) {
//	return static_cast<Tests> (static_cast<int>(x) & static_cast<int>(y));
//}

bool test(const Test& test);

#endif // CUBIQUITY_TEST_H
