#ifndef CUBIQUITY_TESTS_FRAMEWORK_H
#define CUBIQUITY_TESTS_FRAMEWORK_H

#include "rendering.h"

#include <iostream>

#define check(actual, expected) \
do \
{ \
	if (!(actual == expected)) \
	{ \
		std::cerr << std::endl <<\
			"********************************************************************************" \
			<< std::endl << std::endl << "\tCHECK FAILED!" << std::endl << \
			"\t\tActual result   : " << actual << std::endl << \
			"\t\tExpected result : " << expected << std::endl << \
			"\t\tLocation        : Line " << __LINE__ << " of " << __FILE__ << \
			std::endl << std::endl << \
			"********************************************************************************" \
			<< std::endl << std::endl; \
	} \
} while (0)

void saveVisibilityMaskAsImage(Cubiquity::VisibilityMask& visMask, const std::string& filename);

#endif