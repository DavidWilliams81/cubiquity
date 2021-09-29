#ifndef CUBIQUITY_TESTS_FRAMEWORK_H
#define CUBIQUITY_TESTS_FRAMEWORK_H

#include "base/logging.h"

#include "rendering.h"

#include <sstream>
#include <iostream>

#define check(actual, expected) \
do \
{ \
	if (!(actual == expected)) \
	{ \
		std::stringstream ss; \
		ss << std::endl << \
			"********************************************************************************" \
			<< std::endl << "CHECK FAILED!" << std::endl << \
			"\tActual result   : " << actual << std::endl << \
			"\tExpected result : " << expected << std::endl << \
			"\tLocation        : Line " << __LINE__ << " of " << __FILE__ << std::endl << \
			"********************************************************************************" \
			<< std::endl; \
		log(Error, ss.str().c_str()); \
	} \
} while (0)

void saveVisibilityMaskAsImage(Cubiquity::VisibilityMask& visMask, const std::string& filename);

#endif
