#ifndef CUBIQUITY_TESTS_FRAMEWORK_H
#define CUBIQUITY_TESTS_FRAMEWORK_H

#include "base/logging.h"

#include "visibility.h"

#define check(actual, expected) \
do \
{ \
	if (!(actual == expected)) \
	{ \
		log_error("\n********************************************************************************"); \
		log_error("CHECK FAILED!"); \
		log_error("\tActual result   : {}", actual); \
		log_error("\tExpected result : {}", expected); \
		log_error("\tLocation        : Line {} of {}", __LINE__, __FILE__); \
		log_error("********************************************************************************\n"); \
	} \
} while (0)

void saveVisibilityMaskAsImage(Cubiquity::VisibilityMask& visMask, const std::string& filename);

#endif
