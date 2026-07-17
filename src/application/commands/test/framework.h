#ifndef CUBIQUITY_TESTS_FRAMEWORK_H
#define CUBIQUITY_TESTS_FRAMEWORK_H

#include "base/logging.h"

#include "extraction.h"

#define check(name, actual, expected) \
do \
{ \
	if (!(actual == expected)) \
	{ \
		log_error("\n********************************************************************************"); \
		log_error("CHECK FAILED!"); \
		log_error("\tName            : '{}'", name); \
		log_error("\tLocation        : Line {} of {}", __LINE__, __FILE__); \
		log_error("\tActual result   : {}", actual); \
		log_error("\tExpected result : {}", expected); \
		log_error("********************************************************************************\n"); \
	} \
	else \
	{ \
		log_info("Check '{}' passed", name); \
	} \
} while (0)

#endif
