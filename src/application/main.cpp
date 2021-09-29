#include "base/logging.h"
#include "base/progress.h"
#include "commands/export/export.h"
#include "commands/test/test.h"
#include "commands/view/view.h"
#include "commands/voxelise/voxelize.h"

#include "flags.h"

#include "base.h"

#include <map>

using namespace Cubiquity;

void printUsageAndExit() {
	log(Info, "Usage: cubiquity <command> [--quiet] [--verbose] ...");
	exit(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
	typedef bool (*CommandPtr)(const flags::args&);
	typedef std::map<std::string, CommandPtr> CommandMap;

	CommandMap commands = {
		{ "export",   &exportVolume },
		{ "test",     &test },
		{ "view",     &viewVolume },
		{ "voxelise", &voxelizeMesh },
	};

	const flags::args args(argc, argv);
	// Note that flags library does not seem to set exe name as zeroth 'parameter'
	if (args.positional().empty())
	{
		printUsageAndExit();
	}

	// Configure logging verbosity (could also add a 'quiet' mode?).
	const auto quiet = args.get<bool>("quiet", false);
	const auto verbose = args.get<bool>("verbose", false);
	if (quiet && verbose) {
		log(Error, "Quiet and Verbose modes cannot be used together!");
		printUsageAndExit();
	}
	else if (quiet) {
		setVerbosity(Warning);
	}
	else if (verbose) {
		setVerbosity(Debug);
	}
	log(Debug, "Verbose output enabled");

	// Hook up Cubiquity to the applications logging system
	Cubiquity::setLogHandler(&cubiquityLogHandler);

	Cubiquity::setProgressHandler(&cubiquityProgressHandler);

	// Execute the requested command
	const std::string commandName { args.positional().at(0) };
	CommandMap::iterator commandIter = commands.find(commandName);
	if (commandIter != commands.end())
	{
		CommandPtr commandPtr = commandIter->second;
		if (!(*commandPtr)(args)) {
			log(Error, "Failed to execute command \'%s\'.", commandName.c_str());
		}
	}
	else
	{
		log(Error, "Unrecognised command \'%s\'.", commandName.c_str());
	}

	return EXIT_SUCCESS;
}
