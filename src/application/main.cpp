// Cubiquity application includes
#include "base/logging.h"
#include "base/progress.h"
#include "commands/export/export.h"
#include "commands/generate/generate.h"
#include "commands/test/test.h"
#include "commands/voxelise/voxelize.h"

// We can only support the 'view' command if SDL is available.
#ifdef CUBIQUITY_APP_ENABLE_VIEW
	#include "commands/view/view.h"
	// SDL should be included from the file defining main() so that it can
	// rename it to SDL_main, which is the preferred SDL model of operation.
	#include "SDL.h"
#else
	// Fallback implentation.
	bool viewVolume(const flags::args& args)
	{
		log_error("Subcommand 'view' not available (built without SDL support)");
		return true;
	}
#endif // CUBIQUITY_APP_ENABLE_VIEW

// Cubiquity library
#include "base.h"

// External libraries
#include "flags.h"

// System libraries
#include <map>

using namespace Cubiquity;

void log_debug_func(const char* message)
{
    log_debug("{}", message);
}

void log_warning_func(const char* message)
{
    log_warning("{}", message);
}

void printUsageAndExit() {
	print("Usage: cubiquity <command> [--quiet] [--verbose] ...\n");
	exit(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
	typedef bool (*CommandPtr)(const flags::args&);
	typedef std::map<std::string, CommandPtr> CommandMap;

	CommandMap commands = {
		{ "export",     &exportVolume },
		{ "generate",   &generateVolume },
		{ "test",       &test },
		{ "view",       &viewVolume },
		{ "voxelize",   &voxelize },
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
		log_error("Quiet and Verbose modes cannot be used together!");
		printUsageAndExit();
	}
	else if (quiet) {
		set_verbosity(log_level::warning);
	}
	else if (verbose) {
		set_verbosity(log_level::debug);
	}
	log_debug("Verbose output enabled");

	// Hook up Cubiquity to the applications logging system
	Cubiquity::setLogDebugFunc(&log_debug_func);
	Cubiquity::setLogWarningFunc(&log_warning_func);

	Cubiquity::setProgressHandler(&cubiquityProgressHandler);

	// Execute the requested command
	const std::string commandName { args.positional().at(0) };
	CommandMap::iterator commandIter = commands.find(commandName);
	if (commandIter != commands.end())
	{
		CommandPtr commandPtr = commandIter->second;
		if (!(*commandPtr)(args)) {
			log_error("Failed to execute command '{}'", commandName);
		}
	}
	else
	{
		log_error("Unrecognised command '{}'", commandName);
	}

	return EXIT_SUCCESS;
}
