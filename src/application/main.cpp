/*******************************************************************************
Cubiquity - A micro-voxel engine for games and other interactive applications

Written by David Williams

To the extent possible under law, the author(s) have dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

You should have received a copy of the CC0 Public Domain Dedication along with
this software. If not, see http://creativecommons.org/publicdomain/zero/1.0/.
*******************************************************************************/

// Cubiquity application includes
#include "license.h"
#include "base/logging.h"
#include "base/progress.h"
#include "commands/export/export.h"
#include "commands/generate/generate.h"
#include "commands/import/import.h"
#include "commands/test/test.h"
#include "commands/voxelize/voxelize.h"

// We can only support the 'view' command if SDL is available.
#ifdef CUBIQUITY_APP_ENABLE_VIEW
	#include "commands/view/view.h"
	// SDL should be included from the file defining main() so that it can
	// rename it to SDL_main, which is the preferred SDL model of operation.
	#include "SDL.h"
#endif // CUBIQUITY_APP_ENABLE_VIEW

// Cubiquity library
#include "base.h"

// External libraries
#include "CLI11.hpp"

// Cubiquity logging functions
void log_debug_func(const char* message)
{
	log_debug("{}", message);
}

void log_warning_func(const char* message)
{
	log_warning("{}", message);
}

int main(int argc, char** argv)
try
{
	// Connect Cubiquity to this application's logging/progress system
	Cubiquity::setLogDebugFunc(&log_debug_func);
	Cubiquity::setLogWarningFunc(&log_warning_func);
	Cubiquity::setProgressHandler(&cubiquityProgressHandler);

	// ================ Configure top-level application ================

	CLI::App app("Cubiquity micro-voxel engine by David Williams", "cubiquity");
	app.footer("Cubiquity is public domain software with some open source "
		       "dependencies.\n For details run 'cubiquity --license'");

	// Inherited by subcommands, and lets them pass unknown flags up a level.
	// This lets the user write e.g. 'cubiquity subcommand ... --verbose'
	// rather than 'cubiquity --verbose subcommand ...' (which feels clunky).
	app.fallthrough(true);

	// Needed for correct Unicode handling on Windows.
	argv = app.ensure_utf8(argv);

	// We may want to print the command line later as the user may
	// not have seen it (e.g. if run from an external program).
	std::string cmd_line(argv[0]); // Always valid
	for (int i = 1; i < argc; i++) {
		cmd_line = cmd_line + " " + argv[i];
	}

	// 'Pure' flags (no params)
	auto license = app.add_flag("--license", "Print license");
	auto quiet   = app.add_flag("--quiet",   "Only print warnings and errors");
	auto verbose = app.add_flag("--verbose", "Print extra debug information");

	// Options (with params)
	bool color_output{ true };
	app.add_option("--color-output", color_output, "Colorize output");

	// Constraints
	quiet->excludes(verbose);

	// Note that 'parse_complete_callback()' on the app gets
	// run before 'final_callback()' on the subcommands.
	app.parse_complete_callback([&]() {
		// Apply top-level flags
		set_color_enabled(color_output);
		if (*license) {
			printLicense();
			exit(EXIT_SUCCESS);
		}

		if (*quiet)   { set_verbosity(log_level::warning); }
		if (*verbose) { set_verbosity(log_level::debug);   }

		// Print command line in case user didn't type it themself
		log_debug("Running: {}", cmd_line);
		log_debug("Verbose output enabled"); // Only prints if actually enabled
	});

	// ==================== Configure subcommands ====================
	// Note that a command line can include multiple subcommands, e.g.
	// 'cubiquity voxelize ... view...'. Therefore subcommands need their own
	// copy of variables which can appear more than once, such as input/output
	// path, to avoid overwriting each others data.

	// ------------------------ Export ------------------------
	CLI::App* exp_cmd = app.add_subcommand("export",
		"Export a Cubiquity volume as the specified format");

	ExportFormat exp_fmt;
	std::map<std::string, ExportFormat> exp_fmt_map{
		{"bin",  ExportFormat::bin},
		{"pngs", ExportFormat::pngs},
		{"vox",  ExportFormat::vox}
	};
	exp_cmd->add_option("format", exp_fmt, "Output format to export as")
		->required()
		->transform(CLI::CheckedTransformer(exp_fmt_map, CLI::ignore_case));

	std::filesystem::path exp_in_path;
	exp_cmd->add_option("--input,input", exp_in_path,
		"Path to the '.dag' file to export")
		->required()
		->check(CLI::ExistingFile);

	// No default value as it will depend on the export file format. If not
	// specified it is also desirable to derive it from the input filename.
	std::filesystem::path exp_out_path;
	exp_cmd->add_option("--output,output", exp_out_path,
		"Path to the resulting file "
		"(derived from input path if not specified)");

	std::filesystem::path exp_out_meta_path;
	exp_cmd->add_option("--output-metadata", exp_out_meta_path,
		"Path to the resulting metadata file (when applicable)."
		"If not specified then this is derived from output path, "
		"unless writing to stdout (in which case it must be specified).");

	exp_cmd->final_callback([&]() {
		export_as(exp_fmt, exp_in_path, exp_out_path, exp_out_meta_path); });

	// ------------------------ Generate ------------------------
	CLI::App* gen_cmd = app.add_subcommand("generate",
		"Generate volume procedurally");

	Algorithm gen_algo;
	std::map<std::string, Algorithm> gen_algo_map{
		{"fractal_noise",  Algorithm::fractal_noise},
		{"menger_sponge",  Algorithm::menger_sponge},
		{"worley_noise",   Algorithm::worley_noise }
	};
	gen_cmd->add_option("algorithm", gen_algo, "Algorithm to use")
		->required()
		->transform(CLI::CheckedTransformer(gen_algo_map, CLI::ignore_case));

	std::filesystem::path gen_out_path = "output.dag";
	gen_cmd->add_option("--output,output", gen_out_path, "Output file name");

	int size = 500;
	gen_cmd->add_option("--size", size, "Volume side length in voxels");

	gen_cmd->final_callback([&]() {
		generateVolume(gen_algo, gen_out_path, size); });

	// ------------------------ Import ------------------------
	CLI::App* imp_cmd = app.add_subcommand("import",
		"Import a Cubiquity volume from the specified format");

	ImportFormat imp_fmt;
	std::map<std::string, ImportFormat> imp_fmt_map{
		{"bin",  ImportFormat::bin}
	};
	imp_cmd->add_option("format", imp_fmt, "Input format to import from")
		->required()
		->transform(CLI::CheckedTransformer(imp_fmt_map, CLI::ignore_case));

	std::filesystem::path imp_in_path;
	imp_cmd->add_option("--input,input", imp_in_path,
		"Path to the file to import")
		->required()
		->check(CLI::ExistingFile);

	// If not specified it will be derived from the input filename.
	std::filesystem::path imp_out_path;
	imp_cmd->add_option("--output,output", imp_out_path,
		"Path to the resulting '.dag' file "
		"(derived from input path if not specified)");

	imp_cmd->final_callback([&]() {
		import_from(imp_fmt, imp_in_path, imp_out_path); });

	// ------------------------ Test ------------------------
	CLI::App* test_cmd = app.add_subcommand("test", "Run tests");

	Test test_to_run{ Test::all };
	std::map<std::string, Test> test_map{
		{"all", Test::all},
		{"base", Test::base},
		{"volume", Test::volume},
		{"voxelization", Test::voxelization},
		{"raytracing", Test::raytracing}
	};
	test_cmd->add_option("test_to_run", test_to_run, "The test to run")
		    ->transform(CLI::CheckedTransformer(test_map, CLI::ignore_case));

	test_cmd->final_callback([&]() { test(test_to_run); });

	// ------------------------ View ------------------------
	CLI::App* view_cmd = app.add_subcommand("view", "View volume");

#ifdef CUBIQUITY_APP_ENABLE_VIEW
	std::filesystem::path view_in_path;
	view_cmd->add_option("--input,input", view_in_path,
			             "Path to the '.dag' file to view");

	ViewMode view_mode{ ViewMode::gpu_pathtracing };
	std::map<std::string, ViewMode> view_mode_map{
		{"instancing", ViewMode::instancing},
		{"cpu-pathtracing", ViewMode::cpu_pathtracing},
		{"gpu-pathtracing", ViewMode::gpu_pathtracing}
	};
	view_cmd->add_option("--mode", view_mode, "Rendering mode")
			->transform(CLI::CheckedTransformer(
				view_mode_map, CLI::ignore_case));

	view_cmd->final_callback([&]() { viewVolume(view_mode, view_in_path	); });
#else
	view_cmd->final_callback([&]() {
		log_error("Subcommand 'view' not available (built without SDL support)");
	});
#endif // CUBIQUITY_APP_ENABLE_VIEW

	// ------------------------ Voxelize ------------------------
	CLI::App* vox_cmd = app.add_subcommand("voxelize", "Voxelize mesh");

	std::filesystem::path vox_in_path;
	vox_cmd->add_option("--input,input", vox_in_path, "Input file name")
		   ->required()
		   ->check(CLI::ExistingFile);

	// No default value as we derive it from the input filename.
	std::filesystem::path vox_out_path;
	vox_cmd->add_option("--output", vox_out_path, "Output path. If not set "
			            "then this is derived from the input path but "
			            "with a different extension.");

	// Not std::optional as this hides the default value from the help
	// (and it's not really optional anyway as a default is provided).
	int vox_size = VoxelizeDefaultSize;
	auto vox_size_opt =
		vox_cmd->add_option("--size", vox_size,
			                "Desired length of the longest axis of the output "
			                "volume (in voxels)")
		       ->capture_default_str()
		       ->check(CLI::PositiveNumber);


	std::optional<float> vox_scale;
	vox_cmd->add_option("--scale", vox_scale,
		                "Scale factor applied to the input")
		   ->check(CLI::PositiveNumber)
		   ->excludes(vox_size_opt);

	vox_cmd->final_callback([&]() {
		// If the user provides a scale then they are excluded from providing a 
		// size, but a default size is still present. We must forcibly ignore 
		// it to adhere to the requirement that only one be specified.
		voxelize(vox_in_path, vox_out_path, vox_scale, 
			vox_scale ? std::optional<int>{} : vox_size); });

	// ==================== Parse and run ====================

	// If no arguments were provided, display help and exit
	if (argc == 1) {
		std::cout << app.help() << std::endl;
		return EXIT_SUCCESS;
	}

	// Parse the command line (modified expansion of CLI11_PARSE() macro)
	try {
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError& e) {
		// It can be useful to see the failed command line
		log_error_if(
			e.get_exit_code() != static_cast<int>(CLI::ExitCodes::Success),
			"Error parsing command line \"{}\"", cmd_line);

		// Redirect error messages to our own handlers.
		std::stringstream out, err;
		auto result = app.exit(e, out, err);
		log_info ("{}", out.str());
		log_error("{}", err.str());
		return result;
	}

	return EXIT_SUCCESS;
}
catch (const std::exception& e) {
	log_error("{}", e.what());
	return EXIT_FAILURE;
}
catch (...) {
	log_error("Unhandled exception!");
	return EXIT_FAILURE;
}
