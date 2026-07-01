#include "combine.h"

#include "base/logging.h"
#include "base/progress.h"
#include "base/serialize.h"

#include <fstream>

using Cubiquity::Volume;

bool combine_volumes(Combiner combiner,
	const std::filesystem::path& input_a_path,
	const std::filesystem::path& input_b_path,
	std::filesystem::path  output_path)
{
	std::unique_ptr<Cubiquity::Volume> volume;
	Metadata metadata;

	auto [volume_a, metadata_a] = loadVolume(input_a_path.string());
	auto [volume_b, metadata_b] = loadVolume(input_b_path.string());

	volume_a->combine(*volume_b, [](Cubiquity::MaterialId a, Cubiquity::MaterialId b) -> Cubiquity::MaterialId {
		if (a != 0 && b != 0) return (a > b) ? a : b; // Max label intersection
		return 0;
	});

	saveVolume(output_path, *volume_a, metadata_a);

    return true;
}
