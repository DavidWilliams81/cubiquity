#include "import.h"

#include "base/logging.h"
#include "base/progress.h"
#include "base/serialize.h"

#include "utility.h"

#include <fstream>

using Cubiquity::Volume;

std::pair<std::unique_ptr<Cubiquity::Volume>, Metadata> import_from_bin(
	const std::filesystem::path& in_bin_path,
	const std::filesystem::path& metadata_path)
{
	std::unique_ptr<Cubiquity::Volume> volume =
		std::make_unique<Cubiquity::Volume>();

    //std::ifstream file = make_safe_ifstream(
	//	in_bin_path.string(), std::ios::in | std::ios::binary);

	InputHandle file(in_bin_path.string());

	Metadata metadata(metadata_path);

	// Copy dimensions for the loop, after which they are no longer wanted for the DAG.
    ivec3 dims = metadata.dimensions.value();
	metadata.dimensions.reset();

	// FIXME - Fill wioth correctg background colour
	Cubiquity::MaterialId background = 0;

	std::vector< Cubiquity::MaterialId> buffer(dims.x);
    for (int z = 0; z < dims.z; z++) {
        for (int y = 0; y < dims.y; y++) {
			file->read(reinterpret_cast<char*>(buffer.data()),
				buffer.size() * sizeof(Cubiquity::MaterialId));
            for (int x = 0; x < dims.x; x++) {
				if (buffer[x] != background) {
					volume->setVoxel(x, y, z, buffer[x]);
				}
            }
        }

		auto unshared_node_count = volume->mNodeStore.unsharedNodesEnd() - volume->mNodeStore.unsharedNodesBegin();
		if ((z == dims.z - 1) || (unshared_node_count > 100000000)) {
			volume->bake();
		}

        // A bit cheeky, but we can directly call our Cubiquity progress handling code for progress bar.
		cubiquityProgressHandler("Importing volume from raw 3D array",
			0, z, dims.z-1);
    }

	return { std::move(volume), metadata };
}

bool import_from(ImportFormat           format,
           const std::filesystem::path& input_path,
                 std::filesystem::path& input_metadata_path,
                 std::filesystem::path  output_path)
{
	std::unique_ptr<Cubiquity::Volume> volume;
	Metadata metadata;

	switch (format)
	{
	case ImportFormat::bin:
		if(input_path == "-") {
			if (input_metadata_path.empty()) {
				throw std::runtime_error(
					"Input metadata path must be specified when reading from stdin");
			}
		} else {
			if (input_metadata_path.empty()) {
				// We use the .txt extension for metadata associated with a .bin file partly
				// because it avoids a name conflict when importing/exporting a .bin to/from
				// a .dag, and partly to make it obvious that the metadata for a .bin is
				// human-readable (I'm not sure everyone knows what a .toml extension means).
				input_metadata_path = input_path;
				input_metadata_path.replace_extension(".txt");
			}
		}
		if (output_path.empty()) {
			output_path = input_path.filename().replace_extension(".dag");
		}
		std::tie(volume, metadata) = import_from_bin(input_path, input_metadata_path);
		break;
	default:
		log_error("Unknown import format");
	}

	saveVolume(output_path, *volume, metadata);

    return true;
}
