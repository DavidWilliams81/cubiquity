#include "import.h"

#include "base/logging.h"
#include "base/serialize.h"

#include <fstream>

using Cubiquity::Volume;

std::pair<std::unique_ptr<Cubiquity::Volume>, Metadata> import_from_bin(const std::filesystem::path& in_bin_path)
{
	std::unique_ptr<Cubiquity::Volume> volume =
		std::make_unique<Cubiquity::Volume>();

	std::filesystem::path metadata_path = in_bin_path;
	metadata_path.replace_extension(".txt");
	Metadata metadata(metadata_path);

    std::ifstream file(in_bin_path, std::ios::in | std::ios::binary);

    ivec3 dims = metadata.dimensions.value();
    for (int z = 0; z < dims.z; z++) {
        for (int y = 0; y < dims.y; y++) {
            for (int x = 0; x < dims.x; x++) {
                Cubiquity::MaterialId matId;
                file.read(reinterpret_cast<char*>(&matId), sizeof(matId));
                volume->setVoxel(x, y, z, matId);
            }
        }
    }

	return { std::move(volume), metadata };
}

bool import_from(ImportFormat           format,
           const std::filesystem::path& input_path,
                 std::filesystem::path  output_path)
{
	std::unique_ptr<Cubiquity::Volume> volume;
	Metadata metadata;

	switch (format)
	{
	case ImportFormat::bin:
		if (output_path.empty()) {
			output_path = input_path.filename().replace_extension(".dag");
		}
		std::tie(volume, metadata) = import_from_bin(input_path);
		break;
	default:
		log_error("Unknown export format");
	}

	saveVolume(output_path, *volume, metadata);

    return true;
}