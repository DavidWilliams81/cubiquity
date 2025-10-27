#ifndef CUBIQUITY_APP_IMPORT_H
#define CUBIQUITY_APP_IMPORT_H

#include <filesystem>

enum class ImportFormat {
    bin,  // Raw 3D array
};

bool import_from(ImportFormat           format,
           const std::filesystem::path& input_path,
                 std::filesystem::path  output_path);

#endif // CUBIQUITY_APP_IMPORT_H
