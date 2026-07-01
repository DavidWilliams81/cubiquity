#ifndef CUBIQUITY_APP_COMBINE_H
#define CUBIQUITY_APP_COMBINE_H

#include <filesystem>

enum class Combiner {
    priority_union,  // Raw 3D array
};

bool combine_volumes(Combiner combiner,
               const std::filesystem::path& input_a_path,
               const std::filesystem::path& input_b_path,
                     std::filesystem::path  output_path);

#endif // CUBIQUITY_APP_COMBINE_H
