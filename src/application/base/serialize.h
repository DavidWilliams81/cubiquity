#ifndef CUBIQUITY_PATHS_H
#define CUBIQUITY_PATHS_H

#include "storage.h"

#include "base/logging.h"
#include "metadata.h"

#include <filesystem>
#include <fstream>
#include <memory>
#include <utility>

class OutputHandle
{
public:
	OutputHandle(const std::string& path, bool force = false) {
		namespace fs = std::filesystem;

		if (path == "-") {
			// NOTE: Untested on Windows (may require setting cout to binary mode)
			log_debug("Writing to stdout");
			use_stdout = true;
		} else {
			std::ios::openmode mode = std::ios::out | std::ios::binary;
			if (fs::is_fifo(path)) {
				log_debug("Writing to FIFO {} (may block waiting for reader...)", path);
				m_stream = std::make_unique<std::ofstream>(path, mode |= std::ios::app);
				log_debug("Done constructing stream");
			} else if (fs::exists(path) && force == false) {
				throw std::runtime_error(fmt::format(
					"Path '{}' already exists (and force flag not specified): ", path));
			} else {
				log_debug("Writing to regular file {}", path);
				m_stream = std::make_unique<std::ofstream>(path, mode |= std::ios::trunc);
			}

			if(m_stream->bad() || m_stream->is_open() == false) {
				throw std::runtime_error("Error writing to  " + path);
			}
		}
	}

	// Accessor
	std::ostream& get() {
		assert(((use_stdout && m_stream) == false) && "Cannot both be active");
		return use_stdout ? std::cout : *m_stream;
	}

	std::ostream* operator->() { return &get(); }
	std::ostream& operator*()  { return get(); }

private:
	bool use_stdout = false;
	std::unique_ptr<std::ofstream> m_stream;
};

bool checkInputFileIsValid(const std::filesystem::path& inputFile);
bool checkOutputDirIsValid(const std::filesystem::path& outputDir);

// FIXME - I'd rather return the volume by value but
// I need to make a working move constructor first.
std::pair<std::unique_ptr<Cubiquity::Volume>, Metadata>
    loadVolume(const std::filesystem::path& vol_path);

void saveVolume(const std::filesystem::path& volume_path, 
                Cubiquity::Volume& volume, Metadata& metadata);

#endif // CUBIQUITY_PATHS_H


