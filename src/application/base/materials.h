#ifndef CUBIQUITY_APP_TYPES_H
#define CUBIQUITY_APP_TYPES_H

#include "base/logging.h"

#include "storage.h"

#include <array>
#include <filesystem>
#include <vector>

typedef std::array<float, 3> Col; // To be renamed when we can avoid conflicts.;

class MaterialSet
{
public:
	typedef std::array<Col, Cubiquity::Internals::MaterialCount> ColourArray;

public:
	MaterialSet();

	const Col& operator[](int index) const;

	void clear();
	const ColourArray& data() const;
	Cubiquity::MaterialId findOrInsert(const Col& material);

	void load(const std::filesystem::path& materialPath);
	void save(const std::filesystem::path& materialPath);

private:
	ColourArray mData;
	ColourArray::iterator mNextFreeSlot;
};

std::filesystem::path getMaterialsPath(const std::filesystem::path& volumePath);

#endif // CUBIQUITY_APP_TYPES_H
