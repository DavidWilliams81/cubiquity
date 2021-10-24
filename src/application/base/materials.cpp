#include "materials.h"

#include "storage.h"

#include <algorithm>
#include <fstream>

using namespace Cubiquity;

std::filesystem::path getMaterialsPath(const std::filesystem::path& volumePath)
{
	std::filesystem::path materialsPath = volumePath;
	materialsPath.replace_extension(".mat");
	return materialsPath;
}

MaterialSet::MaterialSet()
{
	clear();
}

const Col& MaterialSet::operator[](int index) const
{
	return mData[index];
}

void MaterialSet::clear()
{
	Col white = { 1.0f, 1.0f, 1.0f };
	std::fill(mData.begin(), mData.end(), white);
	mNextFreeSlot = mData.begin() + 1; // Don't write to zeroth slot later.
}

const MaterialSet::ColourArray& MaterialSet::data() const
{
	return mData;
}


// Note - It is likely that succussive calls to this functions will be searching for
// the same material (e.g. from a mesh with many triangles having the same material).
// If performance is an issue it might be worth caching the last returned result.
Cubiquity::MaterialId MaterialSet::findOrInsert(const Col& material)
{
	Cubiquity::MaterialId matId = 0; // Default state indicates failure

	// Compare colours with some tolerance.
	auto matchesInput = [material](const Col& matToTest)
	{
		float threshold = 0.01;
		return
			std::abs(material[0] - matToTest[0]) < threshold &&
			std::abs(material[1] - matToTest[1]) < threshold &&
			std::abs(material[2] - matToTest[2]) < threshold;
	};
	
	// Find an existing matching entry.
	auto iter = std::find_if(mData.begin() + 1, mNextFreeSlot, matchesInput); // Skip zeroth element

	if (iter == mNextFreeSlot) // Not found
	{
		// Insert it if there is space.
		if (mNextFreeSlot < mData.end())
		{
			*mNextFreeSlot = material;
			matId = mNextFreeSlot - mData.begin(); // Convert iterator to index.
			mNextFreeSlot++;
		}
	}
	else // Found material in existing slot
	{
		matId = iter - mData.begin(); // Convert iterator to index.
	}

	return matId; // If matId is still zero this indicates failure
}

void MaterialSet::load(const std::filesystem::path& materialPath)
{
	clear();
	std::string line;
	std::ifstream materialFile(materialPath);
	while (std::getline(materialFile, line))
	{
		std::istringstream iss(line);
		uint32 matId; // Type avoids being treated as ASCII
		iss >> matId;
		iss >> mData[matId][0];
		iss >> mData[matId][1];
		iss >> mData[matId][2];
	}

	// We don't support letting the user insert after loading as it complicates things.
	mNextFreeSlot = mData.end();
}

void MaterialSet::save(const std::filesystem::path& materialPath)
{
	std::ofstream matFile(materialPath);
	for (int i = 1; i < mData.size(); i++)
	{
		matFile << i << " " << mData[i][0] << " "
			<< mData[i][1] << " " << mData[i][2] << std::endl;
	}
}