/***************************************************************************************************
* Cubiquity - A micro-voxel engine for games and other interactive applications                    *
*                                                                                                  *
* Written in 2019 by David Williams                                                                *
*                                                                                                  *
* To the extent possible under law, the author(s) have dedicated all copyright and related and     *
* neighboring rights to this software to the public domain worldwide. This software is distributed *
* without any warranty.                                                                            *
*                                                                                                  *
* You should have received a copy of the CC0 Public Domain Dedication along with this software.    *
* If not, see http://creativecommons.org/publicdomain/zero/1.0/.                                   *
***************************************************************************************************/
#include "utility.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <set>
#include <string>

#include <iostream>
#include <sstream>

using namespace std;

namespace Cubiquity
{
	using namespace Internals;

	uint32_t Image::hash()
	{
		return Internals::murmurHash3(mData.get(), mWidth * mHeight * sizeof(Colour));
	}

	void Image::save(const std::string& filename)
	{
		// Get the file extension
		std::string ext = filename.substr(filename.find_last_of(".") + 1);
		std::transform(ext.begin(), ext.end(), ext.begin(),
			[](unsigned char c) { return std::tolower(c); });

		if (ext == "ppm")
		{
			std::fstream file;
			file = std::fstream(filename, std::ios::out | std::ios::binary);
			if (!file)
			{
				std::cerr << "Failed to open '" << filename << "' for writing" << std::endl;
			}

			//file.write((char*)&a, size * sizeof(unsigned long long));
			file << "P6\n";
			file << mWidth << " " << mHeight << "\n";
			file << "255\n";

			uint8* ptr = reinterpret_cast<uint8*>(mData.get());
			for (uint i = 0; i < mWidth * mHeight; i++)
			{
				file << *ptr;
				ptr++;
				file << *ptr;
				ptr++;
				file << *ptr;
				ptr++;
				ptr++;
			}

			file.close();
		}
		else
		{
			std::cerr << "Unrecognised extension '" << ext << "'" << std::endl;
		}
	}

	MaterialId gOutsideMaterialId;

	class NodeCounter
	{
	public:
		void operator()(NodeArray& /*nodeArray*/, uint32 nodeIndex, Box3i /*bounds*/)
		{
			mUniqueNodes.insert(nodeIndex);
		}
		uint64_t count() { return mUniqueNodes.size(); }
	private:
		std::set<uint32> mUniqueNodes;
	};

	uint64_t countNodes(Volume& volume)
	{
		NodeCounter nodeCounter;
		traverseNodes(volume, nodeCounter);
		return nodeCounter.count();
	}

	bool isInside(MaterialId matId)
	{
		return matId != gOutsideMaterialId;
	}

	class BoundsCalculator
	{
	public:
		BoundsCalculator(bool(*include)(MaterialId))
		{
			mIncludeFunc = include;
			mBounds.invalidate();
		}

		void operator()(NodeArray& /*nodeArray*/, uint32 nodeIndex, Box3i bounds)
		{
			if (nodeIndex < MaterialCount)
			{
				MaterialId matId = static_cast<MaterialId>(nodeIndex);
				if (mIncludeFunc(matId))
				{
					mBounds.accumulate(bounds.lower());
					mBounds.accumulate(bounds.upper());
				}
			}
		}

		Box3i bounds()
		{
			return mBounds;
		}

	private:
		Box3i mBounds;

		bool(*mIncludeFunc)(MaterialId);
	};

	Box3i computeBounds(Cubiquity::Volume& volume, bool(*include)(MaterialId))
	{
		BoundsCalculator boundsCalculator(include);
		traverseNodes(volume, boundsCalculator);
		return boundsCalculator.bounds();

	}

	std::pair<uint16_t, Box3i> estimateBounds(Volume& volume)
	{
		// Take a guess at what the outside material is by looking at the most extreme voxels
		// Note: This line could be shorter if we had a version of setVoxel which took a Vector3i, should add that.
		uint16_t mostNegativeVoxel = volume.voxel(std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest(), std::numeric_limits<int32_t>::lowest());
		uint16_t mostPositiveVoxel = volume.voxel(std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max());

		if (mostNegativeVoxel != mostPositiveVoxel)
		{
			std::cout << "Warning: Unable to accurately determine outside voxel" << std::endl;
		}
		uint16_t outsideMaterialId = mostPositiveVoxel;
		gOutsideMaterialId = outsideMaterialId;

		std::cout << "Outside material = " << outsideMaterialId << std::endl;

		Box3i bounds = computeBounds(volume, isInside);

		// If the bounds are the whole volume then we probably choose the wrong material to compute bounds for.
		// FIXME - Also check for Invalid here (max < min). A volume could be *all* empty or *not at all* empty.
		if (bounds == Box3i::invalid())
		{
			std::cout << "Warning: Bounds are invalid, volume filled with outside material (i.e. it is empty)." << std::endl;
		}

		if (bounds == Box3i::max())
		{
			std::cout << "Warning: Bounds are maxed out, did something go wrong?" << std::endl;
		}

		return std::make_pair(outsideMaterialId, bounds);
	}

	Histogram computeHistogram(Volume& volume, const Box3i& bounds)
	{
		Histogram histogram;
		//std::fill(histogram.begin(), histogram.end(), 0);

		for (int32 z = bounds.lower().z(); z <= bounds.upper().z(); z++)
		{
			for (int32 y = bounds.lower().y(); y <= bounds.upper().y(); y++)
			{
				for (int32 x = bounds.lower().x(); x <= bounds.upper().x(); x++)
				{
					histogram[volume.voxel(x, y, z)]++;
				}
			}
		}

		return histogram;
	}

	void printHistogram(const Histogram& histogram)
	{
		for (const auto& entry : histogram)
		{
			std::cout << "Material " << entry.first << ": " << entry.second << " voxels" << std::endl;
		}
		std::cout << "--------------------------------------------------------------------------------" << std::endl;
	}

	// Should match shader code
	float mix(float x, float y, float a)
	{
		return x * (1.0f - a) + y * a;
	}

	// See http://lolengine.net/blog/2013/07/27/rgb-to-hsv-in-glsl
	Vector3f hsv2rgb(Vector3f c)
	{
		Vector4f K = Vector4f(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
		Vector3f p = abs(fract(c.xxx() + K.xyz()) * 6.0f - K.www());
		return mix(K.xxx(), clamp(p - K.xxx(), 0.0f, 1.0f), c.y()) * c.z();
	}

	Vector3f colourFromMaterialId(MaterialId matId)
	{
		// Use the golden ratio to pick well-spaced colour in HSV space.
		// http://gofiguremath.org/natures-favorite-math/the-golden-ratio/the-golden-angle/
		float golden_ratio = 1.618;
		float angle = golden_ratio / 1 + golden_ratio;

		double intpart;
		float hue = std::modf(matId * angle, &intpart);

		// Saturation and value are fixed.
		float sat = 1.0;
		float val = 1.0;

		// Convert to RGB.
		Vector3f hsv = Vector3f(hue, sat, val);
		return hsv2rgb(hsv);
	}

	Colour rgbFromMaterialId(MaterialId matId)
	{
		uint8 red = (matId >> 8) & 0xF;
		uint8 green = (matId >> 4) & 0xF;
		uint8 blue = (matId) & 0xF;

		return Colour(red * 16.0f, green * 16.0f, blue * 16.0f);
	}

	void saveVolumeAsImages(Volume& volume, const std::string& dirName, ProgressMonitor* progMon)
	{
		Box3i bounds = estimateBounds(volume).second;

		uint32_t width = bounds.sideLength(0);
		uint32_t height = bounds.sideLength(1);

		if (progMon) { progMon->startTask("Saving volume as images"); }
		for (int z = bounds.lower().z(); z <= bounds.upper().z(); z += 1)
		{
			if (progMon) { progMon->setProgress(bounds.lower().z(), z, bounds.upper().z()); }

			// Note that the filenames start at zero (they are never negative). Using +/- symbols in the filenames is problematic,
			// at least because when sorting by name the OS lists '+' before'-', and also larger-magnitiude negative number after
			// smaller-magnitude negative numbers. This makes it more difficult to scroll through the slices.
			char filepath[256];
			std::snprintf(filepath, sizeof(filepath), "%s/%06d.ppm", dirName.c_str(), z - bounds.lower().z());

			Image image(width, height);
			for (int y = bounds.lower().y(); y <= bounds.upper().y(); y++)
			{
				for (int x = bounds.lower().x(); x <= bounds.upper().x(); x++)
				{
					MaterialId matId = volume.voxel(x, y, z);
					Colour colour = rgbFromMaterialId(matId);
					image.setPixel(x - bounds.lower().x(), y - bounds.lower().y(), colour);
				}
			}

			image.save(filepath);
		}

		if (progMon) { progMon->finishTask(); }
	}

	std::map<std::string, MaterialId> loadMtlFile(string filename)
	{
		std::map<std::string, MaterialId> materials;

		ifstream file(filename);

		string currentMaterialName;

		string line;
		while (std::getline(file, line))
		{
			istringstream iss(line);
			string element;
			if (!(iss >> element)) { continue; } // Blank line

			if (element == "newmtl")
			{
				iss >> currentMaterialName;
			}

			if (element == "Kd")
			{
				float r, g, b;
				iss >> r >> g >> b;

				uint16_t red = std::lround(r * 15.0f);
				uint16_t green = std::lround(g * 15.0f);
				uint16_t blue = std::lround(b * 15.0f);
				MaterialId materialId = (red << 8) | (green << 4) | blue;

				materials[currentMaterialName] = materialId;
			}
		}

		return materials;
	}

	// See https://en.wikipedia.org/wiki/Wavefront_.obj_file
	Geometry loadObjFile(const string& path, const string& filename)
	{
		Geometry geometry;

		ifstream file(path + "/" + filename);

		//string objName;
		vector<Vector3f> vertices;
		vertices.push_back(Vector3f(0.0f)); // Dummy value because obj file indices start at 1.
											//TriangleList triangles;

		std::map<std::string, MaterialId> materials;
		//string currentMaterialName;

		string line;
		uint32_t lineNo = 0;
		while (std::getline(file, line))
		{
			lineNo++;

			istringstream iss(line);
			string element;
			if (!(iss >> element)) { continue; } // Blank line

			if (element == "mtllib")
			{
				string materialFilename;
				iss >> materialFilename;
				materials = loadMtlFile(path + "/" + materialFilename);
			}

			if (element == "usemtl")
			{
				Object& object = geometry.back();
				object.subObjects.push_back(SubObject());
				SubObject& subObject = object.subObjects.back();

				std::string materialName;
				iss >> materialName;
				try
				{
					subObject.first = materials.at(materialName);
				}
				catch (const out_of_range& e)
				{
					subObject.first = 0x0fff;
				}
			}

			if (element == "o")
			{
				geometry.push_back(Object());
				Object& object = geometry.back();
				iss >> object.name;
			}

			if (element == "v")
			{
				float x, y, z;
				iss >> x >> y >> z;
				vertices.push_back(Vector3f(x, y, z));
			}

			if (element == "f")
			{
				// I have seen an OBJ file start without an 'o' or a 'usemtl' (or with them in the
				// wrong order). So if we find face definitions without an object/subobject having
				// been created then we create one on demand.
				if (geometry.empty())
				{
					geometry.push_back(Object());
				}
				Object& object = geometry.back();

				if (object.subObjects.empty())
				{
					object.subObjects.push_back(SubObject());
					object.subObjects.back().first = 0x0fff;
				}
				SubObject& subObject = object.subObjects.back();

				string i0Str, i1Str, i2Str;
				iss >> i0Str >> i1Str >> i2Str;

				i0Str = i0Str.substr(0, i0Str.find("/"));
				i1Str = i1Str.substr(0, i1Str.find("/"));
				i2Str = i2Str.substr(0, i2Str.find("/"));

				uint32_t i0 = stoul(i0Str);
				uint32_t i1 = stoul(i1Str);
				uint32_t i2 = stoul(i2Str);

				subObject.second.push_back(Triangle(vertices[i0], vertices[i1], vertices[i2]));

				uint32 temp;
				while (iss >> temp)
				{
					i1 = i2;
					i2 = temp;
					subObject.second.push_back(Triangle(vertices[i0], vertices[i1], vertices[i2]));
				}
			}
		}

		return geometry;
	}

	GaloisLFSR::GaloisLFSR(uint32_t mask, uint32_t startState)
		: mMask(mask), mState(startState) {}

	void GaloisLFSR::next()
	{
		uint32_t lsb = mState & 1;
		mState >>= 1;
		if (lsb)
		{
			mState ^= mMask;
		}
	}

	uint32_t GaloisLFSR::state()
	{
		return mState;
	}

	// Note that the LSFR will never hit zero. The period of a
	// maximal length LSFR is therefore 2^n-1 instead of 2^n.
	uint64_t GaloisLFSR::computePeriod()
	{
		uint64_t period = 0;
		uint32_t startState = mState;

		do
		{
			period++;
			next();
		} while (state() != startState);

		return period;
	}

	ShuffledSequence::ShuffledSequence(uint32 sequenceLength)
		: mSequenceLength(sequenceLength)
		, mGaloisLFSR(maximulLengthMask(logBase2(roundUpToPowerOf2(mSequenceLength + 1))))
	{
	}

	void ShuffledSequence::next()
	{
		do
		{
			mGaloisLFSR.next();
		} while (mGaloisLFSR.state() > mSequenceLength);
	}

	uint32_t ShuffledSequence::state()
	{
		// Subtract one to include zero in outputs.
		return mGaloisLFSR.state() - 1;
	}

	uint32_t ShuffledSequence::maximulLengthMask(int sizeInBits)
	{
		if (sizeInBits < 4 || sizeInBits > 32)
		{
			throw std::out_of_range("Invalid mask selection");
		}
		return MaximulLengthMasks[sizeInBits];
	}

	// These values are taken from here:
	//
	//     https://users.ece.cmu.edu/~koopman/lfsr/index.html
	//
	// For each bitsize the page provides a list of mask values which
	// correspond to maximul-length LSFRs. I've taken the first value
	// from each list and placed it in the array below.
	// 
	// The position in the array corresponds to the size (in bits) of
	// the LFSR. There is no such thing as a zero-bit version so that
	// entry is invalid, and the page above does not contain mask for
	// 1-3 bits. Therefore only elements 4-32 in the array are valid.
	const uint32_t ShuffledSequence::MaximulLengthMasks[33] =
	{
		0x00000000u, // Element zero is invalid

		0x00000000u,0x00000000u,0x00000000u,0x00000009u,
		0x00000012u,0x00000021u,0x00000041u,0x0000008Eu,
		0x00000108u,0x00000204u,0x00000402u,0x00000829u,
		0x0000100Du,0x00002015u,0x00004001u,0x00008016u,
		0x00010004u,0x00020013u,0x00040013u,0x00080004u,
		0x00100002u,0x00200001u,0x00400010u,0x0080000Du,
		0x01000004u,0x02000023u,0x04000013u,0x08000004u,
		0x10000002u,0x20000029u,0x40000004u,0x80000057u,
	};
}