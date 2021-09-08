#include "mesh.h"

#include <fstream>
#include <sstream>

using namespace std;
using namespace Cubiquity;

MaterialId colourToMaterialId(std::array<float, 3> colour)
{
	// Clamp to a valid range. Note that we can't use zero as black because
	// that represents empty space. I just boost the level here by setting
	// the minimum to 1, but could look for a better approch in the future.
	uint16_t r = std::clamp(std::lround(colour[0] * 15.0f), 1L, 15L);
	uint16_t g = std::clamp(std::lround(colour[1] * 15.0f), 1L, 15L);
	uint16_t b = std::clamp(std::lround(colour[2] * 15.0f), 1L, 15L);
	MaterialId result = (r << 8) | (g << 4) | b;
	return result;
}

// https://en.wikipedia.org/wiki/Wavefront_.obj_file#Material_template_library
// and (official spec?) http://paulbourke.net/dataformats/mtl/
std::map<std::string, MaterialId> loadMtlFile(const string& filename)
{
	ifstream file(filename);
	std::map<std::string, MaterialId> materials;
	std::map<std::string, MaterialId>::iterator currentMaterial;
	bool foundDiffuse = true;

	auto finishMaterial = [&]()
	{ 
		if (!foundDiffuse)
		{
			std::cerr << "Warning: No diffuse colour found for material \""
				<< currentMaterial->first << "\"" << std::endl;
		}
	};

	string line;
	uint lineNo = 0;
	while (std::getline(file, line))
	{
		lineNo++;
		istringstream iss(line);
		string element;
		if (!(iss >> element)) { continue; } // Blank line

		if (element == "newmtl")
		{
			// Finish the current material
			finishMaterial();

			// Start the new one
			string materialName;
			iss >> materialName;

			// Build the material with default values and add it to the map.
			// It will be updated later as we encounter relevent tags.
			MaterialId materialId = colourToMaterialId({1.0f, 1.0f, 1.0f});
			foundDiffuse = false;
			auto result = materials.insert(std::make_pair(materialName, materialId));
			if (!result.second)
			{
				std::cerr << "Warning: Material \"" << materialName
					<< "\" is already defined in \"" << filename << "\"" << std::endl;
			}
			// Update even if the material couldn't be inserted. Pointing at the duplicate
			// is probably better than whatever we were pointing at before (if anything).
			currentMaterial = result.first;
		}

		if (element == "Kd")
		{
			std::array<float, 3> diffuse;
			iss >> diffuse[0] >> diffuse[1]>> diffuse[2];
			currentMaterial->second = colourToMaterialId(diffuse);
			foundDiffuse = true;
		}

		if (iss.fail())
		{
			std::cerr << "Error: Failed to parse line " << lineNo
				<< " (\"" << line << "\")" << std::endl;
		}
	}

	finishMaterial();

	if (materials.empty())
	{
		std::cerr << "Warning: No materials found in \"" << filename << "\"" << std::endl;
	}

	return materials;
}

// See https://en.wikipedia.org/wiki/Wavefront_.obj_file 
// and https://www.fileformat.info/format/wavefrontobj/egff.htm
// and (official spec?) http://paulbourke.net/dataformats/obj/ 
std::list<Object> loadObjFile(const string& path, const string& filename)
{
	ifstream file(path + "/" + filename);

	std::list<Object> objects;
	std::vector<Vertex> vertices;
	std::map<std::string, MaterialId> materials;

	// Insert an (unused) dummy vertex because obj file indices start at 1.
	vertices.push_back({ 0,0,0 }); 

	// I have seen an OBJ file start without an 'o' or a 'usemtl' (or 
	// with them in the wrong order). Create a dummy object to hold these.
	MaterialId activeMaterialId = colourToMaterialId({1.0f, 1.0f, 1.0f});
	objects.push_back(Object());

	string line;
	uint lineNo = 0;
	uint tooManyComponentsLineNo = 0; // For reporting warning

	// Parse the file one line at a time
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
			std::map < std::string, MaterialId > localMaterials = loadMtlFile(path + "/" + materialFilename);
			for (const auto& entry : localMaterials)
			{
				if (!materials.insert(entry).second)
				{
					std::cerr << "Warning: Material \"" << entry.first 
						<< "\" was already defined in another material file." << std::endl;
				}
			}
		}

		if (element == "usemtl")
		{
			std::string materialName;;
			iss >> materialName;

			auto materialIter = materials.find(materialName);
			activeMaterialId = materialIter != materials.end() ?
				materialIter->second:  colourToMaterialId({ 1.0f, 1.0f, 1.0f });
		}

		if (element == "o")
		{
			// The spec require objects to always have a name. I'm not 
			// sure if duplicates are allowed, but we can handle them. 
			Object object;
			iss >> object.name;
			objects.push_back(object);
		}

		if (element == "v")
		{
			Vertex vertex;
			iss >> vertex[0] >> vertex[1] >> vertex[2];
			vertices.push_back(vertex);

			if (!tooManyComponentsLineNo && !iss.eof()) // Only warn the first time
			{
				// Note that the user might not have control over any extra components
				// as they might come from their modelling package or exporter/converter.
				// So the warning gets reported once (later) to avoid spamming the user.
				tooManyComponentsLineNo = lineNo;
			}
		}

		if (element == "f")
		{
			std::vector<int> vertexIndices;
			while (!iss.eof())
			{
				std::string indexString;
				iss >> indexString;
				indexString = indexString.substr(0, indexString.find("/"));
				int index = stol(indexString);

				// If we wanted to support relative vertex indices then I think we would
				// have to resolve them here. It's probably not difficult, but I don't
				// have a sample file to test against so won't worry about it for now.
				if (index > 0) // Valid
				{
					vertexIndices.push_back(index);
				}
				else // Not valid
				{
					std::cerr << "Warning: Vertex indices must be positive integers. Relative (or "
						<< "zero) indices are not supported (line " << lineNo << ")." << std::endl;
				}
			}

			if (vertexIndices.size() == 3)
			{
				Tri triangle;
				triangle.materialId = activeMaterialId;
				triangle.vertices[0] = vertices[vertexIndices[0]];
				triangle.vertices[1] = vertices[vertexIndices[1]];
				triangle.vertices[2] = vertices[vertexIndices[2]];

				Object& activeObject = objects.back();
				activeObject.triangles.push_back(triangle);
			}
			else
			{
				std::cerr << "Warning: Face must contain exactly three valid vertex "
					<< "indices for voxelisation (line " << lineNo << ")." << std::endl;
			}
		}
	}

	if (tooManyComponentsLineNo)
	{
		std::cerr << "Warning: The 'v' element on line " << tooManyComponentsLineNo
			<< " contains more than three components. These components (possibly"
			<< " 'w' or a vertex colour?) will be ignored. Later lines may"
			<< " have the same problem (this has not been checked)." << std::endl;
	}

	return objects;
}
