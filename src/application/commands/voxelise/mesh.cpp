#include "mesh.h"

#include "base/logging.h"

#include <fstream>
#include <sstream>

using namespace Cubiquity;

using namespace std;

// https://en.wikipedia.org/wiki/Wavefront_.obj_file#Material_template_library
// and (official spec?) http://paulbourke.net/dataformats/mtl/
std::map<std::string, Col> loadMtlFile(const std::filesystem::path& pathToMat)
{
	ifstream file(pathToMat);
	std::map<std::string, Col> materials;
	std::map<std::string, Col>::iterator currentMaterial;
	bool foundDiffuse = true;

	auto finishMaterial = [&]()
	{ 
		if (!foundDiffuse)
		{
			log(Warning, "No diffuse colour found for material \'", currentMaterial->first, "\'.");
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
			Col colour({1.0f, 1.0f, 1.0f});
			foundDiffuse = false;
			auto result = materials.insert(std::make_pair(materialName, colour));
			if (!result.second)
			{
				log(Warning, "Material \'", materialName, "\' is already defined in \'", pathToMat, "\'.");
			}
			// Update even if the material couldn't be inserted. Pointing at the duplicate
			// is probably better than whatever we were pointing at before (if anything).
			currentMaterial = result.first;
		}

		if (element == "Kd")
		{
			std::array<float, 3> diffuse;
			iss >> diffuse[0] >> diffuse[1]>> diffuse[2];
			currentMaterial->second = diffuse;
			foundDiffuse = true;
		}

		if (iss.fail())
		{
			log(Error, "Failed to parse line ", lineNo, ". (\"", line, "\")");
		}
	}

	finishMaterial();

	if (materials.empty())
	{
		log(Warning, "No materials found in \'", pathToMat, "\'");
	}

	return materials;
}

// See https://en.wikipedia.org/wiki/Wavefront_.obj_file 
// and https://www.fileformat.info/format/wavefrontobj/egff.htm
// and (official spec?) http://paulbourke.net/dataformats/obj/ 
std::list<Object> loadObjFile(const std::filesystem::path& pathToObj)
{
	ifstream file(pathToObj);

	std::list<Object> objects;
	std::vector<Vertex> vertices;
	std::map<std::string, Col> materials;

	// Insert an (unused) dummy vertex because obj file indices start at 1.
	vertices.push_back({ 0,0,0 }); 

	// I have seen an OBJ file start without an 'o' or a 'usemtl' (or 
	// with them in the wrong order). Create a dummy object to hold these.
	Col activeMaterial({1.0f, 1.0f, 1.0f});
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
			std::filesystem::path pathToMat = pathToObj.parent_path().append(materialFilename);
			if (!std::filesystem::exists(pathToMat))
			{
				log(Warning, "Material file \'", pathToMat, "\' does not exist!");
			}
			std::map < std::string, Col > localMaterials = loadMtlFile(pathToMat);

			for (const auto& entry : localMaterials)
			{
				if (!materials.insert(entry).second)
				{
					log(Warning, "Material \'", entry.first, "\' was already defined in another material file.");
				}
			}
		}

		if (element == "usemtl")
		{
			std::string materialName;;
			iss >> materialName;

			auto materialIter = materials.find(materialName);
			activeMaterial = materialIter != materials.end() ?
				materialIter->second : Col({ 1.0f, 1.0f, 1.0f });
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
					log(Warning, "Vertex indices must be positive integers. Relative "
						"(or zero) indices are not supported (line ", lineNo, ").");
				}
			}

			if (vertexIndices.size() == 3)
			{
				Tri triangle;
				triangle.colour = activeMaterial;
				triangle.vertices[0] = vertices[vertexIndices[0]];
				triangle.vertices[1] = vertices[vertexIndices[1]];
				triangle.vertices[2] = vertices[vertexIndices[2]];

				Object& activeObject = objects.back();
				activeObject.triangles.push_back(triangle);
			}
			else
			{
				log(Warning, "Face must contain exactly three valid vertex "
					"indices for voxelisation (line ", lineNo, ").");
			}
		}
	}

	if (tooManyComponentsLineNo)
	{
		log(Warning, "The 'v' element on line ", tooManyComponentsLineNo, " contains more than "
			"three components. These components (possibly 'w' or a vertex colour?) will be ignored. "
			"Later lines may have the same problem (this has not been checked).");
	}

	return objects;
}
