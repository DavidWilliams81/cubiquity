#ifndef CUBIQUITY_UTILITY_MESH_H
#define CUBIQUITY_UTILITY_MESH_H

#include "base.h"

#include "voxelization.h"

#include <array>
#include <list>
#include <string>

typedef std::array<float, 3> Vertex;
struct Tri // Will rename this to 'Triangle' once we can avoid name conflicts with main library.
{
	std::array<Vertex, 3> vertices;
	Cubiquity::MaterialId materialId;
};

struct Object
{
	std::string name;
	std::vector<Tri> triangles;
};

// We return a list of objects (not a map) as the order can affect the voxelisation result.
std::list<Object> loadObjFile(const std::string& path, const std::string& filename);

#endif // CUBIQUITY_UTILITY_MESH_H