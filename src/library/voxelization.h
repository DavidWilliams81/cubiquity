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
#ifndef CUBIQUITY_VOXELIZATION_H
#define CUBIQUITY_VOXELIZATION_H

#include "utility.h"
#include "geometry.h"
#include "storage.h"

#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <span>
#include <unordered_set>

namespace Cubiquity
{
	typedef std::vector<MaterialId> MaterialList;

	// Hierarchical representation of a list of triangles which is optimised
	// for computing generalized winding numbers, as described by Jacobson et al.
	class Patch : private Internals::NonCopyable
	{
	public:
		explicit Patch(TriangleSpan triSpan);

		Box3f bounds;
		TriangleSpan triangles;
		TriangleList closingTriangles;
		std::vector<Patch> children;
	};

	class Mesh
	{
	public:
		Mesh() {}

		// We should note in the API docs that the order of triangles can metter. When
		// multiple triangles write to the same voxel the last one in the list takes effect.
		void addTriangle(const Triangle& tri, MaterialId matId);
		void build();

		// General properties
		std::string name;
		Box3f bounds;

		// Geometry (stored in user-submitted order)
		TriangleList triangles;
		MaterialList materials;

		// Flags describing how well-formed the geomentry is
		bool isClosed = false;
		bool isInsideOut = false;

		// Thin meshes may pass between voxels
		bool isThin = false;

		// Hierarchical representation
		std::optional<Patch> rootPatch; // Root of hierachical representation
		TriangleList patchTriangles;    // Spatially-sorted copy of triangles
	};

	MaterialId findMainMaterial(const Mesh& mesh);

	// Voxelise the mesh into the volume
	void voxelize(Volume& volume, Mesh& mesh, MaterialId fill, MaterialId background);
}

#endif //CUBIQUITY_VOXELIZATION_H
