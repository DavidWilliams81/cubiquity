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
#include <random>

namespace Cubiquity
{
	//! Respresents a set of triangles in a way which is optimised for winding number calculation.
	//! This class provides and implementation of the hierarchical approach to computing generalized winging numbers, as described in [REF]
	class ClosedTriangleTree : private Internals::NonCopyable // Avoid deep vs. shallow copy ambiguity
	{
	public:
		//! Constructs a ClosedTriangleTree from a TriangleList.
		explicit ClosedTriangleTree(const TriangleList &triangles, bool tidyInput = true);

	private:

		Box3f mBounds;
		TriangleList mTriangles;
		TriangleList mClosingTriangles;
		std::unique_ptr<ClosedTriangleTree> mChildren[2];

		// Making this function a friend lets us keep this class opaque to the user.
		friend float computeWindingNumber(const Vector3f& queryPoint, const ClosedTriangleTree& meshNode);
	};

	//! Calculates the winding number for a query point and a list of triangles
	//! This overload uses a simple brute-force approch and is privided mostly for reference.
	float computeWindingNumber(const Vector3f& queryPoint, const TriangleList& triangles);
	//! Calculates the winding number for a query point and a hierarchical triangle structure.
	//! This overload is much faster than using a triangle list and the theorectical result is
	//! exactly the same. In practice there are some small deviations for at least two reasons:
	//!     * Diffeent triangles are processed and in a different order, hence there is
	//!       potential for floating-point incacuracies to creep in.
	//!     * The triangles in the hierarchical struture are (usually) quantized slightly to 
	//!       improve floating-point comparisons (see ClosedTriangleTree).
	float computeWindingNumber(const Vector3f& queryPoint, const ClosedTriangleTree& meshNode);

	typedef std::list<std::pair<MaterialId, ClosedTriangleTree> > ClosedTriangleTreeList;

	enum class Thickness { Separate6, Separate26, Custom };

	void scanConvert3D(const Triangle& triangle, MaterialId matId, Volume& volume, Thickness thickness, float multiplier = 1.0f);
	void scanConvert3DRecursive(const Triangle& triangle, MaterialId matId, Volume& volume, Thickness thickness, float multiplier = 1.0f);

	Geometry mergeSubObjects(const Geometry& input);

	void voxelizeShell(Volume& volume, const Geometry& splitTriangles, bool preserveSurfaceMaterials, ProgressMonitor* progMon = nullptr);

	//! Perform a solid voxelisation of the given geometry into the volume.
	//!
	//! The provided geometry is assumed to be in 'volume space', i.e. it is correctly transformed
	//! to be in the desired location in the the volume. The algorithm works by firstly drawing the
	//! 'shell' of the geometry into the volume using a 3D version of scan conversion, and then 
	//! classifying each part of the volume (either per-voxel or per-octree-node) as inside or
	//! outside. It is generally this second part which dominates the runtime.
	//!
	//! The voxelisation process works best under the following conditions:
	//!
	//!     * All polygons in the scene have consistant winding. There should not be a mix of
	//!       clockwise and counter-clockwise wound polygons. It can cause significant problems if
	//!       this rule is violated.
	//!     * The scene should consist of closed objects ideally with a single material. Single-
	//!       sided objects represent windows or a ground plane can cause significant problem. 
	//!       These should be given some thickness if possible.
	//!     * There should be no hierarchy amounst the objects, as they will be treated as a flat
	//!       list.
	//!
	//! \param preserveSurfaceMaterials Ensures that materials which get written to the volume
	//!        during the scan conversion do not get overwritten during the classification. This can
	//!        be useful when one object contains multiple different materials, as otherwise the
	//!        material which gets chosen to fill the object might overwrite the surface. It is also
	//!        useful in conjunction wth 'internalMaterialOveride' to create objects with a different
	//!        material on the surface compared to the inside, or to represent hollow objects.
	//! \param internalMaterialOveride Causes all objects to be filled with the specified material
	//!        instead of their own. If 'preserveSurfaceMaterials' is disabled then this results in a
	//!        binary voxelisation.
	//! \param progMon An implementation of ProgressMonitor to monitor the voxelisation.
	//! \param useBruteForce Perform the classification per-voxel instead of per-octree-node. This is
	//!        much slower and for a well formed mesh the result should be the same, but it may be
	//!        useful on scenes with open meshes, ground planes, windows, etc.
	void voxelize(Volume& volume, Geometry& splitTriangles,	bool preserveSurfaceMaterials = false,
		uint16 internalMaterialOveride = 0, ProgressMonitor* progMon = nullptr, bool useBruteForce = false);
}

#endif //CUBIQUITY_VOXELIZATION_H
