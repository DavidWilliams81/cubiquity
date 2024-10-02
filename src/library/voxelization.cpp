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
#include "voxelization.h"

#include "geometry.h"
#include "storage.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <execution>

#ifdef _MSC_VER // Not supported on Debian 12 GCC
#include <format>
#endif //_MSC_VER

#include <fstream>
#include <mutex>
#include <set>
#include <thread>

namespace Cubiquity
{
using namespace Internals;

// These compile-time controls are only for validation purposes
const bool evalWindingNumberHierarchy = true; // Evaluate hierarchy rather than flat list (faster)
const bool classifyNodesNotVoxels = true; // Classify octree nodes rather than voxels (faster)

/***************************************************************************************************
*    ______       _                         __  _   _  _    _  _                                   *
*    |  _  \     | |                       / / | | | || |  (_)| |                                  *
*    | | | | ___ | |__   _   _   __ _     / /  | | | || |_  _ | | ___                              *
*    | | | |/ _ \| '_ \ | | | | / _` |   / /   | | | || __|| || |/ __|                             *
*    | |/ /|  __/| |_) || |_| || (_| |  / /    | |_| || |_ | || |\__ \                             *
*    |___/  \___||_.__/  \__,_| \__, | /_/      \___/  \__||_||_||___/                             *
*                                __/ |                                                             *
*                               |___/                                                              *
***************************************************************************************************/

#ifdef _MSC_VER // Not supported on Debian 12 GCC

// Write out the winding number tree to a stream for inspection
void writePatch(const Patch& patch, std::ostream& file, int indent = 0)
{
	file << std::format("{}Patch {}: Triangle count = {}, Closing triangle count = {}",
		std::string(indent, '\t'), reinterpret_cast<uintptr_t>(&patch),
		patch.triangles.size(), patch.closingTriangles.size()) << std::endl;

	for (const auto& child : patch.children) {
		writePatch(child, file, indent + 1);
	}
}	

// Support function for exporting Collada
void exportGeometry(std::ostream& file, const Patch& patch)
{
	// For visualisation purposes it is most useful to only output geometry for leaf nodes
	if (patch.children.empty()) // Leaf node, output triangles
	{
		const int triangle_count = patch.triangles.size();
		const std::string name = std::to_string(reinterpret_cast<uintptr_t>(&patch));		

		file << std::format(
			"\t\t<geometry id=\"{}-mesh\" name=\"{}\">\n"
			"\t\t\t<mesh>\n"
			"\t\t\t\t<source id=\"{}-mesh-positions\">\n"
			"\t\t\t\t\t<float_array id=\"{}-mesh-positions-array\" count=\"{}\">",
			name, name, name, name, triangle_count * 3 * 3);

		for (const auto& triangle : patch.triangles) {
			for (const auto& vertex : triangle.vertices) {
				file << vertex[0] << " " << vertex[1] << " " << vertex[2] << " ";
			}
		}

		file << std::format("</float_array>\n"
			"\t\t\t\t\t<technique_common>\n"
			"\t\t\t\t\t\t<accessor source=\"#{}-mesh-positions-array\" count=\"{}\" stride=\"3\">\n"
			"\t\t\t\t\t\t\t<param name=\"X\" type=\"float\"/>\n"
			"\t\t\t\t\t\t\t<param name=\"Y\" type=\"float\"/>\n"
			"\t\t\t\t\t\t\t<param name=\"Z\" type=\"float\"/>\n"
			"\t\t\t\t\t\t</accessor>\n"
			"\t\t\t\t\t</technique_common>\n"
			"\t\t\t\t</source>\n"
			"\t\t\t\t<vertices id=\"{}-mesh-vertices\">\n"
			"\t\t\t\t\t<input semantic=\"POSITION\" source=\"#{}-mesh-positions\"/>\n"
			"\t\t\t\t</vertices>\n"
			"\t\t\t\t<triangles count=\"{}\">\n"
			"\t\t\t\t\t<input semantic=\"VERTEX\" source=\"#{}-mesh-vertices\" offset=\"0\"/>\n",
			name, triangle_count * 3, name, name, triangle_count, name);

		// Generate indices for our non-indexed triangles (as Collada seems to require them).
		file << "\t\t\t\t\t<p>";
		for (int index = 0; index < triangle_count * 3; index++) {
			file << index << " ";
		}

		// Closing tags
		file << "</p>\n\t\t\t\t</triangles>\n\t\t\t</mesh>\n\t\t</geometry>\n";
	}
	else // Internal node, process children
	{
		for (const auto& child : patch.children) {
			exportGeometry(file, child);
		}
	}
}

// Support function for exporting Collada
void exportNode(std::ostream& file, const Patch& patch, int indent)
{
	const std::string name = std::to_string(reinterpret_cast<uintptr_t>(&patch));
	file << std::format("{}<node id=\"{}-node\" name=\"{}-node\" type=\"NODE\">\n",
		std::string(indent+1, '\t'), name, name);
	
	// For visualisation purposes it is most useful to only output geometry for leaf nodes
	if (patch.children.empty())	{ // Leaf node, output geometry
		file << std::format("{}<instance_geometry url=\"#{}-mesh\" name=\"{}-name\"/>\n",
			std::string(indent+2, '\t'), name, name);
	}
	else { // Internal node, process children
		for (const auto& child : patch.children) {
			exportNode(file, child, indent + 1);
		}
	}

	file << std::format("{}</node>\n", std::string(indent+1, '\t'));
}

// Write the winding number tree as Collada for visualisation purposes.
// Collada is used because it is easy to write, and supports hierarchical
// transformations which are imported into Blender correctly.
void exportCollada(const std::string& filename, const Patch& patch)
{
	std::ofstream file(filename);

	file << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
	file << "<COLLADA xmlns = \"http://www.collada.org/2005/11/COLLADASchema\" "; // No newline
	file << "version = \"1.4.1\" xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\">\n";
	file << "\t<library_geometries>\n";

	exportGeometry(file, patch);

	file << "\t</library_geometries>\n";
	file << "\t<library_visual_scenes>\n";
	file << "\t\t<visual_scene id=\"Scene\" name=\"Scene\">\n";

	exportNode(file, patch, 2);

	file << "\t\t</visual_scene>\n";
	file << "\t</library_visual_scenes>\n";
	file << "\t<scene>\n";
	file << "\t\t<instance_visual_scene url=\"#Scene\"/>\n";
	file << "\t</scene>\n";
	file << "</COLLADA>\n";
	file.close();
}

#endif // _MSC_VER

/***************************************************************************************************
*   _    _ _           _ _               _   _                 _                                   *
*  | |  | (_)         | (_)             | \ | |               | |                                  *
*  | |  | |_ _ __   __| |_ _ __   __ _  |  \| |_   _ _ __ ___ | |__   ___ _ __ ___                 *
*  | |/\| | | '_ \ / _` | | '_ \ / _` | | . ` | | | | '_ ` _ \| '_ \ / _ \ '__/ __|                *
*  \  /\  / | | | | (_| | | | | | (_| | | |\  | |_| | | | | | | |_) |  __/ |  \__ \                *
*   \/  \/|_|_| |_|\__,_|_|_| |_|\__, | \_| \_/\__,_|_| |_| |_|_.__/ \___|_|  |___/                *
*                                 __/ |                                                            *
*                                |___/                                                             *
****************************************************************************************************
* Implements the hierarchical method of winding number evaluation and segmentation as described in *
* "Robust Inside-Outside Segmentation using Generalized Winding Numbers" by Jacobson et al (2013)  *
***************************************************************************************************/

float computeWindingNumber(const Vector3f& queryPoint, ConstTriangleSpan triangles)
{
	float windingNumber = 0.0f;

	for (const Triangle& triangle : triangles)
	{
		const auto& a = triangle.vertices[0];
		const auto& b = triangle.vertices[1];
		const auto& c = triangle.vertices[2];

		Vector3f qa = a - queryPoint;
		Vector3f qb = b - queryPoint;
		Vector3f qc = c - queryPoint;

		const float alength = length(qa);
		const float blength = length(qb);
		const float clength = length(qc);

		// Avoid potential NaNs/Infs in the result. If any of these values are
		// zero then the triangle does not contribute to the winding number.
		if (alength != 0 && blength != 0 && clength != 0)
		{
			// Normalize the vectors
			qa /= alength;
			qb /= blength;
			qc /= clength;

			// Subtracting qa from qb and qc is not strictly required,
			// but gives a more stable result if the triangle is far away.
			const float numerator = dot(qa, cross(qb - qa, qc - qa));
			const float denominator = 1.0f + dot(qa, qb) + dot(qa, qc) + dot(qb, qc);

			// Numerator of zero means we are on the surface (treat as no solid angle).
			if (numerator != 0)
			{
				//assert(std::isnormal(numerator));
				//assert(std::isnormal(denominator));
				windingNumber += 2.0f * ::atan2(numerator, denominator);
			}
		}
	}

	// Normalise to [-1.0f, 1.0f] range.
	windingNumber /= 4.0 * 3.14159265358979323846;

	// Make sure we got a 'valid' result.
	assert(std::isnormal(windingNumber) || (windingNumber == 0));

	return windingNumber;
}

float computeWindingNumber(const Vector3f& queryPoint, const Patch& patch)
{
	// Implementation of Algorithm 2 from Jacobson et al.
	if (patch.children.empty()) // Leaf node
	{
		return computeWindingNumber(queryPoint, patch.triangles);
	}
	else if (!patch.bounds.contains(queryPoint)) // Outside node
	{
		return -computeWindingNumber(queryPoint, patch.closingTriangles); // Note negation
	}
	else // Inside non-leaf node
	{
		float windingNumber = 0.0f;
		for (const auto& child : patch.children) {
			windingNumber += computeWindingNumber(queryPoint, child);
		}
		return windingNumber;
	}
}

TriangleList computeClosingTriangles(TriangleSpan triangles)
{
	typedef std::pair<Vector3f, Vector3f> TriangleEdge;
	std::unordered_map<TriangleEdge, int32_t, Internals::MurmurHash3<TriangleEdge> > edgeCounts;

	edgeCounts.max_load_factor(1.0); // Default anyway, changing it doesn't seem to help?
	edgeCounts.reserve(triangles.size() * 3);

	// Find the edge counts, where a count of non-zero indicates an exterior edge. I think the
	// algorithm handles degenerate triangles gracefully (at least they don't break the result).
	for (const auto& triangle : triangles)
	{
		for (int i = 0; i < 3; i++)
		{
			const Vector3f& v0 = triangle.vertices[i];
			const Vector3f& v1 = triangle.vertices[(i + 1) % 3];

			// Curiously, this condition seems arbitrary (same overall resuilt if we change
			// 'true'  to 'false') and just serves to order vertices so we can find pairs of
			// identical edges running in opposite directions. Jacobson et al state 'we
			// increment count(i, j) if i < j' but don't explain what 'i' and 'j' are. Are
			// they vertex positions or indices? I don't actually think it matters.
			if ((v0 < v1) == true) {
				edgeCounts[TriangleEdge(v0, v1)]++;
			} else {
				edgeCounts[TriangleEdge(v1, v0)]--;
			}
		}
	}

	// Generate closing triangles for exterior edges. I don't believe this solution is 
	// optimal (e.g. if a single triangle is missing then three closing triangles will
	// be generated) but it is simple and the approach described by Jakobson et al.
	TriangleList closingTriangles;
	for (const auto& edgeCount : edgeCounts)
	{
		TriangleEdge edge = edgeCount.first;
		int32_t count = edgeCount.second;

		if (count < 0)
		{
			std::swap(edge.first, edge.second);
			count = -count;
		}

		// Insert 'count' copies of the closing triangle, where count can be zero
		// (the most common case) or greater than one (as noted in the paper).
		const Vector3f& arbitraryVertex = triangles[0].vertices[0]; // Inside bounds
		// Jacobson et al state:
		//     "Finally all edges with count(i, j) != 0 are declared exterior and 
		//     triangulated with some arbitrary vertex k with orientation { i, j, k }
		//     if count(i, j) = c > 0 and {j, i, k} if count(i, j) = -c < 0".
		// But I think this is the wrong way around - if there are too many edges from
		// i->j then the closing triangles need to go in the *other* direction (j->i),
		// which is what we do here and it seems to work correctly.
		const Triangle closingTriangle(edge.second, edge.first, arbitraryVertex);
		closingTriangles.insert(closingTriangles.end(), count, closingTriangle);
	}

	return closingTriangles;
}

// Classify a point as being inside or outside the surface. Jacobson et al state:
// 
//   "In character meshes and CAD models, there may be [...] nearly duplicated patches of the
//    input mesh. These shift the winding number range locally [...]. This disqualifies simply
//    thresholding the winding number for final segmentation, hence our use of a carefully
//    crafted graphcut energy."
//
// However, elsewhere (https://github.com/libigl/libigl/issues/772) the same author says:
//
//   "In my experience, graphcut is often not necessary. It's good for handling some gnarly cases,
//    but I rarely use it."
//
// We use thresholding here as it is a simple and local operation which has so far been sufficient.
bool isInside(const Vector3f& queryPoint, const Mesh& mesh)
{
	float winding_number = evalWindingNumberHierarchy ?
		computeWindingNumber(queryPoint, *mesh.rootPatch) :
		computeWindingNumber(queryPoint, mesh.triangles);

	// Theoretically the threshold should be 0.5f. However, triangles often pass *exactly* through 
	// the centre of a voxel due to e.g. being modelled on a grid in a DCC package, and/or being
	// scaled by interger amounts. Rounding errors mean that such triangles can end up with a noisy 
	// surface, and the epsilon below helps to avoid this.
	return std::abs(winding_number) > 0.5f + 0.001f;
}

Patch::Patch(TriangleSpan triSpan)
	:triangles(triSpan)
{
	bounds = computeBounds(triangles);
	closingTriangles = computeClosingTriangles(triangles);

	// Decide whether the patch should be split further. There are more sophisticated
	// criteria which could be used (e.g. consider volume of bounding boxes) but this
	// simple approach from Jacobson et al seems to give good results.
	const int minTriangles = 100; // Soft limit, leaf tri count can be less
	if (triangles.size() > closingTriangles.size() && triangles.size() > minTriangles)
	{
		// Simply splitting along the longest axis works surprisingly well for most cases.
		// I have tried more advanced approaches such as using mesh connectivity information
		// or evaluating multiple split points but the limited benefits were not worth the
		// additional code complexity and increase in tree construction time. There are
		// still more things I could try (more than two groups, non-axis-aligned split
		// planes, clustering algorithms, etc), but I'm not convinced they're worthwhile.
		Vector3f dims = bounds.upper() - bounds.lower();
		int longestAxis = std::max_element(dims.begin(), dims.end()) - dims.begin();

		// Sort all triangles in this patch (other triangles in the mesh are unchanged).
		std::sort(triangles.begin(), triangles.end(), 
				    [=](const Triangle& a, const Triangle& b)
		{
			return a.centre()[longestAxis] < b.centre()[longestAxis];
		});

		// Split at median. Mean can be better but median is simpler and more robust
		// against degenerate inputs (all triangles cannot fall on same side of split).
		auto split_point = triangles.begin() + (triangles.size() / 2);
		children.emplace_back(TriangleSpan({ triangles.begin(), split_point }));
		children.emplace_back(TriangleSpan({ split_point,   triangles.end() }));
	}
}

/***************************************************************************************************
*   ___________   _____                   _____                               _                    *
*  |____ |  _  \ /  ___|                 /  __ \                             (_)                   *
*      / / | | | \ `--.  ___ __ _ _ __   | /  \/ ___  _ ____   _____ _ __ ___ _  ___  _ __         *
*      \ \ | | |  `--. \/ __/ _` | '_ \  | |    / _ \| '_ \ \ / / _ \ '__/ __| |/ _ \| '_ \        *
*  .___/ / |/ /  /\__/ / (_| (_| | | | | | \__/\ (_) | | | \ V /  __/ |  \__ \ | (_) | | | |       *
*  \____/|___/   \____/ \___\__,_|_| |_|  \____/\___/|_| |_|\_/ \___|_|  |___/_|\___/|_| |_|       *
*                                                                                                  *
****************************************************************************************************
* This functionality scan-converts a 3D triangle into a corresponding set of '3D fragments' (voxel *
* positions). For large triangles it does this by recursively subdividing them into smaller        *
* triangles. For small triangles it iterates over each voxel in the bounding box and tests whether *
* the distance to the trangle is below a threshold, or alternatively whether the triangle touches  *
* the voxel's 'intersection target' described by 'A Topological Approach to Voxelization' (Laine). *
***************************************************************************************************/

// A callback for each voxel which is part of the triangle. 
typedef std::function<void(int32, int32, int32, MaterialId)> DrawVoxelFunc;

void drawSmallTriangle(const Triangle& triangle, MaterialId matId,
	                   float thickness, DrawVoxelFunc drawVoxel)
{
	Box3f triBounds = computeBounds(triangle.vertices);

	// Expand bounds to encompass all candidate voxels.
	if (thickness >= 0.0f) {
		triBounds.dilate(thickness);
	} else {
		triBounds.dilate(0.5f); // Extent of intersection target.
	}

	// Shrink to fit integer bounds (captures all integer positions within float bounds).
	// Note: We could probably make a Box member function for this.
	Box3i triBoundsAsInt = {
		static_cast<Vector3i>(ceil(triBounds.mExtents[0])),
		static_cast<Vector3i>(floor(triBounds.mExtents[1])) };

	for (int32_t z = triBoundsAsInt.lower().z(); z <= triBoundsAsInt.upper().z(); z++)
	{
		for (int32_t y = triBoundsAsInt.lower().y(); y <= triBoundsAsInt.upper().y(); y++)
		{
			for (int32_t x = triBoundsAsInt.lower().x(); x <= triBoundsAsInt.upper().x(); x++)
			{
				const Vector3f pos = { x,y,z };

				if (thickness >= 0.0f) // Draw surface with user-provided thickness.
				{
					// Only voxels which are within required distance of surface. 
					if ((distance(pos, triangle) <= thickness))
					{
						drawVoxel(x, y, z, matId);
					}
				}
				else // Juse a draw a minimal (6-seperating) surface.
				{
					// Test the triangle against the voxel's intersection target. See the paper
					// 'A Topological Approach to Voxelization' (Laine) for more details (note
					// that the slides are easier to read). Possible improvements:
					//   * Use single (diagonal) intersection target (probably more complex).
					//   * Each intersection target is axis-aligned, so we can flatten the triangle
					//     along that axis and do a simpler 'point in 2D triangle' test.
					//   * Pull normal calculation out of ray-triangle intersection test and
					//     compute at higher level (might be more accurate too?)
					for (int axis = 0; axis < 3; axis++)
					{
						Ray3f intersectionTarget = { pos, { 0.0f, 0.0f, 0.0f } };
						intersectionTarget.mOrigin[axis] -= 0.5f;
						intersectionTarget.mDir[axis] = 1.0f;

						float t = 0.0f;
						bool hit = intersect(intersectionTarget, triangle, t);
						if (hit && t <= 1.0f)
						{
							drawVoxel(x, y, z, matId);
							break;
						}
					}
				}
			}
		}
	}
}

void drawLargeTriangle(const Triangle& triangle, MaterialId matId,
	                   float thickness, DrawVoxelFunc drawVoxel)
{
	// Find the longest side
	int longestSide = 0;
	if (triangle.sideLength(1) > triangle.sideLength(0)) { longestSide = 1; }
	if (triangle.sideLength(2) > triangle.sideLength(1)) { longestSide = 2; }

	// Split large triangle in half as bounds may overlap many redundant voxels.
	const float maxSideLength = 16.0f;
	if (triangle.sideLength(longestSide) > maxSideLength) {
		Vector3f midpoint = (triangle.vertices[longestSide] +
			triangle.vertices[(longestSide + 1) % 3]) / 2.0f;
		// Build two small triangles from input triangle
		for (int i = 0; i < 2; i++) {
			Triangle halfTri = triangle; // Copy input triangle
			halfTri.vertices[(longestSide + i) % 3] = midpoint; // Shift one vertex to midpoint
			drawLargeTriangle(halfTri, matId, thickness, drawVoxel);
		}
	}
	else {
		// Draw small triangles directly.
		drawSmallTriangle(triangle, matId, thickness, drawVoxel);
	}
}

void drawTriangles(const TriangleList& triList, const MaterialList& materialList, float thickness, DrawVoxelFunc drawVoxel)
{
	uint32_t trianglesDrawn = 0;
	for (uint i = 0; i < triList.size(); i++)
	{
		drawLargeTriangle(triList[i], materialList[i], thickness, drawVoxel);
		reportProgress("Performing 3D scan conversion", 1, ++trianglesDrawn, triList.size());
	}
}

/***************************************************************************************************
*     _   _           _        _____           _             _   _                                 *
*    | \ | |         | |      |  ___|         | |           | | (_)                                *
*    |  \| | ___   __| | ___  | |____   ____ _| |_   _  __ _| |_ _  ___  _ __                      *
*    | . ` |/ _ \ / _` |/ _ \ |  __\ \ / / _` | | | | |/ _` | __| |/ _ \| '_ \                     *
*    | |\  | (_) | (_| |  __/ | |___\ V / (_| | | |_| | (_| | |_| | (_) | | | |                    *
*    \_| \_/\___/ \__,_|\___| \____/ \_/ \__,_|_|\__,_|\__,_|\__|_|\___/|_| |_|                    *
*                                                                                                  *
****************************************************************************************************
* Evaluate the winding number for each node in the octree to decide whether it is inside the mesh. *
***************************************************************************************************/

struct NodeToTest
{
	uint32_t index;
	uint32_t childId;
	Vector3f centre;
	bool result;
};

class NodeFinder
{
public:
	// Not for general boxes - only works for node bounds which are cubic,
	// powers-of-two sizes, aligned to power-of-two boundaries, etc.
	Box3i childBounds(Box3i nodeBounds, uint childId)
	{
		uint childX = (childId >> 0) & 0x01;
		uint childY = (childId >> 1) & 0x01;
		uint childZ = (childId >> 2) & 0x01;
		Vector3i childOffset({ static_cast<int>(childX), static_cast<int>(childY), static_cast<int>(childZ) }); // childOffset holds zeros or ones.

		// Careful ordering of operations to avoid signed integer overflow. Note that child
		// node dimensions might max-out the signed integer type but should not exceed it.
		const Vector3i childNodeDimsInCells = ((nodeBounds.upper() - Vector3i({ 1, 1, 1 })) / 2) - (nodeBounds.lower() / 2);
		Vector3i childLowerBound = nodeBounds.lower() + (childNodeDimsInCells * childOffset) + childOffset;
		Vector3i childUpperBound = childLowerBound + childNodeDimsInCells;
		return Box3i(childLowerBound, childUpperBound);
	}

	bool operator()(NodeDAG& nodes, uint32 nodeIndex, const Box3i& nodeBounds)
	{
		// Signal to stop traversing parts of the tree which do not overlap our voxelised object.
		if (!overlaps(nodeBounds, mBounds)) { return false; }

		// If this is a non-leaf (internal) node then check if any of its children are leaves.
		// If so, add them to the list of nodes on which the inside/outside test will be performed.
		// FIXME - Think about whether we really need to do it like this. This is currently the
		// only node visitor which uses the 'nodes' input, so if we can avoid it here then maybe
		// we can simplify the interface?
		if (!isMaterialNode(nodeIndex))
		{
			for (unsigned int childId = 0; childId < 8; childId++)
			{
				uint32 childNodeIndex = nodes[nodeIndex][childId];
				if (isMaterialNode(childNodeIndex))
				{
					Box3i childNodeBounds = childBounds(nodeBounds, childId);
					if (overlaps(childNodeBounds, mBounds))
					{
						NodeToTest toTest;
						toTest.index = nodeIndex;
						toTest.childId = childId;
						toTest.centre = static_cast<Vector3f>((childNodeBounds.lower() + childNodeBounds.upper())) * 0.5f;
						mNodes.push_back(toTest);
					}
				}
			}
		}

		return true;
	}

	std::vector<NodeToTest> nodes()
	{
		return mNodes;
	}

private:
	std::vector<NodeToTest> mNodes;

public:
	Box3i mBounds;

};

std::vector<NodeToTest> findNodes(Volume& volume, Box3i bounds)
{
	NodeFinder nodeFinder;
	nodeFinder.mBounds = bounds;
	visitVolumeNodes(volume, nodeFinder);
	return nodeFinder.nodes();
}

void classifyNodes(std::vector<NodeToTest>& nodesToTest, NodeDAG& /*nodeData*/, Mesh& mesh)
{
	int i = 0;
	std::stringstream ss;
	ss << "Classifying " << nodesToTest.size() << " nodes";
	std::mutex m;

	std::for_each(std::execution::par, nodesToTest.begin(), nodesToTest.end(), [&](auto&& nodeToTest)
	{
		nodeToTest.result = isInside(nodeToTest.centre, mesh);
		std::lock_guard<std::mutex> guard(m);

		reportProgress(ss.str().c_str(), 0, i++, nodesToTest.size()-1);
	});
}

uint32 prune(NodeDAG& nodes, uint32 nodeIndex)
{
	Node node = nodes[nodeIndex];
	for (int i = 0; i < 8; i++)
	{
		uint32 nodeChildIndex = node[i];
		if (!isMaterialNode(nodeChildIndex))
		{
			node[i] = prune(nodes, nodeChildIndex);
		}
	}

	return nodes.isPrunable(node) ? node[0] : nodes.insert(node);
}

void prune(Volume& volume)
{
	uint32 rootNodeIndex = getRootNodeIndex(volume);
	NodeDAG& nodes = getNodes(volume);
	uint32 newRootNodeIndex = prune(nodes, rootNodeIndex);
	volume.setRootNodeIndex(newRootNodeIndex);
}

/***************************************************************************************************
*    ___  ___          _       _   _               _ _           _   _                             *
*    |  \/  |         | |     | | | |             | (_)         | | (_)                            *
*    | .  . | ___  ___| |__   | | | | _____  _____| |_ ___  __ _| |_ _  ___  _ __                  *
*    | |\/| |/ _ \/ __| '_ \  | | | |/ _ \ \/ / _ \ | / __|/ _` | __| |/ _ \| '_ \                 *
*    | |  | |  __/\__ \ | | | \ \_/ / (_) >  <  __/ | \__ \ (_| | |_| | (_) | | | |                *
*    \_|  |_/\___||___/_| |_|  \___/ \___/_/\_\___|_|_|___/\__,_|\__|_|\___/|_| |_|                *
*                                                                                                  *
***************************************************************************************************/

// A fast approach to voxelisation which evaluates the winding number for each octree node
// and uses multiple threads. This is the method which is used under normal circumstances.
void doPerNodeVoxelisation(Volume& volume, Mesh& mesh, MaterialId fill, MaterialId background)
{
	Box3i voxelisationBounds = computeBounds(volume, background);

	// Find all the nodes - both the single-voxel nodes which are part of the
	// voxelised surface and (hopefully larger) nodes which are either side of it.
	auto nodesToTest = findNodes(volume, voxelisationBounds);

	// Classify all nodes according to which side of the surface they are on.
	NodeDAG& nodeData = Internals::getNodes(volume);
	classifyNodes(nodesToTest, nodeData, mesh);

	// Apply the results to the volume
	for (auto& toTest : nodesToTest)
	{
		MaterialId resultingMaterial = toTest.result ? fill : background;
		nodeData.nodes().setNodeChild(toTest.index, toTest.childId, resultingMaterial);
	}

	// The voxelisation process can cause the volume to become unpruned, which we consider to be an invalid state.
	// This might be because the shell voxelisaion is too thick, though I think we have avoided that. But even so,
	// the mesh might represent a small box touching eight voxels, all of which are inside. These would be
	// individually classified and would all get the same value, so they should be pruned and replaced by the parent.
	prune(volume);
}

// A very slow approach to voxelisation which evaluates the winding number for every voxel It also
// only runs on a single thread. It is only intended for validation of the hierarchical version.
void doPerVoxelVoxelisation(Volume& volume, Mesh& mesh, MaterialId fill, MaterialId background)
{
	Box3i voxelisationBounds = static_cast<Box3i>(mesh.bounds);
	voxelisationBounds.dilate(2); // Make sure we really cover everything

	Vector3i minBound = voxelisationBounds.lower();
	Vector3i maxBound = voxelisationBounds.upper();

	// Iterate over each voxel within the bounds and classify as inside or outside
	for (int32 volZ = minBound.z(); volZ <= maxBound.z(); volZ++)
	{
		for (int32 volY = minBound.y(); volY <= maxBound.y(); volY++)
		{
			for (int32 volX = minBound.x(); volX <= maxBound.x(); volX++)
			{
				Vector3f queryPoint = { static_cast<float>(volX), static_cast<float>(volY), static_cast<float>(volZ) };
				auto material = isInside(queryPoint, mesh) ? fill : background;
				volume.setVoxel(volX, volY, volZ, material);
			}
		}
		reportProgress("Classifying voxels (Brute Force!)", minBound.z(), volZ, maxBound.z());
	}
}

void voxelize(Volume& volume, Mesh& mesh, MaterialId fill, MaterialId background)
{
	// TODO - Need to think how triangle colours and fill colour should be used if one, 
	// both, or neither are provided. What makes for the simplest and most useful API?

	if (mesh.isInsideOut) {
		log(INF, "Mesh is inside-out, swapping material overrides.");
		std::swap(fill, background);
	}

	volume.fill(background);

	if (mesh.isClosed && (fill != background)) { // Do a proper 'solid' (filled) voxelisation
		if (classifyNodesNotVoxels)
		{
			// When drawing the mesh into the volume we can use a 'checkerboard' pattern of materials.
			// This provides high-frequency detail which prevents the octree from being pruned. A thin
			// surface is high-frequency anyway, but if e.g. two surfaces come close together then a set of
			// eight voxels can all be set which would then be pruned. The checkerboard does not prevent DAG
			// deduplication, but that does not happen automatically anyway.
			drawTriangles(mesh.triangles, mesh.materials, -1.0,
				[&](int32 x, int32 y, int32 z, MaterialId matId) {
					MaterialId checkerboard = ((x & 0x1) ^ (y & 0x1) ^ (z & 0x1));
					volume.setVoxel(x, y, z, background + checkerboard + 1);
				});

			doPerNodeVoxelisation(volume, mesh, fill, background);
		} else {
			// This path is for debug and validation only.
			doPerVoxelVoxelisation(volume, mesh, fill, background);
		}

		// Write the surface voxels with their correct materials as specified in the mesh.
		// Only write those which were found to be inside the mesh, except for meshes marked
		// as thin (as these may pass between voxels and hence not enclose them, resulting in holes).
		// 
		// Note that the triangle order can matter when multiple triangles pass close to a voxel,
		// and by drawing them in the user-supplied order we let the user control the result.
		drawTriangles(mesh.triangles, mesh.materials, 1.0,
			[&](int32 x, int32 y, int32 z, MaterialId matId) {
				if (volume.voxel(x, y, z) != background || mesh.isThin) {
					volume.setVoxel(x, y, z, matId);
				}
			});
	} else { // Fall back to just drawing the shell of the object (in user-supplied order, as above)
		auto draw = [&](int32 x, int32 y, int32 z, MaterialId matId) {
			volume.setVoxel(x, y, z, matId);
			};
		drawTriangles(mesh.triangles, mesh.materials, -1.0, draw);
	}
}

void Mesh::addTriangle(const Triangle& tri, MaterialId matId)
{
	triangles.push_back(tri);
	materials.push_back(matId);
}

MaterialId findMainMaterial(const Mesh& mesh)
{
	// Find the main material of the surface as the material which covers the greatest area. 
	// Attempts to use winding numbers have less obvious behaviour for open/inverted surfaces.
	std::vector<float> areas(MaterialCount);
	std::fill(areas.begin(), areas.end(), 0.0f);
	for (uint i = 0; i < mesh.triangles.size(); i++)
	{
		areas[mesh.materials[i]] += mesh.triangles[i].area();
	}
	return std::distance(areas.begin(), std::max_element(areas.begin(), areas.end()));
}

void Mesh::build()
{
	bounds = computeBounds(triangles);

	// Sample the winding number at various points in and around the object to determine how
	// well-formed the mesh is. Note that this is not a perfect test (e.g. a doubled-up
	// hemisphere will look like valid geometry) but it catches a lot of real-world scenarios.		
	bool allValid = true;
	bool anyPositive = false;
	bool anyNegative = false;
	const int sampleCount = 100;
	const float tolerance = 0.1f;

	Cubiquity::Box3fSampler sampler(bounds);
	for (uint32_t i = 0; i < sampleCount; i++)
	{
		Vector3f point = sampler.next();

		const float windingNumber = computeWindingNumber(point, triangles);
		const float absWindingNumber = std::abs(windingNumber);

		if (absWindingNumber >= tolerance) // Small winding numbers don't tell us anyting.
		{
			anyPositive |= windingNumber > 0.0f;
			anyNegative |= windingNumber < 0.0f;

			const bool isSufficient = absWindingNumber > (1.0f - tolerance);
			const bool isExcessive = absWindingNumber > (1.0f + tolerance);

			if (!isSufficient || isExcessive) {
				log(INF, "Absolute winding number of ", absWindingNumber, " found at position ", point, ".");
				if (!isSufficient) {
					log(INF, "\t(This indicates the mesh is not closed or has inconsistant winding)");
					allValid = false;
					break;
				}
				else { // is excessive
					log(INF, "\t(This indicates the mesh has doubled-up triangles or separate surface details)");
				}
			}
		}
	}

	const bool mixedSigns = anyNegative && anyPositive;
	isClosed = allValid && (!mixedSigns);
	isInsideOut = isClosed && anyNegative; // Mesh can only be inside out if closed.

	if (isInsideOut) {
		log(INF, "Exclusively negative winding numbers indicate the mesh is inside-out.");
	}

	if (isClosed) {
		patchTriangles = triangles; // Copy gets sorted during patch construction
		rootPatch = std::make_optional<Patch>(patchTriangles);
		//std::ofstream file("tree.txt"); writePatch(*rootPatch, file);
		//exportCollada(std::string("collada_") + name + ".dae", *rootPatch);
	}
}

} // namespace Cubiquity
