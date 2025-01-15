Voxelization
============

Introduction
------------
The Cubiquity voxelizer converts a polygon mesh into voxels. Conceptually it does this by testing each voxel position to determine whether or not it is inside the mesh, though in practice a number of optimizations are applied.

If the input mesh is properly constructed then a true 'solid' voxelization can be performed (i.e. the output is filled, it is not just a hollow shell). There is also support for scenes containing multiple meshes, thin meshes, 'inverted' meshes, and meshes using multiple materials. Note that texture-mapped meshes are not currently supported.

The voxelizer is designed for offline use as voxelization may take anywhere between a few seconds to several hours, depending on the complexity of the input mesh and the desired output resolution. The input is a Wavefront Object (.obj) with corresponding materials (.mt). The output is written as a sparse voxel DAG but can then be exported to other formats.

Quick Start
-----------
An example Wavefront Object file called 'shapes.obj' is provided with Cubiquity. This file was created in Blender using some primitive shapes and then exported to the required .obj format. See section XXXXX for more details on authoring .obj files with Blender.

The example input mesh can be voxelized with the following command:
```
cubiquity voxelize shapes.obj
```
This loads 'shapes.obj' (and its corresponding material file, shapes.mtl), performs the voxelization, and writes out the result as 'shapes.dag', an efficient sparse representation of the voxel data.

The output can then be viewed in Cubiquity as follows:
```
cubiquity view shapes.dag
```
Or exported to MagicaVoxel .vox format with:
```
cubiquity export vox shapes.dag
```
The exported files 'shapes.vox' can then be opened with MagicaVoxel or other compatible applications.

Voxelization controls
---------------------

--scale Scales the scene by the specified amount prior to voxelization. The choice of value will depend on the units of the input mesh and the desired size of the resulting voxel expressed in these units. For example, if the input mesh is in meters and the desired voxel size is 1cm, then a scale of 100 would be appropriate.

There is no default value as it is too application-specific. If scale is not specified, then a suitable scale factor is instead computed based on the desired size in voxels (see the --size parameter).


--size Specifies the size of the output in voxels as measured along the longest axis.

The size and scale parameters cannot be used together.

Preparing suitable meshes
-------------------------
Todo...
