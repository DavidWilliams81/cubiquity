#!/usr/bin/env python3

# This is a Python script (rather than a Bash script or batch file)
# simply to save me maintaining two platform-specific versions.
import subprocess
subprocess.check_call(["../build/cubiquity", "voxelise", "axis.obj", "--output=axis.vol", "--scale=32"])
subprocess.check_call(["../build/cubiquity", "voxelise", "building.obj", "--output=building.vol", "--scale=256"])
subprocess.check_call(["../build/cubiquity", "voxelise", "sci-fi.obj", "--output=sci-fi.vol", "--scale=64"])
subprocess.check_call(["../build/cubiquity", "voxelise", "shapes.obj", "--output=shapes.vol", "--scale=64"])
subprocess.check_call(["../build/cubiquity", "voxelise", "terrain.png", "--output=terrain.vol", "--scale=511"])
