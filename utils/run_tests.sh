# This script should be run from the build folder as follows:
# ../utils/run_tests.sh 2>&1 | tee "run_tests_output_$(date +%Y%m%d_%H%M%S).txt"

#!/usr/bin/env bash
set -euo pipefail # exit on errors, unset vars, or failed pipes

echo "========================================================================="
echo " Build Information"
echo "========================================================================="

# Git info
BRANCH=$(git branch --show-current)
COMMIT=$(git rev-parse --short HEAD)
COMMIT_DATE=$(git show -s --format=%ci HEAD)
VERSION=$(git describe --tags --always 2>/dev/null)
BUILD_TIME=$(date +"%Y-%m-%d %H:%M:%S")

echo "Branch:        $BRANCH"
echo "Commit:        $COMMIT"
echo "Commit Date:   $COMMIT_DATE"
echo "Version:       $VERSION"
echo "Build Time:    $BUILD_TIME"
echo

echo "========================================================================="
echo " System Information"
echo "========================================================================="

# OS and kernel
if command -v lsb_release &>/dev/null; then
    echo "OS:            $(lsb_release -d | cut -f2)"
else
    echo "OS:            $(cat /etc/*release | grep PRETTY_NAME | cut -d= -f2 | tr -d '\"')"
fi

echo "Kernel:        $(uname -sr)"
echo "Architecture:  $(uname -m)"
echo "Hostname:      $(hostname)"
echo "Uptime:        $(uptime -p)"

# CPU info
echo "CPU:           $(lscpu | grep 'Model name' | sed 's/Model name:[ \t]*//')"
echo "Cores:         $(nproc)"

# Memory info
echo "Memory:        $(free -h | awk '/Mem:/ {print $2 " total, " $3 " used"}')"

# GPU info
printf "GPU:           " # No newline
if command -v nvidia-smi &>/dev/null; then
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
elif command -v rocm-smi &>/dev/null; then
    rocm-smi --showproductname --showdriverversion --showmeminfo vram | grep -E "Product|VRAM|Driver"
elif command -v glxinfo &>/dev/null; then
    glxinfo | grep "OpenGL renderer" | sed 's/^/OpenGL:        /'
elif command -v lspci &>/dev/null; then
    GPU=$(lspci | grep -E "VGA|3D|Display" | sed 's/^/GPU:           /')
    echo "$GPU"
else
    echo "GPU:           [No GPU tools detected]"
fi

echo "========================================================================="
echo

# Prepare tests
DATA=../data
DST=output # Destination

# Check if the directory already exists
if [ -d "$DST" ]; then
    echo "Directory '$DST' already exists."
    read -r -p "Do you want to delete and recreate it? [y/N] " reply
    reply=${reply,,}   # convert to lowercase

    if [[ "$reply" == "y" || "$reply" == "yes" ]]; then
        echo "Removing '$DST'..."
        rm -rf -- "$DST" || {
            echo "Error: failed to remove '$DST'." >&2
            exit 1
        }
    else
        echo "Keeping existing directory. Exiting."
        exit 0
    fi
fi

# (Re)create the directory
mkdir -p -- "$DST" || {
    echo "Error: failed to create '$DST'." >&2
    exit 1
}

echo "========================================================================="
echo " Running internal tests..."
echo "========================================================================="
echo

./cubiquity test

echo "========================================================================="
echo " Running volume generation tests..."
echo "========================================================================="
echo

# Note these use the Unix 'time' command and the Cubiquity '--quiet' switch. It
# Should be possible to remove both of these once we have a proper progress bar
# working with these generators.

echo "Generating fractal noise..."
echo "---------------------------"
time ./cubiquity generate fractal_noise --quiet --output=$DST/fractal_noise.dag
echo
md5sum $DST/fractal_noise.dag $DST/fractal_noise.toml
echo

echo "Generating Menger sponge..."
echo "---------------------------"
time ./cubiquity generate menger_sponge --quiet --size=243
mv output.dag $DST/menger_sponge.dag
mv output.toml $DST/menger_sponge.toml
echo
md5sum $DST/menger_sponge.dag $DST/menger_sponge.toml
echo

echo "Generating Worley noise..."
echo "--------------------------"
time ./cubiquity generate worley_noise --quiet --output=$DST/worley_noise.dag
echo
md5sum $DST/worley_noise.dag $DST/worley_noise.toml
echo

echo "========================================================================="
echo " Running volume import tests..."
echo "========================================================================="
echo

echo "Importing raw binary array..."
echo "-----------------------------"
md5sum $DATA/tests/fractal_worley_noise.bin $DATA/tests/fractal_worley_noise.txt
./cubiquity import bin $DATA/tests/fractal_worley_noise.bin --output=$DST/fractal_worley_noise.dag
echo
md5sum $DST/fractal_worley_noise.dag $DST/fractal_worley_noise.toml
echo

echo "========================================================================="
echo " Running volume export tests..."
echo "========================================================================="
echo

echo "Exporting raw binary array..."
echo "-----------------------------"
md5sum $DATA/tests/fractal_worley_noise.dag $DATA/tests/fractal_worley_noise.toml
./cubiquity export bin $DATA/tests/fractal_worley_noise.dag --output=$DST/fractal_worley_noise.bin
echo
md5sum $DST/fractal_worley_noise.bin $DST/fractal_worley_noise.txt
echo

echo "Exporting png image slices..."
echo "-----------------------------"
mkdir -p -- "$DST"/images # FIXME - Application should do this (if it doesn't?)
md5sum $DATA/tests/fractal_worley_noise.dag $DATA/tests/fractal_worley_noise.toml
./cubiquity export pngs $DATA/tests/fractal_worley_noise.dag --output="$DST"/images
echo
echo "Number of files is $(ls -1q "$DST"/images | wc -l)"
md5sum "$DST"/images/000000.png
md5sum "$DST"/images/000499.png
echo

echo "Exporting MagicaVoxel file..."
echo "-----------------------------"
md5sum $DATA/tests/fractal_worley_noise.dag $DATA/tests/fractal_worley_noise.toml
./cubiquity export vox $DATA/tests/fractal_worley_noise.dag --output=$DST/fractal_worley_noise.vox
echo
md5sum $DST/fractal_worley_noise.vox
echo

echo "========================================================================="
echo " Running view tests..."
echo "========================================================================="
echo

echo "Viewing fractal worley noise..."
echo "-------------------------------"
md5sum $DATA/tests/fractal_worley_noise.dag $DATA/tests/fractal_worley_noise.toml
./cubiquity view $DATA/tests/fractal_worley_noise.dag --width=1600 --height=1200 --duration=10
echo

echo "========================================================================="
echo " Running voxelization tests..."
echo "========================================================================="
echo

echo "Voxelizing Glycon..."
echo "--------------------"
md5sum $DATA/examples/glycon_decimated_plus_materials.obj $DATA/examples/glycon_decimated_plus_materials.mtl
./cubiquity voxelize $DATA/examples/glycon_decimated_plus_materials.obj --output=$DST/glycon_decimated_plus_materials.dag
echo
md5sum $DST/glycon_decimated_plus_materials.dag $DST/glycon_decimated_plus_materials.toml
echo
