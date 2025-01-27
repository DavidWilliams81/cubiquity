#!/bin/bash

# Temp dirs
build_dir=$(mktemp -d)
working_dir=$(mktemp -d)

# Run CMake
cmake -S ..  \
      -B $build_dir \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_TOOLCHAIN_FILE=utils/mingw-w64-x86_64.cmake \
      -DSDL2_DIR=~/SDL2-2.30.11/cmake/

# Make (in subshell, to preserve current directory)
(cd $build_dir && make -j8)

# Assemble folder with files to be packaged
cp $build_dir/cubiquity.exe $working_dir
cp -r ../src/application/commands/view/glsl/ $working_dir
wget -P $working_dir https://cubiquity.s3.eu-west-2.amazonaws.com/release_template/glycon.obj
wget -P $working_dir https://cubiquity.s3.eu-west-2.amazonaws.com/release_template/ReadMe.md
wget -P $working_dir https://cubiquity.s3.eu-west-2.amazonaws.com/release_template/shapes.obj
wget -P $working_dir https://cubiquity.s3.eu-west-2.amazonaws.com/release_template/shapes.mtl

# Zip files up in a filename based on current date
# https://askubuntu.com/questions/521011/zip-an-archive-without-including-parent-directory
release_name=cubiquity-win64-$(date +%Y-%m-%d)
(cd $working_dir && zip -r - .) > $release_name.zip
sha256sum $release_name.zip > $release_name.sha256

# Upload to S3 release folder
s3_release_dir=s3://cubiquity/releases/
aws s3 cp --content-type application/zip $release_name.zip $s3_release_dir
aws s3 cp --content-type text/plain $release_name.sha256 $s3_release_dir

# Clean up
rm -rf "$build_dir"
rm -rf "$working_dir"

# Print out Markdown for copying into forum posts, etc.
release_url=https://cubiquity.s3.eu-west-2.amazonaws.com/releases/
echo "[Download Cubiquity for Windows]($release_url$release_name.zip) ([sha256]($release_url$release_name.sha256))"
