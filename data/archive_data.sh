#!/bin/sh
tar -cJf cubiquity-data.tar.xz --exclude='archive_data.sh' --exclude='cubiquity-data.tar.xz' --exclude='get_data.py' *
