#!/bin/sh
tar -cJf cubiquity-data.tar.xz --exclude='archive_data.sh' --exclude='create_data.py' --exclude='get_data.py' --exclude='cubiquity-data.tar.xz' *
