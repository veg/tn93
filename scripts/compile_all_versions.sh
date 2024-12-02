#!/usr/bin/env bash
# compile all versions of tn93 into the current working directory
# USAGE: ./compile_all_versions [MAKE_FLAGS]
wget -qO- "https://api.github.com/repos/veg/tn93/releases" | grep '"tarball_url"' | cut -d'"' -f4 | while read url ; do
    version=$(echo "$url" | rev | cut -d'/' -f1 | rev)
    wget -qO- "$url" | tar -zx
    cd veg-tn93-*
    cmake .
    make $@
    mv tn93 "../tn93.$version"
    cd ..
    rm -rf veg-tn93-*
done
