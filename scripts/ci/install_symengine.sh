#!/usr/bin/env bash
set -euo pipefail

revision=fac9314c78f2809570494017efc6603befeb4eda
source_dir=$(mktemp -d)
build_dir=$(mktemp -d)
trap 'rm -rf "${source_dir}" "${build_dir}"' EXIT

git clone --quiet https://github.com/symengine/symengine.git "${source_dir}"
git -C "${source_dir}" checkout --quiet "${revision}"
cmake -S "${source_dir}" -B "${build_dir}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_TESTS=OFF \
    -DBUILD_BENCHMARKS=OFF \
    -DINTEGER_CLASS=gmp \
    -DWITH_FLINT=ON \
    -DWITH_SYMENGINE_RCP=ON
cmake --build "${build_dir}" -j 2
sudo cmake --install "${build_dir}"
sudo ldconfig
