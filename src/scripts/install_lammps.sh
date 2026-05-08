#!/bin/bash
set -e

LAMMPS_VERSION="patch_11Feb2026"
INSTALL_DIR="${CONDA_PREFIX}"

echo "Building LAMMPS ${LAMMPS_VERSION} from source..."

# Clone
git clone --depth 1 --branch ${LAMMPS_VERSION} \
    https://github.com/lammps/lammps.git /tmp/lammps

# Build
mkdir -p /tmp/lammps/build && cd /tmp/lammps/build

cmake ../cmake \
    -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
    -DLAMMPS_EXCEPTIONS=yes \
    -DBUILD_SHARED_LIBS=yes \
    -DBUILD_LIB=yes \
    -DPKG_MOLECULE=yes \
    -DPKG_KSPACE=yes \
    -DPKG_EXTRA-MOLECULE=yes \
    -DPKG_RIGID=yes \
    -DPKG_MANYBODY=yes \
    -DPKG_GRAPHICS=yes \

make -j$(nproc)
make install

# Install Python bindings into active env
cd /tmp/lammps/python
pip install .

echo "LAMMPS installed successfully"