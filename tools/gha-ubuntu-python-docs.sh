#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install wget

# Install conda+deps.
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -y -q -p $deps_dir c-compiler cxx-compiler cmake fmt boost boost-cpp symengine pagmo-devel python=3.10 sphinx nbsphinx pybind11 pybind11-abi
source activate $deps_dir

# Create the build dir and cd into it.
mkdir build

# Install dsyre headers and library
cd build
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Release \
    -DDSYRE_BUILD_TESTS=no \
    -DDSYRE_BUILD_MAIN=no \
    -DDSYRE_BUILD_PYTHON=no \
    ..
make VERBOSE=1 install

# Install py-dsyre
cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DCMAKE_BUILD_TYPE=Release \
    -DDSYRE_BUILD_TESTS=no \
    -DDSYRE_BUILD_MAIN=no \
    -DDSYRE_BUILD_PYTHON=yes \
    ..
make VERBOSE=1 install

# Test python installation
python -c "import pydsyre as dsy; dsy.test.run_test_suite()"

# Switch to gh-pages branch
git checkout gh-pages

# Build docs
cd docs
make html

set +e
set +x
