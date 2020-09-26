#! /bin/sh

set -ex

# Get current directory
DIR="$( cd "$(dirname "$0")" ; pwd -P )"

# Download and unzip SuiteSparse
curl -L https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/v5.8.1.tar.gz -o SuiteSparse-5.8.1.tar.gz
tar -zxvf SuiteSparse-5.8.1.tar.gz

# Copy modified Makefile to reduce compilation time (removes unnecessary libs)
cp scripts/Makefile SuiteSparse-5.8.1/

# Build the library and install in Stuka directory
cd SuiteSparse-5.8.1
make library
make install INSTALL=. \
    INSTALL_DOC=SuiteSparse/doc \
    INSTALL_LIB=${DIR}/../lib \
    INSTALL_INCLUDE=${DIR}/../include/suitesparse
