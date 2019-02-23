#! /bin/sh

set -ex

# Get current directory
DIR="$( cd "$(dirname "$0")" ; pwd -P )"

# Download and unzip SuiteSparse
curl http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.4.0.tar.gz -o SuiteSparse-5.4.0.tar.gz
tar -zxvf SuiteSparse-5.4.0.tar.gz

# Copy modified Makefile to reduce compilation time (removes unnecessary libs)
cp scripts/Makefile SuiteSparse/

# Build the library and install in Stuka directory
cd SuiteSparse
make library
make install INSTALL=. \
    INSTALL_DOC=SuiteSparse/doc \
    INSTALL_LIB=${DIR}/../lib \
    INSTALL_INCLUDE=${DIR}/../include/suitesparse