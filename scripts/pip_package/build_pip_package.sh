#!/usr/bin/env bash
# Copyright 2015 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================


set -e

function is_absolute {
  [[ "$1" = /* ]] || [[ "$1" =~ ^[a-zA-Z]:[/\\].* ]]
}

function real_path() {
  is_absolute "$1" && echo "$1" || echo "$PWD/${1#./}"
}

function cp_external() {
  local src_dir=$1
  local dest_dir=$2

  pushd .
  cd "$src_dir"
  for f in `find . ! -type d ! -name '*.py' ! -path '*stuka*'`; do
    mkdir -p "${dest_dir}/$(dirname ${f})"
    cp "${f}" "${dest_dir}/$(dirname ${f})/"
  done
  popd
}

function move_to_root_if_exists () {
  arg_to_move="$1"
  if [ -e "${arg_to_move}" ]; then
    mv ${arg_to_move} ./
  fi
}

PLATFORM="$(uname -s | tr 'A-Z' 'a-z')"
function is_windows() {
  if [[ "${PLATFORM}" =~ (cygwin|mingw32|mingw64|msys)_nt* ]]; then
    true
  else
    false
  fi
}

function prepare_src() {
  if [ $# -lt 1 ] ; then
    echo "No destination dir provided"
    exit 1
  fi

  TMPDIR="${1%/}"
  mkdir -p "$TMPDIR"
#  EXTERNAL_INCLUDES="${TMPDIR}/include/external"

  echo $(date) : "=== Preparing sources in dir: ${TMPDIR}"

  if [ ! -d bazel-bin ]; then
    echo "Could not find bazel-bin.  Did you run from the root of the build tree?"
    exit 1
  fi

  RUNFILES=bazel-bin/tools/pip_package/build_pip_package.runfiles/stuka
  cp -R bazel-bin/scripts/pip_package/build_pip_package.runfiles/stuka "${TMPDIR}"
#  cp_external bazel-bin/scripts/pip_package/build_pip_package.runfiles "${EXTERNAL_INCLUDES}"
#  # Copy MKL libs over so they can be loaded at runtime
#  so_lib_dir=$(ls $RUNFILES | grep solib) || true
#  if [ -n "${so_lib_dir}" ]; then
#  mkl_so_dir=$(ls ${RUNFILES}/${so_lib_dir} | grep mkl) || true
#  if [ -n "${mkl_so_dir}" ]; then
#  mkdir "${TMPDIR}/${so_lib_dir}"
#  cp -R ${RUNFILES}/${so_lib_dir}/${mkl_so_dir} "${TMPDIR}/${so_lib_dir}"
#  fi
#  fi


  cp scripts/pip_package/MANIFEST.in ${TMPDIR}
  cp scripts/pip_package/README ${TMPDIR}
  cp scripts/pip_package/setup.py ${TMPDIR}


  mkdir -p ${TMPDIR}/build
  mkdir -p ${TMPDIR}/build/lib
  cp -r ${TMPDIR}/stuka/* ${TMPDIR}/build/lib/

#  rm -f ${TMPDIR}/stuka_extension.so
#  rm -f ${TMPDIR}/libstuka.so.[0-9].*
}

function build_wheel() {
  if [ $# -lt 2 ] ; then
    echo "No src and dest dir provided"
    exit 1
  fi

  TMPDIR="$1"
  DEST="$2"
  PKG_NAME_FLAG="$3"

#  # Before we leave the top-level directory, make sure we know how to
#  # call python.
#  if [[ -e tools/python_bin_path.sh ]]; then
#    source tools/python_bin_path.sh
#  fi

  pushd ${TMPDIR} > /dev/null

  rm -f MANIFEST
  echo $(date) : "=== Building wheel"
  python3 setup.py bdist_wheel
  mkdir -p ${DEST}
  cp dist/* ${DEST}
  popd > /dev/null
  echo $(date) : "=== Output wheel file is in: ${DEST}"
}

function usage() {
  echo "Usage:"
  echo "$0 [--src srcdir] [--dst dstdir] [options]"
  echo "$0 dstdir [options]"
  echo ""
  echo "    --src                 prepare sources in srcdir"
  echo "                              will use temporary dir if not specified"
  echo ""
  echo "    --dst                 build wheel in dstdir"
  echo "                              if dstdir is not set do not build, only prepare sources"
  echo ""
  exit 1
}

function main() {
  PKG_NAME_FLAG=""
  PROJECT_NAME="stuka"
  SRCDIR=""
  DSTDIR=""
  CLEANSRC=1
  while true; do
    if [[ "$1" == "--help" ]]; then
      usage
      exit 1
    elif [[ "$1" == "--src" ]]; then
      shift
      SRCDIR="$(real_path $1)"
      CLEANSRC=0
    elif [[ "$1" == "--dst" ]]; then
      shift
      DSTDIR="$(real_path $1)"
    else
      DSTDIR="$(real_path $1)"
    fi
    shift

    if [[ -z "$1" ]]; then
      break
    fi
  done

  if [[ -z "$DSTDIR" ]] && [[ -z "$SRCDIR" ]]; then
    echo "No destination dir provided"
    usage
    exit 1
  fi

  if [[ -z "$SRCDIR" ]]; then
    # make temp srcdir if none set
    SRCDIR="$(mktemp -d -t tmp.XXXXXXXXXX)"
  fi

  prepare_src "$SRCDIR"

  if [[ -z "$DSTDIR" ]]; then
      # only want to prepare sources
      exit
  fi

  if [[ -n ${PROJECT_NAME} ]]; then
    PKG_NAME_FLAG="--project_name ${PROJECT_NAME}"
  fi

  build_wheel "$SRCDIR" "$DSTDIR" "$PKG_NAME_FLAG"

  if [[ $CLEANSRC -ne 0 ]]; then
    rm -rf "${TMPDIR}"
  fi
}

main "$@"

