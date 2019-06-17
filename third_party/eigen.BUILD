# Description:
#   Eigen is a C++ template library for linear algebra: vectors,
#   matrices, and related algorithms.

licenses([
    "notice"
])

exports_files(["COPYING.MPL2"])

# Files known to be under MPL2 license.
EIGEN_HEADER_FILES = glob(["Eigen/**"])

cc_library(
    name = "eigen",
    hdrs = EIGEN_HEADER_FILES,
    defines = [
        # This define (mostly) guarantees we don't link any problematic
        # code. We use it, but we do not rely on it, as evidenced above.
#        "EIGEN_MPL2_ONLY",
        "EIGEN_MAX_ALIGN_BYTES=64",
        "EIGEN_HAS_TYPE_TRAITS=0",
    ],
    includes = ["."],
    visibility = ["//visibility:public"],
)

filegroup(
    name = "eigen_header_files",
    srcs = EIGEN_HEADER_FILES,
    visibility = ["//visibility:public"],
)