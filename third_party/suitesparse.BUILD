
licenses(["notice"])

cc_library(
    name = "metis",
    srcs = glob(["metis-5.1.0/libmetis/*.c", "metis-5.1.0/GKlib/*.c"]),
    hdrs = glob(["metis-5.1.0/libmetis/*.h", "metis-5.1.0/GKlib/*.h"]),
    includes = ["metis-5.1.0/libmetis", "metis-5.1.0/include", "metis-5.1.0/GKlib"],
    copts = ["-fPIC"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "suitesparse_config",
    srcs = glob(["SuiteSparse_config/*.c"]),
    hdrs = glob(["SuiteSparse_config/*.h"]),
    strip_include_prefix = "SuiteSparse_config",
    includes = ["SuiteSparse_config"],
    copts = ["-fPIC"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "amd_i",
    srcs = glob(["AMD/Source/amd_*.c"]),
    hdrs = glob(["AMD/Include/*.h"]),
    strip_include_prefix = "AMD/Include",
    copts = ["-fPIC", "-DDINT"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "amd_l",
    srcs = glob(["AMD/Source/amd_*.c"]),
    hdrs = glob(["AMD/Include/*.h"]),
    strip_include_prefix = "AMD/Include",
    copts = ["-fPIC", "-DDLONG"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "amd",
    deps = [":amd_i", ":amd_l"],
    includes = ["AMD/Include", "AMD/Source"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "camd_i",
    srcs = glob(["CAMD/Source/camd_*.c"]),
    hdrs = glob(["CAMD/Include/*.h"]),
    strip_include_prefix = "CAMD/Include",
    copts = ["-fPIC", "-DDINT"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "camd_l",
    srcs = glob(["CAMD/Source/camd_*.c"]),
    hdrs = glob(["CAMD/Include/*.h"]),
    strip_include_prefix = "CAMD/Include",
    copts = ["-fPIC", "-DDLONG"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "camd",
    deps = [":camd_i", ":camd_l"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "colamd_i",
    srcs = glob(["COLAMD/Source/colamd.c"]),
    hdrs = glob(["COLAMD/Include/*.h"]),
    strip_include_prefix = "COLAMD/Include",
    copts = ["-fPIC", "-DDINT"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "colamd_l",
    srcs = glob(["COLAMD/Source/colamd.c"]),
    hdrs = glob(["COLAMD/Include/*.h"]),
    strip_include_prefix = "COLAMD/Include",
    copts = ["-fPIC", "-DDLONG"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "colamd",
    deps = [":colamd_i", ":colamd_l"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "ccolamd_i",
    srcs = glob(["CCOLAMD/Source/ccolamd.c"]),
    hdrs = glob(["CCOLAMD/Include/*.h"]),
    strip_include_prefix = "CCOLAMD/Include",
    copts = ["-fPIC", "-DDINT"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "ccolamd_l",
    srcs = glob(["CCOLAMD/Source/ccolamd.c"]),
    hdrs = glob(["CCOLAMD/Include/*.h"]),
    strip_include_prefix = "CCOLAMD/Include",
    copts = ["-fPIC", "-DDLONG"],
    deps = [":suitesparse_config"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "ccolamd",
    deps = [":ccolamd_i", ":ccolamd_l"],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "cholmod_i",
    srcs = glob([
        "CHOLMOD/Core/*.c",
        "CHOLMOD/Check/*.c",
        "CHOLMOD/Cholesky/*.c",
        "CHOLMOD/MatrixOps/*.c",
        "CHOLMOD/Partition/*.c",
        "CHOLMOD/Modify/*.c",
        "CHOLMOD/Supernodal/*.c",
    ], exclude = ["CHOLMOD/**/t_*.c"]),
    hdrs = glob(["CHOLMOD/Include/*.h"]),
    strip_include_prefix = "CHOLMOD/Include",
    includes = ["CHOLMOD/Core", "CHOLMOD/Cholesky", "CHOLMOD/MatrixOps", "CHOLMOD/Partition", "CHOLMOD/Modify", "CHOLMOD/Supernodal"],
    copts = ["-fPIC", "-DDINT", "-fexceptions", "-fno-common"],
    deps = [
        ":metis",
        ":suitesparse_config",
        ":amd",
        ":colamd",
        ":camd",
        ":ccolamd"
    ],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "cholmod_l",
    srcs = glob([
        "CHOLMOD/Core/*.c",
        "CHOLMOD/Check/*.c",
        "CHOLMOD/Cholesky/*.c",
        "CHOLMOD/MatrixOps/*.c",
        "CHOLMOD/Partition/*.c",
        "CHOLMOD/Modify/*.c",
        "CHOLMOD/Supernodal/*.c",
    ], exclude = ["CHOLMOD/**/t_*.c"]),
    hdrs = glob(["CHOLMOD/Include/*.h"]),
    strip_include_prefix = "CHOLMOD/Include",
    includes = ["CHOLMOD/Core", "CHOLMOD/Cholesky", "CHOLMOD/MatrixOps", "CHOLMOD/Partition", "CHOLMOD/Modify", "CHOLMOD/Supernodal"],
    copts = ["-fPIC", "-DDLONG", "-fexceptions", "-fno-common"],
    deps = [
        ":metis",
        ":suitesparse_config",
        ":amd",
        ":colamd",
        ":camd",
        ":ccolamd"
    ],
    visibility = ["//visibility:private"],
)

cc_library(
    name = "cholmod",
    deps = [":cholmod_i", ":cholmod_l", "@openblas//:openblas"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "umfpack",
    srcs = glob(["UMFPACK/Source/*.c"]),
    hdrs = glob(["UMFPACK/Include/*.h", "UMFPACK/Source/*.h"]),
#    strip_include_prefix = "UMFPACK/Include",
    includes = ["UMFPACK/Include"],
    copts = ["-fPIC"],
    deps = [
        ":suitesparse_config",
        ":amd",
        ":cholmod",
        "@openblas//:openblas",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "spqr",
    srcs = glob(["SPQR/Source/*.cpp"]),
    hdrs = glob(["SPQR/Include/*.h", "SPQR/Include/*.hpp"]),
    strip_include_prefix = "SPQR/Include",
    copts = ["-fPIC"],
    deps = [
        ":suitesparse_config",
        ":amd",
        ":colamd",
        ":cholmod",
        "@openblas//:openblas",
    ],
    visibility = ["//visibility:public"],
)