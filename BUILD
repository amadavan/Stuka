
config_setting(
    name = "darwin",
    values = {"cpu": "darwin"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "windows",
    values = {"cpu": "x64_windows"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "linux_x86_64",
    values = {"cpu": "k8"},
    visibility = ["//visibility:public"],
)

load(
    "//third_party/gurobi:build_defs.bzl",
    "if_gurobi",
    "gurobi_deps"
)

cc_library(
    name = "stuka",
    srcs = glob(["src/**/*.cpp"], exclude=["src/example/**/*.cc", "src/example/**/*.cpp", "src/test.cc"]),
    hdrs = glob(["include/stuka/**/*.h"]),
    strip_include_prefix = "include",
    defines = if_gurobi(["ENABLE_GUROBI"]),
    deps = [
        "@eigen//:eigen",
        "@suitesparse//:cholmod",
        "@suitesparse//:spqr",
        "@hdf5//:hdf5",
        "@hdf5//:hdf5_cpp",
        ] + gurobi_deps(),
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "tester",
    srcs = glob(["src/example/**/*.cc", "src/example/**/*.cpp", "src/example/**/*.h"]) + ["src/test.cc"],
    deps = [":stuka"],
    visibility = ["//visibility:public"],
)

load("@stuka//:stuka.bzl", "pybind_extension")

pybind_extension(
    name = "stuka_extension",
    module_name = "stukapy",
    srcs = glob(["bindings/python/src/*.c"]),
    hdrs = glob(["bindings/python/src/*.h"]),
    deps = [":stuka"],
    visibility = ["//visibility:public"],
)

py_library(
    name = "stukapy",
    srcs_version = "PY2AND3",
    deps = [":stuka_extension"],
    visibility = ["//visibility:public"],
)