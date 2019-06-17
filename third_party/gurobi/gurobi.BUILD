
cc_library(
    name = "gurobi_headers",
    srcs = glob(["include/*"]),
    includes = ["include"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "gurobi_libs_linux",
    srcs = [
        "lib/libgurobi81.so",
        "lib/libgurobi_g++5.2.a",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "gurobi_libs_darwin",
    srcs = [
        "lib/libgurobi81.dylib",
        "lib/libgurobi_c++.a",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "gurobi_libs_windows",
    srcs = [
        "lib/gurobi.lib",
        "lib/gurobi_c++md2017.lib",
    ],
    visibility = ["//visibility:public"],
)