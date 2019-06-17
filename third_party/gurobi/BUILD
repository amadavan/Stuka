licenses(["notice"])  # 3-Clause BSD

config_setting(
    name = "enable_gurobi",
    define_values = {
        "enable_gurobi": "true",
    },
    visibility = ["//visibility:public"],
)

load("//third_party/gurobi:build_defs.bzl", "if_gurobi")

cc_library(
    name = "binary_blob",
    deps = select({
        "@stuka//:linux_x86_64": [
            "@gurobi//:gurobi_headers",
            "@gurobi//:gurobi_libs_linux",
        ],
        "@stuka//:darwin": [
            "@gurobi//:gurobi_headers",
            "@gurobi//:gurobi_libs_darwin",
        ],
        "@stuka//:windows": [
            "@gurobi//:gurobi_headers",
            "@gurobi//:gurobi_libs_windows",
        ],
        "//conditions:default": [],
    }),
    visibility = ["//visibility:public"],
)