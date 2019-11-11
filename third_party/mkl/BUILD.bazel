licenses(["notice"])  # 3-Clause BSD

config_setting(
    name = "enable_mkl",
    define_values = {
        "enable_mkl": "true",
    },
    visibility = ["//visibility:public"],
)

load(
    "//third_party/mkl:build_defs.bzl",
    "if_mkl",
)

filegroup(
    name = "LICENSE",
    srcs = ["MKL_LICENSE"] + select({
        "@stuka//:linux_x86_64": [
            "@mkl_linux//:LICENSE",
        ],
        "@stuka//:darwin": [
            "@mkl_darwin//:LICENSE",
        ],
        "@stuka//:windows": [
            "@mkl_windows//:LICENSE",
        ],
        "//conditions:default": [],
    }),
    visibility = ["//visibility:public"],
)

cc_library(
    name = "intel_binary_blob",
    visibility = ["//visibility:public"],
    deps = select({
        "@stuka//:linux_x86_64": [
            "@mkl_linux//:mkl_headers",
            "@mkl_linux//:mkl_libs_linux",
        ],
        "@stuka//:darwin": [
            "@mkl_darwin//:mkl_headers",
            "@mkl_darwin//:mkl_libs_darwin",
        ],
        "@stuka//:windows": [
            "@mkl_windows//:mkl_headers",
            "@mkl_windows//:mkl_libs_windows",
        ],
        "//conditions:default": [],
    }),
)