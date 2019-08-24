licenses(["notice"])

load("@rules_foreign_cc//tools/build_defs:configure.bzl", "configure_make")

filegroup(name = "all", srcs = glob(["**"]), visibility = ["//visibility:private"])

configure_make(
    name = "blis",
    lib_name = "libblis",
    lib_source = ":all",
    configure_options = [
        "--enable-cblas",
        "auto",
    ],
    configure_env_vars = {},
    visibility = ["//visibility:public"]
)

## Copies environment to tmp and runs ./configure so we can get the
## configuration header file
#genrule(
#  name = "configure_blis",
#  srcs = glob(["**/*"]),
#  outs = ["include/blis/blis.h", "include/blis/cblas.h", "lib/libblas.a"],
#  cmd = "".join([
#      "pushd external/hdf5/;",
#      "workdir=$$(mktemp -d -t tmp.XXXXXXXXXX);",
#      "cp -r * $$workdir; pushd $$workdir;",
#      "./configure -p . --enable-cblas auto;",
#      "make -j4; make -j4 install;",
#      "popd; popd;",
#      "cp $$workdir/include/blis/blis.h $(@D)/include/blis/blis.h;",
#      "cp $$workdir/include/blis/cblas.h $(@D)/include/blis/cblas.h;",
#      "cp $$workdir/lib/libblis.a $(@D)/lib/libblis.a;",
#      "rm -rf $$workdir;"
#      ]),
#  tools = ["configure"],
#  message = "Configuring BLIS",
#  visibility = ["//visibility:private"],
#)
#
#cc_library(
#    name = "blis",
#    srcs = ":configure_blis",
#    includes = glob(["include/blis/*"]),
#    strip_include_prefix = "include/blis",
#)