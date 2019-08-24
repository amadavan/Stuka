licenses(["notice"])

load("@rules_foreign_cc//tools/build_defs:configure.bzl", "configure_make")

filegroup(name = "all", srcs = glob(["**"]), visibility = ["//visibility:private"])

configure_make(
    name = "libflame",
    lib_name = "libflame",
    lib_source = ":all",
    configure_options = [
        "--enable-lapack2flame",
        "--enable-cblas-interfaces",
        "--enable-max-arg-list-hack",
    ],
    deps = [
        "@blis//:blis"
    ],
    visibility = ["//visibility:public"]
)

## Copies environment to tmp and runs ./configure so we can get the
## configuration header file
#genrule(
#  name = "configure_flame",
#  srcs = glob(["**/*"]),
#  outs = ["include/flame/flame.h", "include/flame/cblas.h", "lib/libblas.a"],
#  cmd = "".join([
#      "pushd external/hdf5/;",
#      "workdir=$$(mktemp -d -t tmp.XXXXXXXXXX);",
#      "cp -r * $$workdir; pushd $$workdir;",
#      "./configure --prefix=$(pwd) --enable-cblas --enable-lapack2flame;",
#      "make -j4; make -j4 install;",
#      "popd; popd;",
#      "cp $$workdir/include/flame/flame.h $(@D)/include/flame/flame.h;",
#      "cp $$workdir/include/flame/cblas.h $(@D)/include/flame/cblas.h;",
#      "cp $$workdir/lib/libflame.a $(@D)/lib/libflame.a;",
#      "rm -rf $$workdir;"
#      ]),
#  tools = ["configure"],
#  message = "Configuring Flame",
#  visibility = ["//visibility:private"],
#)
#
#cc_library(
#    name = "flame",
#    srcs = ":configure_flame",
#    includes = glob(["include/flame/*"]),
#    strip_include_prefix = "include/flame",
#)