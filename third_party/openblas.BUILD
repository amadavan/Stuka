
licenses(["notice"])

genrule(
    name = "build_openblas",
    srcs = glob(["**/*"]),
    outs = ["libopenblas.a"],
    cmd = "".join([
        "pushd external/openblas/; workdir=$$(mktemp -d -t tmp.XXXXXXXXXX); cp -r * $$workdir; pushd $$workdir;",
        "make;",
        "make PREFIX=. install;",
        "popd; popd;",
        "cp $$(find $$workdir/lib -name 'libopenblas*.a') $(@D)/;",
        "rm -rf $$workdir;"
    ]),
    message = "Building OpenBLAS",
    visibility = ["//visibility:private"],
)

cc_library(
    name = "openblas",
    srcs = ["libopenblas.a"],
    visibility = ["//visibility:public"],
)