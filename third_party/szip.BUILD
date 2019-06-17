
licenses(["notice"])

genrule(
    name = "configure_szip",
    srcs = glob(["**/*"]),
    outs = ["src/SZconfig.h"],
    cmd = "pushd external/szip/; workdir=$$(mktemp -d -t tmp.XXXXXXXXXX); cp -r * $$workdir; pushd $$workdir; ./configure; popd; popd; cp $$workdir/src/SZconfig.h $(@D)/SZconfig.h; rm -rf $$workdir;",
    tools = ["configure"],
    message = "Configuring SZIP",
    visibility = ["//visibility:private"],
)

cc_library(
    name = "szip",
    srcs = glob(["src/*.c"]),
    hdrs = glob(["src/*.h"], exclude=["src/SZconfig.h"]) + [":configure_szip"],
    copts = ["-O6", "-Wall", "-fomit-frame-pointer", "-funroll-loops"],
    includes = ["src"],
    visibility = ["//visibility:public"],
)