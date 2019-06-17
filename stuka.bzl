
def pybind_extension(
        name,
        srcs,
        module_name,
        hdrs = [],
        features = [],
        srcs_version = "PY2AND3",
        data = [],
        copts = None,
        nocopts = None,
        linkopts = [],
        deps = [],
        visibility = None,
        testonly = None,
        licenses = None,
        compatible_with = None,
        restricted_to = None,
        deprecation = None):
    """Builds a Python extension module."""
    _ignore = [module_name]
    p = name.rfind("/")
    if p == -1:
        sname = name
        prefix = ""
    else:
        sname = name[p + 1:]
        prefix = name[:p + 1]
    so_file = "%s%s.so" % (prefix, sname)
    pyd_file = "%s%s.pyd" % (prefix, sname)
    symbol = "init%s" % sname
    symbol2 = "init_%s" % sname
    symbol3 = "PyInit_%s" % sname
    exported_symbols_file = "%s-exported-symbols.lds" % name
    version_script_file = "%s-version-script.lds" % name
    native.genrule(
        name = name + "_exported_symbols",
        outs = [exported_symbols_file],
        cmd = "echo '_%s\n_%s\n_%s' >$@" % (symbol, symbol2, symbol3),
        output_licenses = ["unencumbered"],
        visibility = ["//visibility:private"],
        testonly = testonly,
    )

    native.genrule(
        name = name + "_version_script",
        outs = [version_script_file],
        cmd = "echo '{global:\n %s;\n %s;\n %s;\n local: *;};' >$@" % (symbol, symbol2, symbol3),
        output_licenses = ["unencumbered"],
        visibility = ["//visibility:private"],
        testonly = testonly,
    )
    native.cc_binary(
        name = so_file,
        srcs = srcs + hdrs,
        data = data,
        copts = copts,
        nocopts = nocopts,
        linkopts = linkopts + select({
            "//:darwin": [
                "-Wl,-exported_symbols_list,$(location %s)" % exported_symbols_file,
            ],
            "//:windows": [],
            "//conditions:default": [
                "-Wl,--version-script",
                "$(location %s)" % version_script_file,
            ],
        }),
        deps = deps + [
            exported_symbols_file,
            version_script_file,
        ] + ["@pybind11//:pybind11"],
        features = features,
        linkshared = 1,
        testonly = testonly,
        licenses = licenses,
        visibility = visibility,
        deprecation = deprecation,
        restricted_to = restricted_to,
        compatible_with = compatible_with,
    )
    native.genrule(
        name = name + "_pyd_copy",
        srcs = [so_file],
        outs = [pyd_file],
        cmd = "cp $< $@",
        output_to_bindir = True,
        visibility = visibility,
        deprecation = deprecation,
        restricted_to = restricted_to,
        compatible_with = compatible_with,
    )
    native.py_library(
        name = name,
        data = select({
            "//:windows": [pyd_file],
            "//conditions:default": [so_file],
        }),
        srcs_version = srcs_version,
        licenses = licenses,
        testonly = testonly,
        visibility = visibility,
        deprecation = deprecation,
        restricted_to = restricted_to,
        compatible_with = compatible_with,
    )