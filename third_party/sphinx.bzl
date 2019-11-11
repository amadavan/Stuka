def _sphinx_impl(ctx):
    output_dir = ctx.expand_make_variables("substitutions", ctx.attr.output_dir, ctx.var)
    tree = ctx.actions.declare_directory(output_dir, sibling=ctx.file.doxygen_output)

    doxygen_output = ctx.attr.doxygen_output[DefaultInfo].data_runfiles.files
    srcs = ctx.attr.srcs[DefaultInfo].data_runfiles.files

    # The doxygen directory is sent to the sphinx build file as a relative path.
    # Create an appropriately-sized back-path string to compensate.
    precursor = ''
    for src in srcs.to_list():
        if src.basename == "conf.py":
            precursor = '../' * len(src.dirname.split('/'))

    data = depset(srcs.to_list() + [ctx.file.doxygen_output])

    for otype in ctx.attr.output_types:
        ctx.actions.run(
            inputs = data,
            outputs = [tree],
            arguments = [
                "-b", otype,
                "-Dbreathe_projects."+ctx.attr.breathe_project_name+"="+precursor+ctx.file.doxygen_output.path+"/xml",
                ctx.attr.source_dir,
                tree.path
            ],
            executable = "sphinx-build",
            use_default_shell_env=True
        )

    return [ DefaultInfo(files = depset([ tree ])) ]

sphinx_rule = rule(
    implementation = _sphinx_impl,
    attrs = {
        "breathe_project_name": attr.string(
            mandatory = True
        ),
        "doxygen_output": attr.label(default = "//docs:doxygen", allow_single_file=True),
        "output_dir": attr.string(default = "sphinx.out"),
        "source_dir": attr.string(default = ".", mandatory = True),
        "srcs": attr.label(mandatory = True),
        "output_types": attr.string_list(default = ["html"])
    },
)