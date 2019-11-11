def _doxygen_impl(ctx):

    dir = ctx.expand_make_variables("substitutions", ctx.attr.output_dir, ctx.var)
    tree = ctx.actions.declare_directory(dir, sibling=ctx.files.doxyfile[0])

    ctx.actions.run(
        inputs = ctx.attr.source[DefaultInfo].data_runfiles.files,
        outputs = [tree],
        arguments = [ctx.files.doxyfile[0].path],
        executable = "doxygen",
        use_default_shell_env=True
    )

    runfiles = ctx.runfiles(ctx.files.doxyfile)

    return [ DefaultInfo(files = depset([ tree ]), runfiles = runfiles) ]

doxygen_rule = rule(
    implementation = _doxygen_impl,
    attrs = {
        "doxyfile": attr.label(
            mandatory = True,
            allow_single_file = True),
        "output_dir": attr.string(default = "doxygen.out"),
        "source": attr.label(mandatory = True),
    },
)