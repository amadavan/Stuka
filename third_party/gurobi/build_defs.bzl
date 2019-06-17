
_GUROBI_HOME = "GUROBI_HOME"

def if_gurobi(if_true, if_false = []):
    """Shorthand for select()'ing on whether we're building with GUROBI.
    Args:
      if_true: expression to evaluate if building with GUROBI.
      if_false: expression to evaluate if building without GUROBI.
    Returns:
      a select evaluating to either if_true or if_false as appropriate.
    """
    return select({
        str(Label("//third_party/gurobi:enable_gurobi")): if_true,
        "//conditions:default": if_false,
    })

def gurobi_deps():
    """Shorthand for select() to pull in the correct set of GUROBI library deps.
    Returns:
      a select evaluating to a list of library dependencies, suitable for
      inclusion in the deps attribute of rules.
    """
    return select({
        str(Label("//third_party/gurobi:enable_gurobi")): ["//third_party/gurobi:binary_blob"],
        "//conditions:default": [],
    })

def _enable_local_gurobi(repository_ctx):
    return _GUROBI_HOME in repository_ctx.os.environ

def _gurobi_autoconf_impl(repository_ctx):
    """Implementation of the local_gurobi_autoconf repository rule."""

    if _enable_local_gurobi(repository_ctx):
        # Symlink lib and include local folders.
        gurobi_root = repository_ctx.os.environ[_GUROBI_HOME]
        gurobi_lib_path = "%s/lib" % gurobi_root
        repository_ctx.symlink(gurobi_lib_path, "lib")
        gurobi_include_path = "%s/include" % gurobi_root
        repository_ctx.symlink(gurobi_include_path, "include")

    # Also setup BUILD file.
    repository_ctx.symlink(repository_ctx.attr.build_file, "BUILD")

gurobi_repository = repository_rule(
    implementation = _gurobi_autoconf_impl,
    environ = [
        _GUROBI_HOME
    ],
    attrs = {
        "build_file": attr.label(),
        "strip_prefix": attr.string(default = ""),
    },
)