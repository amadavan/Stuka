# Rule for simple expansion of template files. This performs a simple
# search over the template file for the keys in substitutions,
# and replaces them with the corresponding values.
#
# Typical usage:
#   load("/tools/build_rules/template_rule", "expand_header_template")
#   template_rule(
#       name = "ExpandMyTemplate",
#       src = "my.template",
#       out = "my.txt",
#       substitutions = {
#         "$VAR1": "foo",
#         "$VAR2": "bar",
#       }
#   )
#
# Args:
#   name: The name of the rule.
#   template: The template file to expand
#   out: The destination of the expanded file
#   substitutions: A dictionary mapping strings to their substitutions

def template_rule_impl(ctx):
    subs = dict(ctx.attr.substitutions)

    if ctx.attr.expand_make:
        for key in subs.keys():
            subs[key] = ctx.expand_make_variables("substitutions", subs[key], ctx.var)
            subs[key] = ctx.expand_location(subs[key])

    ctx.actions.expand_template(
        template = ctx.file.src,
        output = ctx.outputs.out,
        substitutions = subs,
    )

template_rule = rule(
    attrs = {
        "src": attr.label(
            mandatory = True,
            allow_single_file = True,
        ),
        "substitutions": attr.string_dict(mandatory = True),
        "expand_make": attr.bool(default = False),
        "out": attr.output(mandatory = True),
    },
    # output_to_genfiles is required for header files.
    output_to_genfiles = True,
    implementation = template_rule_impl,
)