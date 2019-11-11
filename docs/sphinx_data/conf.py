# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Stuka'
copyright = '2019, Avinash N. Madavan'
author = 'Avinash N. Madavan'

# The full version, including alpha/beta/rc tags
release = '1.0.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "breathe",
    "recommonmark",
    'sphinx.ext.mathjax',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = ['.rst', '.md']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for imgmath -----------------------------------------------------

imgmath_image_format = 'svg'

imgmath_latex_preamble = r'''
    \usepackage{amsmath,amssymb,amsfonts,mathrsfs,mathtools,dsfont,bm,bbm}
    \DeclareMathAlphabet{\mathpzc}{OT1}{pzc}{m}{it}
    \renewcommand{\v}[1]{{\bm{#1}}}
'''

# -- Options for MathJAX -----------------------------------------------------

mathjax_config = {
    'TeX':
        {
            'Macros': {
                'bm': ["\\boldsymbol{#1}", 1],
                'v': ['{\\bm{#1}}', 1]
            }
        }
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'press'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Options for breathe -----------------------------------------------------
breathe_default_project = "stuka"
breathe_default_members = ('members', 'undoc-members')
