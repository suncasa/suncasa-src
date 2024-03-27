# conf.py
import os
import sys

sys.path.insert(0, os.path.abspath('../suncasa'))

# -- Project information -----------------------------------------------------

project = 'suncasa'
author = 'EOVSA Team'
coypright = '2024, EOVSA Team'

# The full version, including alpha/beta/rc tags
release = '1.0.6'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',  # Automatically document your code
    'sphinx.ext.napoleon',  # Support for Google-style docstrings
    'sphinx.ext.viewcode',  # Add links to source code from documentation
    'sphinx.ext.mathjax',  # Render math via JavaScript
    'sphinx.ext.autosummary',  # Automatically generates summary tables from the docstrings.
    'sphinx.ext.githubpages',
    'sphinx.ext.graphviz',
    'sphinx.ext.imgmath',
    'sphinx_gallery.gen_gallery',
    'sphinx_gallery',
    'autoapi.extension'
    # Add any other Sphinx extensions here.
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages. See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using a default sidebar with
# the links to the documentation's roots, contents & search, plus a link to
# the Python.org website.
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
        'donate.html',
    ]
}

# -- Options for sphinx.ext.autodoc -------------------------------------------
autoapi_type = 'python'
autoapi_dirs = ['../suncasa']
autoapi_ignore = [
    'tests/*',
    'suncasatasks/private/*',
    'suncasatasks/gotasks/*',
    'utils/grff/*'
]
