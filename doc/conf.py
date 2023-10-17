# scikit-learn documentation build configuration file, created by
# sphinx-quickstart on Fri Jan  8 09:13:42 2010.
#
# This file is execfile()d with the current directory set to its containing
# dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import os
import re
import sys
import warnings
from datetime import datetime
from io import StringIO
from pathlib import Path

from sklearn.externals._packaging.version import parse

# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory
# is relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
sys.path.insert(0, os.path.abspath("sphinxext"))

import sphinx_gallery
from github_link import make_linkcode_resolve
from sphinx_gallery.notebook import add_code_cell, add_markdown_cell
from sphinx_gallery.sorting import ExampleTitleSortKey

try:
    # Configure plotly to integrate its output into the HTML pages generated by
    # sphinx-gallery.
    import plotly.io as pio

    pio.renderers.default = "sphinx_gallery"
except ImportError:
    # Make it possible to render the doc when not running the examples
    # that need plotly.
    pass

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.linkcode",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.imgconverter",
    "sphinx_gallery.gen_gallery",
    "sphinx_issues",
    "add_toctree_functions",
    "sphinx-prompt",
    "sphinx_copybutton",
    "sphinxext.opengraph",
    "doi_role",
    "allow_nan_estimators",
    "matplotlib.sphinxext.plot_directive",
]

# Specify how to identify the prompt when copying code snippets
copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True
copybutton_exclude = "style"

try:
    import jupyterlite_sphinx  # noqa: F401

    extensions.append("jupyterlite_sphinx")
    with_jupyterlite = True
except ImportError:
    # In some cases we don't want to require jupyterlite_sphinx to be installed,
    # e.g. the doc-min-dependencies build
    warnings.warn(
        "jupyterlite_sphinx is not installed, you need to install it "
        "if you want JupyterLite links to appear in each example"
    )
    with_jupyterlite = False

# Produce `plot::` directives for examples that contain `import matplotlib` or
# `from matplotlib import`.
numpydoc_use_plots = True

# Options for the `::plot` directive:
# https://matplotlib.org/stable/api/sphinxext_plot_directive_api.html
plot_formats = ["png"]
plot_include_source = True
plot_html_show_formats = False
plot_html_show_source_link = False

# this is needed for some reason...
# see https://github.com/numpy/numpydoc/issues/69
numpydoc_class_members_toctree = False


# For maths, use mathjax by default and svg if NO_MATHJAX env variable is set
# (useful for viewing the doc offline)
if os.environ.get("NO_MATHJAX"):
    extensions.append("sphinx.ext.imgmath")
    imgmath_image_format = "svg"
    mathjax_path = ""
else:
    extensions.append("sphinx.ext.mathjax")
    mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"

autodoc_default_options = {"members": True, "inherited-members": True}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["templates"]

# generate autosummary even if no references
autosummary_generate = True

# The suffix of source filenames.
source_suffix = ".rst"

# The encoding of source files.
# source_encoding = 'utf-8'

# The main toctree document.
root_doc = "contents"

# General information about the project.
project = "scikit-learn"
copyright = f"2007 - {datetime.now().year}, scikit-learn developers (BSD License)"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
import sklearn

parsed_version = parse(sklearn.__version__)
version = ".".join(parsed_version.base_version.split(".")[:2])
# The full version, including alpha/beta/rc tags.
# Removes post from release name
if parsed_version.is_postrelease:
    release = parsed_version.base_version
else:
    release = sklearn.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
# language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
# today = ''
# Else, today_fmt is used as the format for a strftime call.
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "templates", "includes", "themes"]

# The reST default role (used for this markup: `text`) to use for all
# documents.
default_role = "literal"

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  Major themes that come with
# Sphinx are currently 'default' and 'sphinxdoc'.
html_theme = "scikit-learn-modern"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "legacy_google_analytics": True,
    "analytics": True,
    "mathjax_path": mathjax_path,
    "link_to_live_contributing_page": not parsed_version.is_devrelease,
}

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = ["themes"]


# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = "scikit-learn"

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "logos/scikit-learn-logo-small.png"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "logos/favicon.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["images"]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
# html_last_updated_fmt = '%b %d, %Y'

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
html_additional_pages = {"index": "index.html"}

# If false, no module index is generated.
html_domain_indices = False

# If false, no index is generated.
html_use_index = False

# If true, the index is split into individual pages for each letter.
# html_split_index = False

# If true, links to the reST sources are added to the pages.
# html_show_sourcelink = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
# html_use_opensearch = ''

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = ''

# Output file base name for HTML help builder.
htmlhelp_basename = "scikit-learndoc"

# If true, the reST sources are included in the HTML build as _sources/name.
html_copy_source = True

# Adds variables into templates
html_context = {}
# finds latest release highlights and places it into HTML context for
# index.html
release_highlights_dir = Path("..") / "examples" / "release_highlights"
# Finds the highlight with the latest version number
latest_highlights = sorted(release_highlights_dir.glob("plot_release_highlights_*.py"))[
    -1
]
latest_highlights = latest_highlights.with_suffix("").name
html_context["release_highlights"] = (
    f"auto_examples/release_highlights/{latest_highlights}"
)

# get version from highlight name assuming highlights have the form
# plot_release_highlights_0_22_0
highlight_version = ".".join(latest_highlights.split("_")[-3:-1])
html_context["release_highlights_version"] = highlight_version


# redirects dictionary maps from old links to new links
redirects = {
    "documentation": "index",
    "auto_examples/feature_selection/plot_permutation_test_for_classification": (
        "auto_examples/model_selection/plot_permutation_tests_for_classification"
    ),
    "modules/model_persistence": "model_persistence",
    "auto_examples/linear_model/plot_bayesian_ridge": (
        "auto_examples/linear_model/plot_ard"
    ),
    "auto_examples/model_selection/grid_search_text_feature_extraction.py": (
        "auto_examples/model_selection/plot_grid_search_text_feature_extraction.py"
    ),
    "auto_examples/miscellaneous/plot_changed_only_pprint_parameter": (
        "auto_examples/miscellaneous/plot_estimator_representation"
    ),
    "auto_examples/decomposition/plot_beta_divergence": (
        "auto_examples/applications/plot_topics_extraction_with_nmf_lda"
    ),
    "auto_examples/svm/plot_svm_nonlinear": "auto_examples/svm/plot_svm_kernels",
    "auto_examples/ensemble/plot_adaboost_hastie_10_2": (
        "auto_examples/ensemble/plot_adaboost_multiclass"
    ),
}
html_context["redirects"] = redirects
for old_link in redirects:
    html_additional_pages[old_link] = "redirects.html"

# Not showing the search summary makes the search page load faster.
html_show_search_summary = False


# The "summary-anchor" IDs will be overwritten via JavaScript to be unique.
# See `doc/theme/scikit-learn-modern/static/js/details-permalink.js`.
rst_prolog = """
.. |details-start| raw:: html

    <details id="summary-anchor">
    <summary class="btn btn-light">

.. |details-split| raw:: html

    <span class="tooltiptext">Click for more details</span>
    <a class="headerlink" href="#summary-anchor" title="Permalink to this heading">¶</a>
    </summary>
    <div class="card">

.. |details-end| raw:: html

    </div>
    </details>

"""

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    "preamble": r"""
        \usepackage{amsmath}\usepackage{amsfonts}\usepackage{bm}
        \usepackage{morefloats}\usepackage{enumitem} \setlistdepth{10}
        \let\oldhref\href
        \renewcommand{\href}[2]{\oldhref{#1}{\hbox{#2}}}
        """
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass
# [howto/manual]).
latex_documents = [
    (
        "contents",
        "user_guide.tex",
        "scikit-learn user guide",
        "scikit-learn developers",
        "manual",
    ),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
latex_logo = "logos/scikit-learn-logo.png"

# Documents to append as an appendix to all manuals.
# latex_appendices = []

# If false, no module index is generated.
latex_domain_indices = False

trim_doctests_flags = True

# intersphinx configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/{.major}".format(sys.version_info), None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest/", None),
    "seaborn": ("https://seaborn.pydata.org/", None),
    "skops": ("https://skops.readthedocs.io/en/stable/", None),
}

v = parse(release)
if v.release is None:
    raise ValueError(
        "Ill-formed version: {!r}. Version should follow PEP440".format(version)
    )

if v.is_devrelease:
    binder_branch = "main"
else:
    major, minor = v.release[:2]
    binder_branch = "{}.{}.X".format(major, minor)


class SubSectionTitleOrder:
    """Sort example gallery by title of subsection.

    Assumes README.txt exists for all subsections and uses the subsection with
    dashes, '---', as the adornment.
    """

    def __init__(self, src_dir):
        self.src_dir = src_dir
        self.regex = re.compile(r"^([\w ]+)\n-", re.MULTILINE)

    def __repr__(self):
        return "<%s>" % (self.__class__.__name__,)

    def __call__(self, directory):
        src_path = os.path.normpath(os.path.join(self.src_dir, directory))

        # Forces Release Highlights to the top
        if os.path.basename(src_path) == "release_highlights":
            return "0"

        readme = os.path.join(src_path, "README.txt")

        try:
            with open(readme, "r") as f:
                content = f.read()
        except FileNotFoundError:
            return directory

        title_match = self.regex.search(content)
        if title_match is not None:
            return title_match.group(1)
        return directory


class SKExampleTitleSortKey(ExampleTitleSortKey):
    """Sorts release highlights based on version number."""

    def __call__(self, filename):
        title = super().__call__(filename)
        prefix = "plot_release_highlights_"

        # Use title to sort if not a release highlight
        if not filename.startswith(prefix):
            return title

        major_minor = filename[len(prefix) :].split("_")[:2]
        version_float = float(".".join(major_minor))

        # negate to place the newest version highlights first
        return -version_float


def notebook_modification_function(notebook_content, notebook_filename):
    notebook_content_str = str(notebook_content)
    warning_template = "\n".join(
        [
            "<div class='alert alert-{message_class}'>",
            "",
            "# JupyterLite warning",
            "",
            "{message}",
            "</div>",
        ]
    )

    message_class = "warning"
    message = (
        "Running the scikit-learn examples in JupyterLite is experimental and you may"
        " encounter some unexpected behavior.\n\nThe main difference is that imports"
        " will take a lot longer than usual, for example the first `import sklearn` can"
        " take roughly 10-20s.\n\nIf you notice problems, feel free to open an"
        " [issue](https://github.com/scikit-learn/scikit-learn/issues/new/choose)"
        " about it."
    )

    markdown = warning_template.format(message_class=message_class, message=message)

    dummy_notebook_content = {"cells": []}
    add_markdown_cell(dummy_notebook_content, markdown)

    code_lines = []

    if "seaborn" in notebook_content_str:
        code_lines.append("%pip install seaborn")
    if "plotly.express" in notebook_content_str:
        code_lines.append("%pip install plotly")
    if "skimage" in notebook_content_str:
        code_lines.append("%pip install scikit-image")
    if "fetch_" in notebook_content_str:
        code_lines.extend(
            [
                "%pip install pyodide-http",
                "import pyodide_http",
                "pyodide_http.patch_all()",
            ]
        )
    # always import matplotlib and pandas to avoid Pyodide limitation with
    # imports inside functions
    code_lines.extend(["import matplotlib", "import pandas"])

    if code_lines:
        code_lines = ["# JupyterLite-specific code"] + code_lines
        code = "\n".join(code_lines)
        add_code_cell(dummy_notebook_content, code)

    notebook_content["cells"] = (
        dummy_notebook_content["cells"] + notebook_content["cells"]
    )


default_global_config = sklearn.get_config()


def reset_sklearn_config(gallery_conf, fname):
    """Reset sklearn config to default values."""
    sklearn.set_config(**default_global_config)


sphinx_gallery_conf = {
    "doc_module": "sklearn",
    "backreferences_dir": os.path.join("modules", "generated"),
    "show_memory": False,
    "reference_url": {"sklearn": None},
    "examples_dirs": ["../examples"],
    "gallery_dirs": ["auto_examples"],
    "subsection_order": SubSectionTitleOrder("../examples"),
    "within_subsection_order": SKExampleTitleSortKey,
    "binder": {
        "org": "scikit-learn",
        "repo": "scikit-learn",
        "binderhub_url": "https://mybinder.org",
        "branch": binder_branch,
        "dependencies": "./binder/requirements.txt",
        "use_jupyter_lab": True,
    },
    # avoid generating too many cross links
    "inspect_global_variables": False,
    "remove_config_comments": True,
    "plot_gallery": "True",
    "reset_modules": ("matplotlib", "seaborn", reset_sklearn_config),
}
if with_jupyterlite:
    sphinx_gallery_conf["jupyterlite"] = {
        "notebook_modification_function": notebook_modification_function
    }


# The following dictionary contains the information used to create the
# thumbnails for the front page of the scikit-learn home page.
# key: first image in set
# values: (number of plot in set, height of thumbnail)
carousel_thumbs = {"sphx_glr_plot_classifier_comparison_001.png": 600}


# enable experimental module so that experimental estimators can be
# discovered properly by sphinx
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.experimental import enable_halving_search_cv  # noqa


def make_carousel_thumbs(app, exception):
    """produces the final resized carousel images"""
    if exception is not None:
        return
    print("Preparing carousel images")

    image_dir = os.path.join(app.builder.outdir, "_images")
    for glr_plot, max_width in carousel_thumbs.items():
        image = os.path.join(image_dir, glr_plot)
        if os.path.exists(image):
            c_thumb = os.path.join(image_dir, glr_plot[:-4] + "_carousel.png")
            sphinx_gallery.gen_rst.scale_image(image, c_thumb, max_width, 190)


def filter_search_index(app, exception):
    if exception is not None:
        return

    # searchindex only exist when generating html
    if app.builder.name != "html":
        return

    print("Removing methods from search index")

    searchindex_path = os.path.join(app.builder.outdir, "searchindex.js")
    with open(searchindex_path, "r") as f:
        searchindex_text = f.read()

    searchindex_text = re.sub(r"{__init__.+?}", "{}", searchindex_text)
    searchindex_text = re.sub(r"{__call__.+?}", "{}", searchindex_text)

    with open(searchindex_path, "w") as f:
        f.write(searchindex_text)


def generate_min_dependency_table(app):
    """Generate min dependency table for docs."""
    from sklearn._min_dependencies import dependent_packages

    # get length of header
    package_header_len = max(len(package) for package in dependent_packages) + 4
    version_header_len = len("Minimum Version") + 4
    tags_header_len = max(len(tags) for _, tags in dependent_packages.values()) + 4

    output = StringIO()
    output.write(
        " ".join(
            ["=" * package_header_len, "=" * version_header_len, "=" * tags_header_len]
        )
    )
    output.write("\n")
    dependency_title = "Dependency"
    version_title = "Minimum Version"
    tags_title = "Purpose"

    output.write(
        f"{dependency_title:<{package_header_len}} "
        f"{version_title:<{version_header_len}} "
        f"{tags_title}\n"
    )

    output.write(
        " ".join(
            ["=" * package_header_len, "=" * version_header_len, "=" * tags_header_len]
        )
    )
    output.write("\n")

    for package, (version, tags) in dependent_packages.items():
        output.write(
            f"{package:<{package_header_len}} {version:<{version_header_len}} {tags}\n"
        )

    output.write(
        " ".join(
            ["=" * package_header_len, "=" * version_header_len, "=" * tags_header_len]
        )
    )
    output.write("\n")
    output = output.getvalue()

    with (Path(".") / "min_dependency_table.rst").open("w") as f:
        f.write(output)


def generate_min_dependency_substitutions(app):
    """Generate min dependency substitutions for docs."""
    from sklearn._min_dependencies import dependent_packages

    output = StringIO()

    for package, (version, _) in dependent_packages.items():
        package = package.capitalize()
        output.write(f".. |{package}MinVersion| replace:: {version}")
        output.write("\n")

    output = output.getvalue()

    with (Path(".") / "min_dependency_substitutions.rst").open("w") as f:
        f.write(output)


# Config for sphinx_issues

# we use the issues path for PRs since the issues URL will forward
issues_github_path = "scikit-learn/scikit-learn"


def disable_plot_gallery_for_linkcheck(app):
    if app.builder.name == "linkcheck":
        sphinx_gallery_conf["plot_gallery"] = "False"


def setup(app):
    # do not run the examples when using linkcheck by using a small priority
    # (default priority is 500 and sphinx-gallery using builder-inited event too)
    app.connect("builder-inited", disable_plot_gallery_for_linkcheck, priority=50)
    app.connect("builder-inited", generate_min_dependency_table)
    app.connect("builder-inited", generate_min_dependency_substitutions)

    # to hide/show the prompt in code examples:
    app.connect("build-finished", make_carousel_thumbs)
    app.connect("build-finished", filter_search_index)


# The following is used by sphinx.ext.linkcode to provide links to github
linkcode_resolve = make_linkcode_resolve(
    "sklearn",
    (
        "https://github.com/scikit-learn/"
        "scikit-learn/blob/{revision}/"
        "{package}/{path}#L{lineno}"
    ),
)

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=(
        "Matplotlib is currently using agg, which is a"
        " non-GUI backend, so cannot show the figure."
    ),
)

# maps functions with a class name that is indistinguishable when case is
# ignore to another filename
autosummary_filename_map = {
    "sklearn.cluster.dbscan": "dbscan-function",
    "sklearn.covariance.oas": "oas-function",
    "sklearn.decomposition.fastica": "fastica-function",
}


# Config for sphinxext.opengraph

ogp_site_url = "https://scikit-learn/stable/"
ogp_image = "https://scikit-learn.org/stable/_static/scikit-learn-logo-small.png"
ogp_use_first_image = True
ogp_site_name = "scikit-learn"

# Config for linkcheck that checks the documentation for broken links

# ignore all links in 'whats_new' to avoid doing many github requests and
# hitting the github rate threshold that makes linkcheck take a lot of time
linkcheck_exclude_documents = [r"whats_new/.*"]

# default timeout to make some sites links fail faster
linkcheck_timeout = 10

# Allow redirects from doi.org
linkcheck_allowed_redirects = {r"https://doi.org/.+": r".*"}
linkcheck_ignore = [
    # ignore links to local html files e.g. in image directive :target: field
    r"^..?/",
    # ignore links to specific pdf pages because linkcheck does not handle them
    # ('utf-8' codec can't decode byte error)
    r"http://www.utstat.toronto.edu/~rsalakhu/sta4273/notes/Lecture2.pdf#page=.*",
    (
        "https://www.fordfoundation.org/media/2976/roads-and-bridges"
        "-the-unseen-labor-behind-our-digital-infrastructure.pdf#page=.*"
    ),
    # links falsely flagged as broken
    (
        "https://www.researchgate.net/publication/"
        "233096619_A_Dendrite_Method_for_Cluster_Analysis"
    ),
    (
        "https://www.researchgate.net/publication/221114584_Random_Fourier"
        "_Approximations_for_Skewed_Multiplicative_Histogram_Kernels"
    ),
    (
        "https://www.researchgate.net/publication/4974606_"
        "Hedonic_housing_prices_and_the_demand_for_clean_air"
    ),
    (
        "https://www.researchgate.net/profile/Anh-Huy-Phan/publication/220241471_Fast_"
        "Local_Algorithms_for_Large_Scale_Nonnegative_Matrix_and_Tensor_Factorizations"
    ),
    "https://doi.org/10.13140/RG.2.2.35280.02565",
    (
        "https://www.microsoft.com/en-us/research/uploads/prod/2006/01/"
        "Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf"
    ),
    "https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/tr-99-87.pdf",
    "https://microsoft.com/",
    "https://www.jstor.org/stable/2984099",
    "https://stat.uw.edu/sites/default/files/files/reports/2000/tr371.pdf",
    # Broken links from testimonials
    "http://www.bestofmedia.com",
    "http://www.data-publica.com/",
    "https://livelovely.com",
    "https://www.mars.com/global",
    "https://www.yhat.com",
    # Ignore some dynamically created anchors. See
    # https://github.com/sphinx-doc/sphinx/issues/9016 for more details about
    # the github example
    r"https://github.com/conda-forge/miniforge#miniforge",
    r"https://github.com/joblib/threadpoolctl/"
    "#setting-the-maximum-size-of-thread-pools",
    r"https://stackoverflow.com/questions/5836335/"
    "consistently-create-same-random-numpy-array/5837352#comment6712034_5837352",
]

# Use a browser-like user agent to avoid some "403 Client Error: Forbidden for
# url" errors. This is taken from the variable navigator.userAgent inside a
# browser console.
user_agent = (
    "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:100.0) Gecko/20100101 Firefox/100.0"
)

# Use Github token from environment variable to avoid Github rate limits when
# checking Github links
github_token = os.getenv("GITHUB_TOKEN")

if github_token is None:
    linkcheck_request_headers = {}
else:
    linkcheck_request_headers = {
        "https://github.com/": {"Authorization": f"token {github_token}"},
    }
