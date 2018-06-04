# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
Sphinx-Gallery Generator
========================

Attaches Sphinx-Gallery to Sphinx in order to generate the galleries
when building the documentation.
"""


from __future__ import division, print_function, absolute_import
import codecs
import copy
import re
import os

from . import sphinx_compatibility, glr_path_static, __version__ as _sg_version
from .gen_rst import generate_dir_rst, SPHX_GLR_SIG
from .docs_resolv import embed_code_links
from .downloads import generate_zipfiles
from .sorting import NumberOfCodeLinesSortKey

try:
    FileNotFoundError
except NameError:
    # Python2
    FileNotFoundError = IOError

DEFAULT_GALLERY_CONF = {
    'filename_pattern': re.escape(os.sep) + 'plot',
    'examples_dirs': os.path.join('..', 'examples'),
    'subsection_order': None,
    'within_subsection_order': NumberOfCodeLinesSortKey,
    'gallery_dirs': 'auto_examples',
    'backreferences_dir': None,
    'doc_module': (),
    'reference_url': {},
    # Build options
    # -------------
    # We use a string for 'plot_gallery' rather than simply the Python boolean
    # `True` as it avoids a warning about unicode when controlling this value
    # via the command line switches of sphinx-build
    'plot_gallery': 'True',
    'download_all_examples': True,
    'abort_on_example_error': False,
    'failing_examples': {},
    'expected_failing_examples': set(),
    'thumbnail_size': (400, 280),  # Default CSS does 0.4 scaling (160, 112)
    'min_reported_time': 0,
}

logger = sphinx_compatibility.getLogger('sphinx-gallery')


def clean_gallery_out(build_dir):
    """Deletes images under the sphx_glr namespace in the build directory"""
    # Sphinx hack: sphinx copies generated images to the build directory
    #  each time the docs are made.  If the desired image name already
    #  exists, it appends a digit to prevent overwrites.  The problem is,
    #  the directory is never cleared.  This means that each time you build
    #  the docs, the number of images in the directory grows.
    #
    # This question has been asked on the sphinx development list, but there
    #  was no response: https://git.net/ml/sphinx-dev/2011-02/msg00123.html
    #
    # The following is a hack that prevents this behavior by clearing the
    #  image build directory from gallery images each time the docs are built.
    #  If sphinx changes their layout between versions, this will not
    #  work (though it should probably not cause a crash).
    # Tested successfully on Sphinx 1.0.7

    build_image_dir = os.path.join(build_dir, '_images')
    if os.path.exists(build_image_dir):
        filelist = os.listdir(build_image_dir)
        for filename in filelist:
            if filename.startswith('sphx_glr') and filename.endswith('png'):
                os.remove(os.path.join(build_image_dir, filename))


def parse_config(app):
    """Process the Sphinx Gallery configuration"""
    try:
        plot_gallery = eval(app.builder.config.plot_gallery)
    except TypeError:
        plot_gallery = bool(app.builder.config.plot_gallery)

    gallery_conf = copy.deepcopy(DEFAULT_GALLERY_CONF)
    gallery_conf.update(app.config.sphinx_gallery_conf)
    gallery_conf.update(plot_gallery=plot_gallery)
    gallery_conf.update(
        abort_on_example_error=app.builder.config.abort_on_example_error)
    gallery_conf['src_dir'] = app.builder.srcdir

    backreferences_warning = """\n========
Sphinx-Gallery now requires you to set the configuration variable
'backreferences_dir' in your config to activate the
backreferences. That is mini galleries clustered by the functions used
in the example scripts. Have a look at it in sphinx-gallery

https://sphinx-gallery.readthedocs.io/en/stable/index.html#examples-using-numpy-linspace
"""

    if gallery_conf.get("mod_example_dir", False):
        update_msg = """\nFor a quick fix try replacing 'mod_example_dir'
by 'backreferences_dir' in your conf.py file. If that does not solve the
present issue read carefully how to update in the online documentation

https://sphinx-gallery.readthedocs.io/en/latest/advanced_configuration.html#references-to-examples"""

        gallery_conf['backreferences_dir'] = gallery_conf['mod_example_dir']
        logger.warning(
            "Old configuration for backreferences detected \n"
            "using the configuration variable `mod_example_dir`\n"
            "%s%s",
            backreferences_warning,
            update_msg,
            type=DeprecationWarning)

    elif gallery_conf['backreferences_dir'] is None:
        no_care_msg = """
If you don't care about this features set in your conf.py
'backreferences_dir': False\n"""

        logger.warning(backreferences_warning + no_care_msg)

        gallery_conf['backreferences_dir'] = os.path.join(
            'modules', 'generated')
        logger.warning(
            "Using old default 'backreferences_dir':'%s'.\n"
            "This will be disabled in future releases\n",
            gallery_conf['backreferences_dir'],
            type=DeprecationWarning)

    # this assures I can call the config in other places
    app.config.sphinx_gallery_conf = gallery_conf
    app.config.html_static_path.append(glr_path_static())

    return gallery_conf


def get_subsections(srcdir, examples_dir, sortkey):
    """Returns the list of subsections of a gallery

    Parameters
    ----------
    srcdir : str
        absolute path to directory containing conf.py
    examples_dir : str
        path to the examples directory relative to conf.py
    sortkey : callable
        The sort key to use.

    Returns
    -------
    out : list
        sorted list of gallery subsection folder names

    """
    subfolders = [subfolder for subfolder in os.listdir(examples_dir)
                  if os.path.exists(os.path.join(examples_dir, subfolder, 'README.txt'))]
    base_examples_dir_path = os.path.relpath(examples_dir, srcdir)
    subfolders_with_path = [os.path.join(base_examples_dir_path, item)
                            for item in subfolders]
    sorted_subfolders = sorted(subfolders_with_path, key=sortkey)

    return [subfolders[i] for i in [subfolders_with_path.index(item)
                                    for item in sorted_subfolders]]


def _prepare_sphx_glr_dirs(gallery_conf, srcdir):
    """Creates necessary folders for sphinx_gallery files """
    examples_dirs = gallery_conf['examples_dirs']
    gallery_dirs = gallery_conf['gallery_dirs']

    if not isinstance(examples_dirs, list):
        examples_dirs = [examples_dirs]

    if not isinstance(gallery_dirs, list):
        gallery_dirs = [gallery_dirs]

    for outdir in gallery_dirs:
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    if bool(gallery_conf['backreferences_dir']):
        backreferences_dir = os.path.join(
            srcdir, gallery_conf['backreferences_dir'])
        if not os.path.exists(backreferences_dir):
            os.makedirs(backreferences_dir)

    return zip(examples_dirs, gallery_dirs)


def generate_gallery_rst(app):
    """Generate the Main examples gallery reStructuredText

    Start the sphinx-gallery configuration and recursively scan the examples
    directories in order to populate the examples gallery
    """
    logger.info('generating gallery...', color='white')
    gallery_conf = parse_config(app)

    clean_gallery_out(app.builder.outdir)

    seen_backrefs = set()

    computation_times = []
    workdirs = _prepare_sphx_glr_dirs(gallery_conf,
                                      app.builder.srcdir)

    for examples_dir, gallery_dir in workdirs:

        examples_dir = os.path.join(app.builder.srcdir, examples_dir)
        gallery_dir = os.path.join(app.builder.srcdir, gallery_dir)

        if not os.path.exists(os.path.join(examples_dir, 'README.txt')):
            raise FileNotFoundError("Main example directory {0} does not "
                                    "have a README.txt file. Please write "
                                    "one to introduce your gallery."
                                    .format(examples_dir))

        # Here we don't use an os.walk, but we recurse only twice: flat is
        # better than nested.

        this_fhindex, this_computation_times = generate_dir_rst(
            examples_dir, gallery_dir, gallery_conf, seen_backrefs)

        computation_times += this_computation_times

        # we create an index.rst with all examples
        with codecs.open(os.path.join(gallery_dir, 'index.rst'), 'w',
                         encoding='utf-8') as fhindex:
            # :orphan: to suppress "not included in TOCTREE" sphinx warnings
            fhindex.write(":orphan:\n\n" + this_fhindex)

            for subsection in get_subsections(app.builder.srcdir, examples_dir, gallery_conf['subsection_order']):
                src_dir = os.path.join(examples_dir, subsection)
                target_dir = os.path.join(gallery_dir, subsection)
                this_fhindex, this_computation_times = generate_dir_rst(src_dir, target_dir, gallery_conf,
                                                                        seen_backrefs)
                fhindex.write(this_fhindex)
                computation_times += this_computation_times

            if gallery_conf['download_all_examples']:
                download_fhindex = generate_zipfiles(gallery_dir)
                fhindex.write(download_fhindex)

            fhindex.write(SPHX_GLR_SIG)

    if gallery_conf['plot_gallery']:
        logger.info("computation time summary:", color='white')
        for time_elapsed, fname in sorted(computation_times, reverse=True):
            if time_elapsed is not None:
                if time_elapsed >= gallery_conf['min_reported_time']:
                    logger.info("\t- %s: %.2g sec", fname, time_elapsed)
            else:
                logger.info("\t- %s: not run", fname)


def touch_empty_backreferences(app, what, name, obj, options, lines):
    """Generate empty back-reference example files

    This avoids inclusion errors/warnings if there are no gallery
    examples for a class / module that is being parsed by autodoc"""

    if not bool(app.config.sphinx_gallery_conf['backreferences_dir']):
        return

    examples_path = os.path.join(app.srcdir,
                                 app.config.sphinx_gallery_conf[
                                     "backreferences_dir"],
                                 "%s.examples" % name)

    if not os.path.exists(examples_path):
        # touch file
        open(examples_path, 'w').close()


def sumarize_failing_examples(app, exception):
    """Collects the list of falling examples during build and prints them with the traceback

    Raises ValueError if there where failing examples
    """
    if exception is not None:
        return

    # Under no-plot Examples are not run so nothing to summarize
    if not app.config.sphinx_gallery_conf['plot_gallery']:
        return

    gallery_conf = app.config.sphinx_gallery_conf
    failing_examples = set(gallery_conf['failing_examples'].keys())
    expected_failing_examples = set([os.path.normpath(os.path.join(app.srcdir, path))
                                     for path in
                                     gallery_conf['expected_failing_examples']])

    examples_expected_to_fail = failing_examples.intersection(
        expected_failing_examples)
    if examples_expected_to_fail:
        logger.info("Examples failing as expected:", color='brown')
        for fail_example in examples_expected_to_fail:
            logger.info('%s failed leaving traceback:', fail_example)
            logger.info(gallery_conf['failing_examples'][fail_example])

    examples_not_expected_to_fail = failing_examples.difference(
        expected_failing_examples)
    fail_msgs = []
    if examples_not_expected_to_fail:
        fail_msgs.append("Unexpected failing examples:")
        for fail_example in examples_not_expected_to_fail:
            fail_msgs.append(fail_example + ' failed leaving traceback:\n' +
                             gallery_conf['failing_examples'][fail_example] +
                             '\n')

    examples_not_expected_to_pass = expected_failing_examples.difference(
        failing_examples)
    if examples_not_expected_to_pass:
        fail_msgs.append("Examples expected to fail, but not failling:\n" +
                         "Please remove these examples from\n" +
                         "sphinx_gallery_conf['expected_failing_examples']\n" +
                         "in your conf.py file"
                         "\n".join(examples_not_expected_to_pass))

    if fail_msgs:
        raise ValueError("Here is a summary of the problems encountered when "
                         "running the examples\n\n" + "\n".join(fail_msgs) +
                         "\n" + "-" * 79)


def get_default_config_value(key):
    def default_getter(conf):
        return conf['sphinx_gallery_conf'].get(key, DEFAULT_GALLERY_CONF[key])
    return default_getter


def setup(app):
    """Setup sphinx-gallery sphinx extension"""
    sphinx_compatibility._app = app

    app.add_config_value('sphinx_gallery_conf', DEFAULT_GALLERY_CONF, 'html')
    for key in ['plot_gallery', 'abort_on_example_error']:
        app.add_config_value(key, get_default_config_value(key), 'html')

    app.add_stylesheet('gallery.css')

    # Sphinx < 1.6 calls it `_extensions`, >= 1.6 is `extensions`.
    extensions_attr = '_extensions' if hasattr(
        app, '_extensions') else 'extensions'
    if 'sphinx.ext.autodoc' in getattr(app, extensions_attr):
        app.connect('autodoc-process-docstring', touch_empty_backreferences)

    app.connect('builder-inited', generate_gallery_rst)

    app.connect('build-finished', sumarize_failing_examples)
    app.connect('build-finished', embed_code_links)
    metadata = {'parallel_read_safe': True,
                'version': _sg_version}
    return metadata


def setup_module():
    # HACK: Stop nosetests running setup() above
    pass
