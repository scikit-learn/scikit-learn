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
from datetime import timedelta, datetime
from difflib import get_close_matches
from importlib import import_module
import re
import os
import pathlib
from xml.sax.saxutils import quoteattr, escape

from sphinx.errors import ConfigError, ExtensionError
from sphinx.util.console import red
from . import sphinx_compatibility, glr_path_static, __version__ as _sg_version
from .utils import _replace_md5, _has_optipng, _has_pypandoc
from .backreferences import _finalize_backreferences
from .gen_rst import (generate_dir_rst, SPHX_GLR_SIG, _get_memory_base,
                      _get_readme)
from .scrapers import _scraper_dict, _reset_dict, _import_matplotlib
from .docs_resolv import embed_code_links
from .downloads import generate_zipfiles
from .sorting import NumberOfCodeLinesSortKey
from .binder import copy_binder_files, check_binder_conf
from .directives import MiniGallery, ImageSg, imagesg_addnode


_KNOWN_CSS = ('sg_gallery', 'sg_gallery-binder', 'sg_gallery-dataframe',
              'sg_gallery-rendered-html')


class DefaultResetArgv:
    def __repr__(self):
        return "DefaultResetArgv"

    def __call__(self, gallery_conf, script_vars):
        return []


DEFAULT_GALLERY_CONF = {
    'filename_pattern': re.escape(os.sep) + 'plot',
    'ignore_pattern': r'__init__\.py',
    'examples_dirs': os.path.join('..', 'examples'),
    'reset_argv': DefaultResetArgv(),
    'subsection_order': None,
    'within_subsection_order': NumberOfCodeLinesSortKey,
    'gallery_dirs': 'auto_examples',
    'backreferences_dir': None,
    'doc_module': (),
    'reference_url': {},
    'capture_repr': ('_repr_html_', '__repr__'),
    'ignore_repr_types': r'',
    # Build options
    # -------------
    # 'plot_gallery' also accepts strings that evaluate to a bool, e.g. "True",
    # "False", "1", "0" so that they can be easily set via command line
    # switches of sphinx-build
    'plot_gallery': True,
    'download_all_examples': True,
    'abort_on_example_error': False,
    'only_warn_on_example_error': False,
    'failing_examples': {},
    'passing_examples': [],
    'stale_examples': [],  # ones that did not need to be run due to md5sum
    'run_stale_examples': False,
    'expected_failing_examples': set(),
    'thumbnail_size': (400, 280),  # Default CSS does 0.4 scaling (160, 112)
    'min_reported_time': 0,
    'binder': {},
    'image_scrapers': ('matplotlib',),
    'compress_images': (),
    'reset_modules': ('matplotlib', 'seaborn'),
    'reset_modules_order': 'before',
    'first_notebook_cell': '%matplotlib inline',
    'last_notebook_cell': None,
    'notebook_images': False,
    'pypandoc': False,
    'remove_config_comments': False,
    'show_memory': False,
    'show_signature': True,
    'junit': '',
    'log_level': {'backreference_missing': 'warning'},
    'inspect_global_variables': True,
    'css': _KNOWN_CSS,
    'matplotlib_animations': False,
    'image_srcset': [],
    'default_thumb_file': None,
    'line_numbers': False,
}

logger = sphinx_compatibility.getLogger('sphinx-gallery')


def _bool_eval(x):
    if isinstance(x, str):
        try:
            x = eval(x)
        except TypeError:
            pass
    return bool(x)


def parse_config(app, check_keys=True):
    """Process the Sphinx Gallery configuration."""
    plot_gallery = _bool_eval(app.builder.config.plot_gallery)
    src_dir = app.builder.srcdir
    abort_on_example_error = _bool_eval(
        app.builder.config.abort_on_example_error)
    lang = app.builder.config.highlight_language
    gallery_conf = _complete_gallery_conf(
        app.config.sphinx_gallery_conf, src_dir, plot_gallery,
        abort_on_example_error, lang, app.builder.name, app,
        check_keys)

    # this assures I can call the config in other places
    app.config.sphinx_gallery_conf = gallery_conf
    app.config.html_static_path.append(glr_path_static())
    return gallery_conf


def _complete_gallery_conf(sphinx_gallery_conf, src_dir, plot_gallery,
                           abort_on_example_error, lang='python',
                           builder_name='html', app=None, check_keys=True):
    gallery_conf = copy.deepcopy(DEFAULT_GALLERY_CONF)
    options = sorted(gallery_conf)
    extra_keys = sorted(set(sphinx_gallery_conf) - set(options))
    if extra_keys and check_keys:
        msg = 'Unknown key(s) in sphinx_gallery_conf:\n'
        for key in extra_keys:
            options = get_close_matches(key, options, cutoff=0.66)
            msg += repr(key)
            if len(options) == 1:
                msg += ', did you mean %r?' % (options[0],)
            elif len(options) > 1:
                msg += ', did you mean one of %r?' % (options,)
            msg += '\n'
        raise ConfigError(msg.strip())
    gallery_conf.update(sphinx_gallery_conf)
    if sphinx_gallery_conf.get('find_mayavi_figures', False):
        logger.warning(
            "Deprecated image scraping variable `find_mayavi_figures`\n"
            "detected, use `image_scrapers` instead as:\n\n"
            "   image_scrapers=('matplotlib', 'mayavi')",
            type=DeprecationWarning)
        gallery_conf['image_scrapers'] += ('mayavi',)
    gallery_conf.update(plot_gallery=plot_gallery)
    gallery_conf.update(abort_on_example_error=abort_on_example_error)
    # XXX anything that can only be a bool (rather than str) should probably be
    # evaluated this way as it allows setting via -D on the command line
    for key in ('run_stale_examples',):
        gallery_conf[key] = _bool_eval(gallery_conf[key])
    gallery_conf['src_dir'] = src_dir
    gallery_conf['app'] = app

    # Check capture_repr
    capture_repr = gallery_conf['capture_repr']
    supported_reprs = ['__repr__', '__str__', '_repr_html_']
    if isinstance(capture_repr, tuple):
        for rep in capture_repr:
            if rep not in supported_reprs:
                raise ConfigError("All entries in 'capture_repr' must be one "
                                  "of %s, got: %s" % (supported_reprs, rep))
    else:
        raise ConfigError("'capture_repr' must be a tuple, got: %s"
                          % (type(capture_repr),))
    # Check ignore_repr_types
    if not isinstance(gallery_conf['ignore_repr_types'], str):
        raise ConfigError("'ignore_repr_types' must be a string, got: %s"
                          % (type(gallery_conf['ignore_repr_types']),))

    # deal with show_memory
    gallery_conf['memory_base'] = 0.
    if gallery_conf['show_memory']:
        if not callable(gallery_conf['show_memory']):  # True-like
            try:
                from memory_profiler import memory_usage  # noqa
            except ImportError:
                logger.warning("Please install 'memory_profiler' to enable "
                               "peak memory measurements.")
                gallery_conf['show_memory'] = False
            else:
                def call_memory(func):
                    mem, out = memory_usage(func, max_usage=True, retval=True,
                                            multiprocess=True)
                    try:
                        mem = mem[0]  # old MP always returned a list
                    except TypeError:  # 'float' object is not subscriptable
                        pass
                    return mem, out
                gallery_conf['call_memory'] = call_memory
                gallery_conf['memory_base'] = _get_memory_base(gallery_conf)
        else:
            gallery_conf['call_memory'] = gallery_conf['show_memory']
    if not gallery_conf['show_memory']:  # can be set to False above
        def call_memory(func):
            return 0., func()
        gallery_conf['call_memory'] = call_memory
    assert callable(gallery_conf['call_memory'])

    # deal with scrapers
    scrapers = gallery_conf['image_scrapers']
    if not isinstance(scrapers, (tuple, list)):
        scrapers = [scrapers]
    scrapers = list(scrapers)
    for si, scraper in enumerate(scrapers):
        if isinstance(scraper, str):
            if scraper in _scraper_dict:
                scraper = _scraper_dict[scraper]
            else:
                orig_scraper = scraper
                try:
                    scraper = import_module(scraper)
                    scraper = getattr(scraper, '_get_sg_image_scraper')
                    scraper = scraper()
                except Exception as exp:
                    raise ConfigError('Unknown image scraper %r, got:\n%s'
                                      % (orig_scraper, exp))
            scrapers[si] = scraper
        if not callable(scraper):
            raise ConfigError('Scraper %r was not callable' % (scraper,))
    gallery_conf['image_scrapers'] = tuple(scrapers)
    del scrapers
    # Here we try to set up matplotlib but don't raise an error,
    # we will raise an error later when we actually try to use it
    # (if we do so) in scrapers.py.
    # In principle we could look to see if there is a matplotlib scraper
    # in our scrapers list, but this would be backward incompatible with
    # anyone using or relying on our Agg-setting behavior (e.g., for some
    # custom matplotlib SVG scraper as in our docs).
    # Eventually we can make this a config var like matplotlib_agg or something
    # if people need us not to set it to Agg.
    try:
        _import_matplotlib()
    except (ImportError, ValueError):
        pass

    # compress_images
    compress_images = gallery_conf['compress_images']
    if isinstance(compress_images, str):
        compress_images = [compress_images]
    elif not isinstance(compress_images, (tuple, list)):
        raise ConfigError('compress_images must be a tuple, list, or str, '
                          'got %s' % (type(compress_images),))
    compress_images = list(compress_images)
    allowed_values = ('images', 'thumbnails')
    pops = list()
    for ki, kind in enumerate(compress_images):
        if kind not in allowed_values:
            if kind.startswith('-'):
                pops.append(ki)
                continue
            raise ConfigError('All entries in compress_images must be one of '
                              '%s or a command-line switch starting with "-", '
                              'got %r' % (allowed_values, kind))
    compress_images_args = [compress_images.pop(p) for p in pops[::-1]]
    if len(compress_images) and not _has_optipng():
        logger.warning(
            'optipng binaries not found, PNG %s will not be optimized'
            % (' and '.join(compress_images),))
        compress_images = ()
    gallery_conf['compress_images'] = compress_images
    gallery_conf['compress_images_args'] = compress_images_args

    # deal with resetters
    resetters = gallery_conf['reset_modules']
    if not isinstance(resetters, (tuple, list)):
        resetters = [resetters]
    resetters = list(resetters)
    for ri, resetter in enumerate(resetters):
        if isinstance(resetter, str):
            if resetter not in _reset_dict:
                raise ConfigError('Unknown module resetter named %r'
                                  % (resetter,))
            resetters[ri] = _reset_dict[resetter]
        elif not callable(resetter):
            raise ConfigError('Module resetter %r was not callable'
                              % (resetter,))
    gallery_conf['reset_modules'] = tuple(resetters)

    if not isinstance(gallery_conf['reset_modules_order'], str):
        raise ConfigError('reset_modules_order must be a str, '
                          'got %r' % gallery_conf['reset_modules_order'])
    if gallery_conf['reset_modules_order'] not in ['before', 'after', 'both']:
        raise ConfigError("reset_modules_order must be in"
                          "['before', 'after', 'both'], "
                          'got %r' % gallery_conf['reset_modules_order'])

    lang = lang if lang in ('python', 'python3', 'default') else 'python'
    gallery_conf['lang'] = lang
    del resetters

    # Ensure the first cell text is a string if we have it
    first_cell = gallery_conf.get("first_notebook_cell")
    if (not isinstance(first_cell, str)) and (first_cell is not None):
        raise ConfigError("The 'first_notebook_cell' parameter must be type "
                          "str or None, found type %s" % type(first_cell))
    # Ensure the last cell text is a string if we have it
    last_cell = gallery_conf.get("last_notebook_cell")
    if (not isinstance(last_cell, str)) and (last_cell is not None):
        raise ConfigError("The 'last_notebook_cell' parameter must be type str"
                          " or None, found type %s" % type(last_cell))
    # Check pypandoc
    pypandoc = gallery_conf['pypandoc']
    if not isinstance(pypandoc, (dict, bool)):
        raise ConfigError("'pypandoc' parameter must be of type bool or dict,"
                          "got: %s." % type(pypandoc))
    gallery_conf['pypandoc'] = dict() if pypandoc is True else pypandoc
    has_pypandoc, version = _has_pypandoc()
    if isinstance(gallery_conf['pypandoc'], dict) and has_pypandoc is None:
        logger.warning("'pypandoc' not available. Using Sphinx-Gallery to "
                       "convert rst text blocks to markdown for .ipynb files.")
        gallery_conf['pypandoc'] = False
    elif isinstance(gallery_conf['pypandoc'], dict):
        logger.info("Using pandoc version: %s to convert rst text blocks to "
                    "markdown for .ipynb files" % (version,))
    else:
        logger.info("Using Sphinx-Gallery to convert rst text blocks to "
                    "markdown for .ipynb files.")
    if isinstance(pypandoc, dict):
        accepted_keys = ('extra_args', 'filters')
        for key in pypandoc:
            if key not in accepted_keys:
                raise ConfigError("'pypandoc' only accepts the following key "
                                  "values: %s, got: %s."
                                  % (accepted_keys, key))

    # Make it easy to know which builder we're in
    gallery_conf['builder_name'] = builder_name
    gallery_conf['titles'] = {}
    # Ensure 'backreferences_dir' is str, pathlib.Path or None
    backref = gallery_conf['backreferences_dir']
    if (not isinstance(backref, (str, pathlib.Path))) and \
            (backref is not None):
        raise ConfigError("The 'backreferences_dir' parameter must be of type "
                          "str, pathlib.Path or None, "
                          "found type %s" % type(backref))
    # if 'backreferences_dir' is pathlib.Path, make str for Python <=3.5
    # compatibility
    if isinstance(backref, pathlib.Path):
        gallery_conf['backreferences_dir'] = str(backref)

    # binder
    gallery_conf['binder'] = check_binder_conf(gallery_conf['binder'])

    if not isinstance(gallery_conf['css'], (list, tuple)):
        raise ConfigError('gallery_conf["css"] must be list or tuple, got %r'
                          % (gallery_conf['css'],))
    for css in gallery_conf['css']:
        if css not in _KNOWN_CSS:
            raise ConfigError('Unknown css %r, must be one of %r'
                              % (css, _KNOWN_CSS))
        if gallery_conf['app'] is not None:  # can be None in testing
            gallery_conf['app'].add_css_file(css + '.css')

    return gallery_conf


def get_subsections(srcdir, examples_dir, gallery_conf):
    """Return the list of subsections of a gallery.

    Parameters
    ----------
    srcdir : str
        absolute path to directory containing conf.py
    examples_dir : str
        path to the examples directory relative to conf.py
    gallery_conf : dict
        The gallery configuration.

    Returns
    -------
    out : list
        sorted list of gallery subsection folder names
    """
    sortkey = gallery_conf['subsection_order']
    subfolders = [subfolder for subfolder in os.listdir(examples_dir)
                  if _get_readme(os.path.join(examples_dir, subfolder),
                                 gallery_conf, raise_error=False) is not None]
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

    if bool(gallery_conf['backreferences_dir']):
        backreferences_dir = os.path.join(
            srcdir, gallery_conf['backreferences_dir'])
        if not os.path.exists(backreferences_dir):
            os.makedirs(backreferences_dir)

    return list(zip(examples_dirs, gallery_dirs))


def generate_gallery_rst(app):
    """Generate the Main examples gallery reStructuredText

    Start the Sphinx-Gallery configuration and recursively scan the examples
    directories in order to populate the examples gallery
    """
    logger.info('generating gallery...', color='white')
    gallery_conf = parse_config(app)

    seen_backrefs = set()

    costs = []
    workdirs = _prepare_sphx_glr_dirs(gallery_conf,
                                      app.builder.srcdir)

    # Check for duplicate filenames to make sure linking works as expected
    examples_dirs = [ex_dir for ex_dir, _ in workdirs]
    files = collect_gallery_files(examples_dirs, gallery_conf)
    check_duplicate_filenames(files)
    check_spaces_in_filenames(files)

    for examples_dir, gallery_dir in workdirs:

        examples_dir = os.path.join(app.builder.srcdir, examples_dir)
        gallery_dir = os.path.join(app.builder.srcdir, gallery_dir)

        # Here we don't use an os.walk, but we recurse only twice: flat is
        # better than nested.
        this_fhindex, this_costs = generate_dir_rst(
            examples_dir, gallery_dir, gallery_conf, seen_backrefs)

        costs += this_costs
        write_computation_times(gallery_conf, gallery_dir, this_costs)

        # we create an index.rst with all examples
        index_rst_new = os.path.join(gallery_dir, 'index.rst.new')
        with codecs.open(index_rst_new, 'w', encoding='utf-8') as fhindex:
            # :orphan: to suppress "not included in TOCTREE" sphinx warnings
            fhindex.write(":orphan:\n\n" + this_fhindex)

            for subsection in get_subsections(
                    app.builder.srcdir, examples_dir, gallery_conf):
                src_dir = os.path.join(examples_dir, subsection)
                target_dir = os.path.join(gallery_dir, subsection)
                this_fhindex, this_costs = \
                    generate_dir_rst(src_dir, target_dir, gallery_conf,
                                     seen_backrefs)
                fhindex.write(this_fhindex)
                costs += this_costs
                write_computation_times(gallery_conf, target_dir, this_costs)

            if gallery_conf['download_all_examples']:
                download_fhindex = generate_zipfiles(
                    gallery_dir, app.builder.srcdir)
                fhindex.write(download_fhindex)

            if (app.config.sphinx_gallery_conf['show_signature']):
                fhindex.write(SPHX_GLR_SIG)
        _replace_md5(index_rst_new, mode='t')
    _finalize_backreferences(seen_backrefs, gallery_conf)

    if gallery_conf['plot_gallery']:
        logger.info("computation time summary:", color='white')
        lines, lens = _format_for_writing(
            costs, os.path.normpath(gallery_conf['src_dir']), kind='console')
        for name, t, m in lines:
            text = ('    - %s:   ' % (name,)).ljust(lens[0] + 10)
            if t is None:
                text += '(not run)'
                logger.info(text)
            else:
                t_float = float(t.split()[0])
                if t_float >= gallery_conf['min_reported_time']:
                    text += t.rjust(lens[1]) + '   ' + m.rjust(lens[2])
                    logger.info(text)
        # Also create a junit.xml file, useful e.g. on CircleCI
        write_junit_xml(gallery_conf, app.builder.outdir, costs)


SPHX_GLR_COMP_TIMES = """
:orphan:

.. _{0}:

Computation times
=================
"""


def _sec_to_readable(t):
    """Convert a number of seconds to a more readable representation."""
    # This will only work for < 1 day execution time
    # And we reserve 2 digits for minutes because presumably
    # there aren't many > 99 minute scripts, but occasionally some
    # > 9 minute ones
    t = datetime(1, 1, 1) + timedelta(seconds=t)
    t = '{0:02d}:{1:02d}.{2:03d}'.format(
        t.hour * 60 + t.minute, t.second,
        int(round(t.microsecond / 1000.)))
    return t


def cost_name_key(cost_name):
    cost, name = cost_name
    # sort by descending computation time, descending memory, alphabetical name
    return (-cost[0], -cost[1], name)


def _format_for_writing(costs, path, kind='rst'):
    lines = list()
    for cost in sorted(costs, key=cost_name_key):
        if kind == 'rst':  # like in sg_execution_times
            name = ':ref:`sphx_glr_{0}_{1}` (``{1}``)'.format(
                path, os.path.basename(cost[1]))
            t = _sec_to_readable(cost[0][0])
        else:  # like in generate_gallery
            assert kind == 'console'
            name = os.path.relpath(cost[1], path)
            t = '%0.2f sec' % (cost[0][0],)
        m = '{0:.1f} MB'.format(cost[0][1])
        lines.append([name, t, m])
    lens = [max(x) for x in zip(*[[len(item) for item in cost]
                                  for cost in lines])]
    return lines, lens


def write_computation_times(gallery_conf, target_dir, costs):
    total_time = sum(cost[0][0] for cost in costs)
    if total_time == 0:
        return
    target_dir_clean = os.path.relpath(
        target_dir, gallery_conf['src_dir']).replace(os.path.sep, '_')
    new_ref = 'sphx_glr_%s_sg_execution_times' % target_dir_clean
    with codecs.open(os.path.join(target_dir, 'sg_execution_times.rst'), 'w',
                     encoding='utf-8') as fid:
        fid.write(SPHX_GLR_COMP_TIMES.format(new_ref))
        fid.write('**{0}** total execution time for **{1}** files:\n\n'
                  .format(_sec_to_readable(total_time), target_dir_clean))
        lines, lens = _format_for_writing(costs, target_dir_clean)
        del costs
        hline = ''.join(('+' + '-' * (length + 2)) for length in lens) + '+\n'
        fid.write(hline)
        format_str = ''.join('| {%s} ' % (ii,)
                             for ii in range(len(lines[0]))) + '|\n'
        for line in lines:
            line = [ll.ljust(len_) for ll, len_ in zip(line, lens)]
            text = format_str.format(*line)
            assert len(text) == len(hline)
            fid.write(text)
            fid.write(hline)


def write_junit_xml(gallery_conf, target_dir, costs):
    if not gallery_conf['junit'] or not gallery_conf['plot_gallery']:
        return
    failing_as_expected, failing_unexpectedly, passing_unexpectedly = \
        _parse_failures(gallery_conf)
    n_tests = 0
    n_failures = 0
    n_skips = 0
    elapsed = 0.
    src_dir = gallery_conf['src_dir']
    output = ''
    for cost in costs:
        (t, _), fname = cost
        if not any(fname in x for x in (gallery_conf['passing_examples'],
                                        failing_unexpectedly,
                                        failing_as_expected,
                                        passing_unexpectedly)):
            continue  # not subselected by our regex
        title = gallery_conf['titles'][fname]
        output += (
            u'<testcase classname={0!s} file={1!s} line="1" '
            u'name={2!s} time="{3!r}">'
            .format(quoteattr(os.path.splitext(os.path.basename(fname))[0]),
                    quoteattr(os.path.relpath(fname, src_dir)),
                    quoteattr(title), t))
        if fname in failing_as_expected:
            output += u'<skipped message="expected example failure"></skipped>'
            n_skips += 1
        elif fname in failing_unexpectedly or fname in passing_unexpectedly:
            if fname in failing_unexpectedly:
                traceback = gallery_conf['failing_examples'][fname]
            else:  # fname in passing_unexpectedly
                traceback = 'Passed even though it was marked to fail'
            n_failures += 1
            output += (u'<failure message={0!s}>{1!s}</failure>'
                       .format(quoteattr(traceback.splitlines()[-1].strip()),
                               escape(traceback)))
        output += u'</testcase>'
        n_tests += 1
        elapsed += t
    output += u'</testsuite>'
    output = (u'<?xml version="1.0" encoding="utf-8"?>'
              u'<testsuite errors="0" failures="{0}" name="sphinx-gallery" '
              u'skipped="{1}" tests="{2}" time="{3}">'
              .format(n_failures, n_skips, n_tests, elapsed)) + output
    # Actually write it
    fname = os.path.normpath(os.path.join(target_dir, gallery_conf['junit']))
    junit_dir = os.path.dirname(fname)
    if not os.path.isdir(junit_dir):
        os.makedirs(junit_dir)
    with codecs.open(fname, 'w', encoding='utf-8') as fid:
        fid.write(output)


def touch_empty_backreferences(app, what, name, obj, options, lines):
    """Generate empty back-reference example files.

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


def _expected_failing_examples(gallery_conf):
    return set(
        os.path.normpath(os.path.join(gallery_conf['src_dir'], path))
        for path in gallery_conf['expected_failing_examples'])


def _parse_failures(gallery_conf):
    """Split the failures."""
    failing_examples = set(gallery_conf['failing_examples'].keys())
    expected_failing_examples = _expected_failing_examples(gallery_conf)
    failing_as_expected = failing_examples.intersection(
        expected_failing_examples)
    failing_unexpectedly = failing_examples.difference(
        expected_failing_examples)
    passing_unexpectedly = expected_failing_examples.difference(
        failing_examples)
    # filter from examples actually run
    passing_unexpectedly = [
        src_file for src_file in passing_unexpectedly
        if re.search(gallery_conf.get('filename_pattern'), src_file)]
    return failing_as_expected, failing_unexpectedly, passing_unexpectedly


def summarize_failing_examples(app, exception):
    """Collects the list of falling examples and prints them with a traceback.

    Raises ValueError if there where failing examples.
    """
    if exception is not None:
        return

    # Under no-plot Examples are not run so nothing to summarize
    if not app.config.sphinx_gallery_conf['plot_gallery']:
        logger.info('Sphinx-Gallery gallery_conf["plot_gallery"] was '
                    'False, so no examples were executed.', color='brown')
        return

    gallery_conf = app.config.sphinx_gallery_conf
    failing_as_expected, failing_unexpectedly, passing_unexpectedly = \
        _parse_failures(gallery_conf)

    if failing_as_expected:
        logger.info("Examples failing as expected:", color='brown')
        for fail_example in failing_as_expected:
            logger.info('%s failed leaving traceback:', fail_example,
                        color='brown')
            logger.info(gallery_conf['failing_examples'][fail_example],
                        color='brown')

    fail_msgs = []
    if failing_unexpectedly:
        fail_msgs.append(red("Unexpected failing examples:"))
        for fail_example in failing_unexpectedly:
            fail_msgs.append(fail_example + ' failed leaving traceback:\n' +
                             gallery_conf['failing_examples'][fail_example] +
                             '\n')

    if passing_unexpectedly:
        fail_msgs.append(red("Examples expected to fail, but not failing:\n") +
                         "Please remove these examples from\n" +
                         "sphinx_gallery_conf['expected_failing_examples']\n" +
                         "in your conf.py file"
                         "\n".join(passing_unexpectedly))

    # standard message
    n_good = len(gallery_conf['passing_examples'])
    n_tot = len(gallery_conf['failing_examples']) + n_good
    n_stale = len(gallery_conf['stale_examples'])
    logger.info('\nSphinx-Gallery successfully executed %d out of %d '
                'file%s subselected by:\n\n'
                '    gallery_conf["filename_pattern"] = %r\n'
                '    gallery_conf["ignore_pattern"]   = %r\n'
                '\nafter excluding %d file%s that had previously been run '
                '(based on MD5).\n'
                % (n_good, n_tot, 's' if n_tot != 1 else '',
                   gallery_conf['filename_pattern'],
                   gallery_conf['ignore_pattern'],
                   n_stale, 's' if n_stale != 1 else '',
                   ),
                color='brown')

    if fail_msgs:
        fail_message = ("Here is a summary of the problems encountered "
                        "when running the examples\n\n" +
                        "\n".join(fail_msgs) + "\n" + "-" * 79)
        if gallery_conf['only_warn_on_example_error']:
            logger.warning(fail_message)
        else:
            raise ExtensionError(fail_message)


def collect_gallery_files(examples_dirs, gallery_conf):
    """Collect python files from the gallery example directories."""
    files = []
    for example_dir in examples_dirs:
        for root, dirnames, filenames in os.walk(example_dir):
            for filename in filenames:
                if filename.endswith('.py'):
                    if re.search(gallery_conf['ignore_pattern'],
                                 filename) is None:
                        files.append(os.path.join(root, filename))
    return files


def check_duplicate_filenames(files):
    """Check for duplicate filenames across gallery directories."""
    # Check whether we'll have duplicates
    used_names = set()
    dup_names = list()

    for this_file in files:
        this_fname = os.path.basename(this_file)
        if this_fname in used_names:
            dup_names.append(this_file)
        else:
            used_names.add(this_fname)

    if len(dup_names) > 0:
        logger.warning(
            'Duplicate example file name(s) found. Having duplicate file '
            'names will break some links. '
            'List of files: {}'.format(sorted(dup_names),))


def check_spaces_in_filenames(files):
    """Check for spaces in filenames across example directories."""
    regex = re.compile(r'[\s]')
    files_with_space = list(filter(regex.search, files))
    if files_with_space:
        logger.warning(
            'Example file name(s) with space(s) found. Having space(s) in '
            'file names will break some links. '
            'List of files: {}'.format(sorted(files_with_space),))


def get_default_config_value(key):
    def default_getter(conf):
        return conf['sphinx_gallery_conf'].get(key, DEFAULT_GALLERY_CONF[key])
    return default_getter


def setup(app):
    """Setup Sphinx-Gallery sphinx extension"""
    sphinx_compatibility._app = app

    app.add_config_value('sphinx_gallery_conf', DEFAULT_GALLERY_CONF, 'html')
    for key in ['plot_gallery', 'abort_on_example_error']:
        app.add_config_value(key, get_default_config_value(key), 'html')

    if 'sphinx.ext.autodoc' in app.extensions:
        app.connect('autodoc-process-docstring', touch_empty_backreferences)

    # Add the custom directive
    app.add_directive('minigallery', MiniGallery)
    app.add_directive("image-sg", ImageSg)

    imagesg_addnode(app)

    app.connect('builder-inited', generate_gallery_rst)
    app.connect('build-finished', copy_binder_files)
    app.connect('build-finished', summarize_failing_examples)
    app.connect('build-finished', embed_code_links)
    metadata = {'parallel_read_safe': True,
                'parallel_write_safe': True,
                'version': _sg_version}
    return metadata


def setup_module():
    # HACK: Stop nosetests running setup() above
    pass
