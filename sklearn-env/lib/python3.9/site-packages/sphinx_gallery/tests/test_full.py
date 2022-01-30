# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
Test the SG pipeline used with Sphinx
"""
from __future__ import division, absolute_import, print_function

import codecs
from io import StringIO
import os
import os.path as op
import re
import shutil
import sys
import time

from sphinx.application import Sphinx
from sphinx.errors import ExtensionError
from sphinx.util.docutils import docutils_namespace
from sphinx_gallery.utils import (_get_image, scale_image, _has_optipng,
                                  _has_pypandoc)

import pytest

N_TOT = 13

N_FAILING = 2
N_GOOD = N_TOT - N_FAILING
N_RST = 15 + N_TOT
N_RST = '(%s|%s)' % (N_RST, N_RST - 1)  # AppVeyor weirdness


@pytest.fixture(scope='module')
def sphinx_app(tmpdir_factory, req_mpl, req_pil):
    # Skip if numpy not installed
    pytest.importorskip("numpy")

    temp_dir = (tmpdir_factory.getbasetemp() / 'root').strpath
    src_dir = op.join(op.dirname(__file__), 'tinybuild')

    def ignore(src, names):
        return ('_build', 'gen_modules', 'auto_examples')

    shutil.copytree(src_dir, temp_dir, ignore=ignore)
    # For testing iteration, you can get similar behavior just doing `make`
    # inside the tinybuild directory
    src_dir = temp_dir
    conf_dir = temp_dir
    out_dir = op.join(temp_dir, '_build', 'html')
    toctrees_dir = op.join(temp_dir, '_build', 'toctrees')
    # Avoid warnings about re-registration, see:
    # https://github.com/sphinx-doc/sphinx/issues/5038
    with docutils_namespace():
        app = Sphinx(src_dir, conf_dir, out_dir, toctrees_dir,
                     buildername='html', status=StringIO(),
                     warning=StringIO())
        # need to build within the context manager
        # for automodule and backrefs to work
        app.build(False, [])
    return app


def test_timings(sphinx_app):
    """Test that a timings page is created."""
    out_dir = sphinx_app.outdir
    src_dir = sphinx_app.srcdir
    # local folder
    timings_rst = op.join(src_dir, 'auto_examples',
                          'sg_execution_times.rst')
    assert op.isfile(timings_rst)
    with codecs.open(timings_rst, 'r', 'utf-8') as fid:
        content = fid.read()
    assert ':ref:`sphx_glr_auto_examples_plot_numpy_matplotlib.py`' in content
    parenthetical = '(``%s``)' % ('plot_numpy_matplotlib.py',)
    assert parenthetical in content
    # HTML output
    timings_html = op.join(out_dir, 'auto_examples',
                           'sg_execution_times.html')
    assert op.isfile(timings_html)
    with codecs.open(timings_html, 'r', 'utf-8') as fid:
        content = fid.read()
    assert 'href="plot_numpy_matplotlib.html' in content
    # printed
    status = sphinx_app._status.getvalue()
    fname = op.join('examples', 'plot_numpy_matplotlib.py')
    assert ('- %s: ' % fname) in status


def test_optipng(sphinx_app):
    """Test that optipng is detected."""
    status = sphinx_app._status.getvalue()
    w = sphinx_app._warning.getvalue()
    substr = 'will not be optimized'
    if _has_optipng():
        assert substr not in w
    else:
        assert substr in w
    assert 'optipng version' not in status.lower()  # catch the --version


def test_junit(sphinx_app, tmpdir):
    out_dir = sphinx_app.outdir
    junit_file = op.join(out_dir, 'sphinx-gallery', 'junit-results.xml')
    assert op.isfile(junit_file)
    with codecs.open(junit_file, 'r', 'utf-8') as fid:
        contents = fid.read()
    assert contents.startswith('<?xml')
    assert 'errors="0" failures="0"' in contents
    assert 'tests="%d"' % (N_TOT,) in contents
    assert 'local_module' not in contents  # it's not actually run as an ex
    assert 'expected example failure' in contents
    assert '<failure message' not in contents
    src_dir = sphinx_app.srcdir
    new_src_dir = op.join(str(tmpdir), 'src')
    shutil.copytree(src_dir, new_src_dir)
    del src_dir
    new_out_dir = op.join(new_src_dir, '_build', 'html')
    new_toctree_dir = op.join(new_src_dir, '_build', 'toctrees')
    passing_fname = op.join(new_src_dir, 'examples',
                            'plot_numpy_matplotlib.py')
    failing_fname = op.join(new_src_dir, 'examples', 'future',
                            'plot_future_imports_broken.py')
    shutil.move(passing_fname, passing_fname + '.temp')
    shutil.move(failing_fname, passing_fname)
    shutil.move(passing_fname + '.temp', failing_fname)
    with docutils_namespace():
        app = Sphinx(new_src_dir, new_src_dir, new_out_dir,
                     new_toctree_dir,
                     buildername='html', status=StringIO())
        # need to build within the context manager
        # for automodule and backrefs to work
        with pytest.raises(ExtensionError, match='Here is a summary of the '):
            app.build(False, [])
    junit_file = op.join(new_out_dir, 'sphinx-gallery', 'junit-results.xml')
    assert op.isfile(junit_file)
    with codecs.open(junit_file, 'r', 'utf-8') as fid:
        contents = fid.read()
    assert 'errors="0" failures="2"' in contents
    # this time we only ran the stale files
    assert 'tests="%s"' % (N_FAILING + 1,) in contents
    assert '<failure message="RuntimeError: Forcing' in contents
    assert 'Passed even though it was marked to fail' in contents


def test_run_sphinx(sphinx_app):
    """Test basic outputs."""
    out_dir = sphinx_app.outdir
    out_files = os.listdir(out_dir)
    assert 'index.html' in out_files
    assert 'auto_examples' in out_files
    generated_examples_dir = op.join(out_dir, 'auto_examples')
    assert op.isdir(generated_examples_dir)
    status = sphinx_app._status.getvalue()
    assert 'executed %d out of %d' % (N_GOOD, N_TOT) in status
    assert 'after excluding 0' in status
    # intentionally have a bad URL in references
    warning = sphinx_app._warning.getvalue()
    want = '.*fetching .*wrong_url.*404.*'
    assert re.match(want, warning, re.DOTALL) is not None, warning


def test_thumbnail_path(sphinx_app, tmpdir):
    """Test sphinx_gallery_thumbnail_path."""
    import numpy as np
    # Make sure our thumbnail matches what it should be
    fname_orig = op.join(
        sphinx_app.srcdir, '_static_nonstandard', 'demo.png')
    fname_thumb = op.join(
        sphinx_app.outdir, '_images',
        'sphx_glr_plot_second_future_imports_thumb.png')
    fname_new = str(tmpdir.join('new.png'))
    scale_image(fname_orig, fname_new,
                *sphinx_app.config.sphinx_gallery_conf["thumbnail_size"])
    Image = _get_image()
    orig = np.asarray(Image.open(fname_thumb))
    new = np.asarray(Image.open(fname_new))
    assert new.shape[:2] == orig.shape[:2]
    assert new.shape[2] in (3, 4)  # optipng can strip the alpha channel
    corr = np.corrcoef(new[..., :3].ravel(), orig[..., :3].ravel())[0, 1]
    assert corr > 0.99


def test_negative_thumbnail_config(sphinx_app, tmpdir):
    """Test 'sphinx_gallery_thumbnail_number' config works correctly for
    negative numbers."""
    import numpy as np
    # Make sure our thumbnail is the 2nd (last) image
    fname_orig = op.join(
        sphinx_app.outdir, '_images',
        'sphx_glr_plot_matplotlib_alt_002.png')
    fname_thumb = op.join(
        sphinx_app.outdir, '_images',
        'sphx_glr_plot_matplotlib_alt_thumb.png')
    fname_new = str(tmpdir.join('new.png'))
    scale_image(fname_orig, fname_new,
                *sphinx_app.config.sphinx_gallery_conf["thumbnail_size"])
    Image = _get_image()
    orig = np.asarray(Image.open(fname_thumb))
    new = np.asarray(Image.open(fname_new))
    assert new.shape[:2] == orig.shape[:2]
    assert new.shape[2] in (3, 4)  # optipng can strip the alpha channel
    corr = np.corrcoef(new[..., :3].ravel(), orig[..., :3].ravel())[0, 1]
    assert corr > 0.99


def test_command_line_args_img(sphinx_app):
    generated_examples_dir = op.join(sphinx_app.outdir, 'auto_examples')
    thumb_fname = '../_images/sphx_glr_plot_command_line_args_thumb.png'
    file_fname = op.join(generated_examples_dir, thumb_fname)
    assert op.isfile(file_fname), file_fname


def test_image_formats(sphinx_app):
    """Test Image format support."""
    generated_examples_dir = op.join(sphinx_app.outdir, 'auto_examples')
    generated_examples_index = op.join(generated_examples_dir, 'index.html')
    with codecs.open(generated_examples_index, 'r', 'utf-8') as fid:
        html = fid.read()
    thumb_fnames = ['../_images/sphx_glr_plot_svg_thumb.svg',
                    '../_images/sphx_glr_plot_numpy_matplotlib_thumb.png',
                    '../_images/sphx_glr_plot_animation_thumb.gif',
                    ]
    for thumb_fname in thumb_fnames:
        file_fname = op.join(generated_examples_dir, thumb_fname)
        assert op.isfile(file_fname), file_fname
        want_html = 'src="%s"' % (thumb_fname,)
        assert want_html in html
    # the original GIF does not get copied because it's not used in the
    # RST/HTML, so can't add it to this check
    for ex, ext, nums, extra in (
            ('plot_svg', 'svg', [1], None),
            ('plot_numpy_matplotlib', 'png', [1], None),
            ('plot_animation', 'png', [1, 3], 'function Animation')):
        html_fname = op.join(generated_examples_dir, '%s.html' % ex)
        with codecs.open(html_fname, 'r', 'utf-8') as fid:
            html = fid.read()
        for num in nums:
            img_fname0 = '../_images/sphx_glr_%s_%03d.%s' % (ex, num, ext)
            file_fname = op.join(generated_examples_dir, img_fname0)
            assert op.isfile(file_fname), file_fname
            want_html = 'src="%s"' % (img_fname0,)
            assert want_html in html
            img_fname2 = ('../_images/sphx_glr_%s_%03d_2_0x.%s' %
                          (ex, num, ext))
            file_fname2 = op.join(generated_examples_dir, img_fname2)
            want_html = 'srcset="%s, %s 2.0x"' % (img_fname0, img_fname2)
            if ext in ('png', 'jpg', 'svg'):  # check 2.0x (tests directive)
                assert op.isfile(file_fname2), file_fname2
                assert want_html in html

        if extra is not None:
            assert extra in html


def test_repr_html_classes(sphinx_app):
    """Test appropriate _repr_html_ classes."""
    example_file = op.join(
        sphinx_app.outdir, 'auto_examples', 'plot_repr.html')
    with codecs.open(example_file, 'r', 'utf-8') as fid:
        lines = fid.read()
    assert 'div class="output_subarea output_html rendered_html output_result"' in lines  # noqa: E501
    assert 'gallery-rendered-html.css' in lines


def test_embed_links_and_styles(sphinx_app):
    """Test that links and styles are embedded properly in doc."""
    out_dir = sphinx_app.outdir
    src_dir = sphinx_app.srcdir
    examples_dir = op.join(out_dir, 'auto_examples')
    assert op.isdir(examples_dir)
    example_files = os.listdir(examples_dir)
    assert 'plot_numpy_matplotlib.html' in example_files
    example_file = op.join(examples_dir, 'plot_numpy_matplotlib.html')
    with codecs.open(example_file, 'r', 'utf-8') as fid:
        lines = fid.read()
    # ensure we've linked properly
    assert '#module-matplotlib.colors' in lines
    assert 'matplotlib.colors.is_color_like' in lines
    assert 'class="sphx-glr-backref-module-matplotlib-colors sphx-glr-backref-type-py-function">' in lines  # noqa: E501
    assert '#module-numpy' in lines
    assert 'numpy.arange.html' in lines
    assert 'class="sphx-glr-backref-module-numpy sphx-glr-backref-type-py-function">' in lines  # noqa: E501
    assert '#module-matplotlib.pyplot' in lines
    assert 'pyplot.html' in lines
    assert '.html#matplotlib.figure.Figure.tight_layout' in lines
    assert 'matplotlib.axes.Axes.plot.html#matplotlib.axes.Axes.plot' in lines
    assert 'matplotlib_configuration_api.html#matplotlib.RcParams' in lines
    assert 'stdtypes.html#list' in lines
    assert 'warnings.html#warnings.warn' in lines
    assert 'itertools.html#itertools.compress' in lines
    assert 'numpy.ndarray.html' in lines
    # see issue 617
    id_names = re.search(
        r'sphinx_gallery.backreferences.html#sphinx[_,-]gallery[.,-]backreferences[.,-]identify[_,-]names',  # noqa: E501
        lines
    )
    assert id_names is not None
    # instances have an extra CSS class
    assert 'class="sphx-glr-backref-module-matplotlib-figure sphx-glr-backref-type-py-class sphx-glr-backref-instance"><span class="n">x</span></a>' in lines  # noqa: E501
    assert 'class="sphx-glr-backref-module-matplotlib-figure sphx-glr-backref-type-py-class"><span class="n">Figure</span></a>' in lines  # noqa: E501
    # gh-587: no classes that are only marked as module without type
    assert re.search(r'"sphx-glr-backref-module-\S*"', lines) is None
    assert 'class="sphx-glr-backref-module-sphinx_gallery-backreferences sphx-glr-backref-type-py-function"><span class="n">sphinx_gallery</span><span class="o">.</span><span class="n">backreferences</span><span class="o">.</span><span class="n">identify_names</span></a>' in lines  # noqa: E501
    # gh-587: np.random.RandomState links properly
    # NumPy has had this linked as numpy.random.RandomState and
    # numpy.random.mtrand.RandomState so we need regex...
    assert re.search(r'\.html#numpy\.random\.(mtrand\.?)?RandomState" title="numpy\.random\.(mtrand\.?)?RandomState" class="sphx-glr-backref-module-numpy-random(-mtrand?)? sphx-glr-backref-type-py-class"><span class="n">np</span>', lines) is not None  # noqa: E501
    assert re.search(r'\.html#numpy\.random\.(mtrand\.?)?RandomState" title="numpy\.random\.(mtrand\.?)?RandomState" class="sphx-glr-backref-module-numpy-random(-mtrand?)? sphx-glr-backref-type-py-class sphx-glr-backref-instance"><span class="n">rng</span></a>', lines) is not None  # noqa: E501
    # gh-587: methods of classes in the module currently being documented
    # issue 617 (regex '-'s)
    # instance
    dummy_class_inst = re.search(
        r'sphinx_gallery.backreferences.html#sphinx[_,-]gallery[.,-]backreferences[.,-][D,d]ummy[C,c]lass" title="sphinx_gallery.backreferences.DummyClass" class="sphx-glr-backref-module-sphinx_gallery-backreferences sphx-glr-backref-type-py-class sphx-glr-backref-instance"><span class="n">dc</span>',  # noqa: E501
        lines
    )
    assert dummy_class_inst is not None
    # class
    dummy_class_class = re.search(
        r'sphinx_gallery.backreferences.html#sphinx[_,-]gallery[.,-]backreferences[.,-][D,d]ummy[C,c]lass" title="sphinx_gallery.backreferences.DummyClass" class="sphx-glr-backref-module-sphinx_gallery-backreferences sphx-glr-backref-type-py-class"><span class="n">sphinx_gallery</span><span class="o">.</span><span class="n">backreferences</span><span class="o">.</span><span class="n">DummyClass</span>',  # noqa: E501
        lines
    )
    assert dummy_class_class is not None
    # method
    dummy_class_meth = re.search(
        r'sphinx_gallery.backreferences.html#sphinx[_,-]gallery[.,-]backreferences[.,-][D,d]ummy[C,c]lass[.,-]run" title="sphinx_gallery.backreferences.DummyClass.run" class="sphx-glr-backref-module-sphinx_gallery-backreferences sphx-glr-backref-type-py-method"><span class="n">dc</span><span class="o">.</span><span class="n">run</span>',  # noqa: E501
        lines
    )
    assert dummy_class_meth is not None
    # property (Sphinx 2+ calls it a method rather than attribute, so regex)
    dummy_class_prop = re.compile(r'sphinx_gallery.backreferences.html#sphinx[_,-]gallery[.,-]backreferences[.,-][D,d]ummy[C,c]lass[.,-]prop" title="sphinx_gallery.backreferences.DummyClass.prop" class="sphx-glr-backref-module-sphinx_gallery-backreferences sphx-glr-backref-type-py-(attribute|method|property)"><span class="n">dc</span><span class="o">.</span><span class="n">prop</span>')  # noqa: E501
    assert dummy_class_prop.search(lines) is not None

    try:
        import memory_profiler  # noqa, analysis:ignore
    except ImportError:
        assert "memory usage" not in lines
    else:
        assert "memory usage" in lines

    # CSS styles
    assert 'class="sphx-glr-signature"' in lines
    assert 'class="sphx-glr-timing"' in lines
    for kind in ('python', 'jupyter'):
        assert 'class="sphx-glr-download sphx-glr-download-%s docutils container"' % kind in lines  # noqa:E501

    # highlight language
    fname = op.join(src_dir, 'auto_examples', 'plot_numpy_matplotlib.rst')
    assert op.isfile(fname)
    with codecs.open(fname, 'r', 'utf-8') as fid:
        rst = fid.read()
    assert '.. code-block:: python3\n' in rst

    # warnings
    want_warn = (r'.*plot_numpy_matplotlib\.py:[0-9][0-9]: RuntimeWarning: '
                 r'This warning should show up in the output.*')
    assert re.match(want_warn, lines, re.DOTALL) is not None
    sys.stdout.write(lines)

    example_file = op.join(examples_dir, 'plot_pickle.html')
    with codecs.open(example_file, 'r', 'utf-8') as fid:
        lines = fid.read()
    assert 'joblib.Parallel.html' in lines


def test_backreferences(sphinx_app):
    """Test backreferences."""
    out_dir = sphinx_app.outdir
    mod_file = op.join(out_dir, 'gen_modules', 'sphinx_gallery.sorting.html')
    with codecs.open(mod_file, 'r', 'utf-8') as fid:
        lines = fid.read()
    assert 'ExplicitOrder' in lines  # in API doc
    assert 'plot_second_future_imports.html' in lines  # backref via code use
    assert 'FileNameSortKey' in lines  # in API doc
    assert 'plot_numpy_matplotlib.html' in lines  # backref via :class: in str
    mod_file = op.join(out_dir, 'gen_modules',
                       'sphinx_gallery.backreferences.html')
    with codecs.open(mod_file, 'r', 'utf-8') as fid:
        lines = fid.read()
    assert 'NameFinder' in lines  # in API doc
    assert 'plot_future_imports.html' in lines  # backref via doc block


@pytest.mark.parametrize('rst_file, example_used_in', [
    ('sphinx_gallery.backreferences.identify_names.examples',
     'plot_numpy_matplotlib'),
    ('sphinx_gallery.sorting.ExplicitOrder.examples',
     'plot_second_future_imports'),
])
def test_backreferences_examples(sphinx_app, rst_file, example_used_in):
    """Test linking to mini-galleries using backreferences_dir."""
    backref_dir = sphinx_app.srcdir
    examples_rst = op.join(backref_dir, 'gen_modules', 'backreferences',
                           rst_file)
    with codecs.open(examples_rst, 'r', 'utf-8') as fid:
        lines = fid.read()
    assert example_used_in in lines


def test_logging_std_nested(sphinx_app):
    """Test that nested stdout/stderr uses within a given script work."""
    log_rst = op.join(
        sphinx_app.srcdir, 'auto_examples', 'plot_log.rst')
    with codecs.open(log_rst, 'r', 'utf-8') as fid:
        lines = fid.read()
    assert '.. code-block:: none\n\n    is in the same cell' in lines
    assert '.. code-block:: none\n\n    is not in the same cell' in lines


def _assert_mtimes(list_orig, list_new, different=(), ignore=()):
    """Assert that the correct set of files were changed based on mtime."""
    import numpy as np
    from numpy.testing import assert_allclose

    assert ([op.basename(x) for x in list_orig] ==
            [op.basename(x) for x in list_new])
    for orig, new in zip(list_orig, list_new):
        check_name = op.splitext(op.basename(orig))[0]
        if check_name.endswith('_codeobj'):
            check_name = check_name[:-8]
        if check_name in different:
            assert np.abs(op.getmtime(orig) - op.getmtime(new)) > 0.1
        elif check_name not in ignore:
            assert_allclose(op.getmtime(orig), op.getmtime(new),
                            atol=1e-3, rtol=1e-20,
                            err_msg=op.basename(orig))


def test_rebuild(tmpdir_factory, sphinx_app):
    # Make sure that examples that haven't been changed aren't run twice.

    #
    # First run completes in the fixture.
    #
    status = sphinx_app._status.getvalue()
    lines = [line for line in status.split('\n') if 'removed' in line]
    want = '.*%s added, 0 changed, 0 removed.*' % (N_RST,)
    assert re.match(want, status, re.MULTILINE | re.DOTALL) is not None, lines
    want = '.*targets for 3 source files that are out of date$.*'
    lines = [line for line in status.split('\n') if 'out of date' in line]
    assert re.match(want, status, re.MULTILINE | re.DOTALL) is not None, lines
    lines = [line for line in status.split('\n') if 'on MD5' in line]
    want = ('.*executed %d out of %d.*after excluding 0 files.*based on MD5.*'
            % (N_GOOD, N_TOT))
    assert re.match(want, status, re.MULTILINE | re.DOTALL) is not None, lines
    old_src_dir = (tmpdir_factory.getbasetemp() / 'root_old').strpath
    shutil.copytree(sphinx_app.srcdir, old_src_dir)
    generated_modules_0 = sorted(
        op.join(old_src_dir, 'gen_modules', f)
        for f in os.listdir(op.join(old_src_dir, 'gen_modules'))
        if op.isfile(op.join(old_src_dir, 'gen_modules', f)))
    generated_backrefs_0 = sorted(
        op.join(old_src_dir, 'gen_modules', 'backreferences', f)
        for f in os.listdir(op.join(old_src_dir, 'gen_modules',
                                    'backreferences')))
    generated_rst_0 = sorted(
        op.join(old_src_dir, 'auto_examples', f)
        for f in os.listdir(op.join(old_src_dir, 'auto_examples'))
        if f.endswith('.rst'))
    generated_pickle_0 = sorted(
        op.join(old_src_dir, 'auto_examples', f)
        for f in os.listdir(op.join(old_src_dir, 'auto_examples'))
        if f.endswith('.pickle'))
    copied_py_0 = sorted(
        op.join(old_src_dir, 'auto_examples', f)
        for f in os.listdir(op.join(old_src_dir, 'auto_examples'))
        if f.endswith('.py'))
    copied_ipy_0 = sorted(
        op.join(old_src_dir, 'auto_examples', f)
        for f in os.listdir(op.join(old_src_dir, 'auto_examples'))
        if f.endswith('.ipynb'))
    assert len(generated_modules_0) > 0
    assert len(generated_backrefs_0) > 0
    assert len(generated_rst_0) > 0
    assert len(generated_pickle_0) > 0
    assert len(copied_py_0) > 0
    assert len(copied_ipy_0) > 0
    assert len(sphinx_app.config.sphinx_gallery_conf['stale_examples']) == 0
    assert op.isfile(op.join(sphinx_app.outdir, '_images',
                             'sphx_glr_plot_numpy_matplotlib_001.png'))

    #
    # run a second time, no files should be updated
    #

    src_dir = sphinx_app.srcdir
    del sphinx_app  # don't accidentally use it below
    conf_dir = src_dir
    out_dir = op.join(src_dir, '_build', 'html')
    toctrees_dir = op.join(src_dir, '_build', 'toctrees')
    time.sleep(0.1)
    with docutils_namespace():
        new_app = Sphinx(src_dir, conf_dir, out_dir, toctrees_dir,
                         buildername='html', status=StringIO())
        new_app.build(False, [])
    status = new_app._status.getvalue()
    lines = [line for line in status.split('\n') if '0 removed' in line]
    # XXX on AppVeyor this can be 13
    if sys.platform.startswith('win'):
        assert re.match('.*[0|1] added, [1-9][0-3]? changed, 0 removed$.*',
                        status, re.MULTILINE | re.DOTALL) is not None, lines
    else:
        assert re.match('.*[0|1] added, [1-9] changed, 0 removed$.*',
                        status, re.MULTILINE | re.DOTALL) is not None, lines
    want = ('.*executed 0 out of %s.*after excluding %s files.*based on MD5.*'
            % (N_FAILING, N_GOOD))
    assert re.match(want, status, re.MULTILINE | re.DOTALL) is not None
    n_stale = len(new_app.config.sphinx_gallery_conf['stale_examples'])
    assert n_stale == N_GOOD
    assert op.isfile(op.join(new_app.outdir, '_images',
                             'sphx_glr_plot_numpy_matplotlib_001.png'))

    generated_modules_1 = sorted(
        op.join(new_app.srcdir, 'gen_modules', f)
        for f in os.listdir(op.join(new_app.srcdir, 'gen_modules'))
        if op.isfile(op.join(new_app.srcdir, 'gen_modules', f)))
    generated_backrefs_1 = sorted(
        op.join(new_app.srcdir, 'gen_modules', 'backreferences', f)
        for f in os.listdir(op.join(new_app.srcdir, 'gen_modules',
                                    'backreferences')))
    generated_rst_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.rst'))
    generated_pickle_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.pickle'))
    copied_py_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.py'))
    copied_ipy_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.ipynb'))

    # mtimes for modules
    _assert_mtimes(generated_modules_0, generated_modules_1)

    # mtimes for backrefs (gh-394)
    _assert_mtimes(generated_backrefs_0, generated_backrefs_1)

    # generated RST files
    ignore = (
        # these two should almost always be different, but in case we
        # get extremely unlucky and have identical run times
        # on the one script that gets re-run (because it's a fail)...
        'sg_execution_times',
        'plot_future_imports_broken',
        'plot_scraper_broken'
    )
    _assert_mtimes(generated_rst_0, generated_rst_1, ignore=ignore)

    # mtimes for pickles
    _assert_mtimes(generated_pickle_0, generated_pickle_1)

    # mtimes for .py files (gh-395)
    _assert_mtimes(copied_py_0, copied_py_1)

    # mtimes for .ipynb files
    _assert_mtimes(copied_ipy_0, copied_ipy_1)

    #
    # run a third and a fourth time, changing one file or running one stale
    #

    for how in ('run_stale', 'modify'):
        # modify must be last as this rerun setting tries to run the
        # broken example (subsequent tests depend on it)
        _rerun(how, src_dir, conf_dir, out_dir, toctrees_dir,
               generated_modules_0, generated_backrefs_0, generated_rst_0,
               generated_pickle_0, copied_py_0, copied_ipy_0)


def _rerun(how, src_dir, conf_dir, out_dir, toctrees_dir,
           generated_modules_0, generated_backrefs_0, generated_rst_0,
           generated_pickle_0, copied_py_0, copied_ipy_0):
    """Rerun the sphinx build and check that the right files were changed."""
    time.sleep(0.1)
    confoverrides = dict()
    if how == 'modify':
        fname = op.join(src_dir, 'examples', 'plot_numpy_matplotlib.py')
        with codecs.open(fname, 'r', 'utf-8') as fid:
            lines = fid.readlines()
        with codecs.open(fname, 'w', 'utf-8') as fid:
            for line in lines:
                if 'FYI this' in line:
                    line = 'A ' + line
                fid.write(line)
        out_of, excluding = N_FAILING + 1, N_GOOD - 1
        n_stale = N_GOOD - 1
    else:
        assert how == 'run_stale'
        confoverrides['sphinx_gallery_conf.run_stale_examples'] = 'True'
        confoverrides['sphinx_gallery_conf.filename_pattern'] = 'plot_numpy_ma'
        out_of, excluding = 1, 0
        n_stale = 0
    with docutils_namespace():
        new_app = Sphinx(src_dir, conf_dir, out_dir, toctrees_dir,
                         buildername='html', status=StringIO(),
                         confoverrides=confoverrides)
        new_app.build(False, [])
    status = new_app._status.getvalue()
    lines = [line for line in status.split('\n') if 'source files tha' in line]
    lines = '\n'.join([how] + lines)
    flags = re.MULTILINE | re.DOTALL
    # for some reason, setting "confoverrides" above causes Sphinx to show
    # all targets out of date, even though they haven't been modified...
    want = '.*targets for %s source files that are out of date$.*' % N_RST
    assert re.match(want, status, flags) is not None, lines
    # ... but then later detects that only some have actually changed:
    # Linux: 8 changed when how='run_stale', 9 when how='modify'.
    # Windows: always 9 for some reason
    lines = [line for line in status.split('\n') if 'changed,' in line]
    lines = '\n'.join([how] + lines)
    n_ch = '[8|9]'
    want = '.*updating environment:.*0 added, %s changed, 0 removed.*' % n_ch
    assert re.match(want, status, flags) is not None, lines
    want = ('.*executed 1 out of %s.*after excluding %s files.*based on MD5.*'
            % (out_of, excluding))
    assert re.match(want, status, flags) is not None
    got_stale = len(new_app.config.sphinx_gallery_conf['stale_examples'])
    assert got_stale == n_stale
    assert op.isfile(op.join(new_app.outdir, '_images',
                             'sphx_glr_plot_numpy_matplotlib_001.png'))

    generated_modules_1 = sorted(
        op.join(new_app.srcdir, 'gen_modules', f)
        for f in os.listdir(op.join(new_app.srcdir, 'gen_modules'))
        if op.isfile(op.join(new_app.srcdir, 'gen_modules', f)))
    generated_backrefs_1 = sorted(
        op.join(new_app.srcdir, 'gen_modules', 'backreferences', f)
        for f in os.listdir(op.join(new_app.srcdir, 'gen_modules',
                                    'backreferences')))
    generated_rst_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.rst'))
    generated_pickle_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.pickle'))
    copied_py_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.py'))
    copied_ipy_1 = sorted(
        op.join(new_app.srcdir, 'auto_examples', f)
        for f in os.listdir(op.join(new_app.srcdir, 'auto_examples'))
        if f.endswith('.ipynb'))

    # mtimes for modules
    _assert_mtimes(generated_modules_0, generated_modules_1)

    # mtimes for backrefs (gh-394)
    _assert_mtimes(generated_backrefs_0, generated_backrefs_1)

    # generated RST files
    different = ('plot_numpy_matplotlib',)
    ignore = (
        # this one should almost always be different, but in case we
        # get extremely unlucky and have identical run times
        # on the one script above that changes...
        'sg_execution_times',
        # this one will not change even though it was retried
        'plot_future_imports_broken',
        'plot_scraper_broken',
    )
    # not reliable on Windows and one Ubuntu run
    bad = sys.platform.startswith('win') or os.getenv('BAD_MTIME', '0') == '1'
    if not bad:
        _assert_mtimes(generated_rst_0, generated_rst_1, different, ignore)

        # mtimes for pickles
        use_different = () if how == 'run_stale' else different
        _assert_mtimes(generated_pickle_0, generated_pickle_1, ignore=ignore)

        # mtimes for .py files (gh-395)
        _assert_mtimes(copied_py_0, copied_py_1, different=use_different)

        # mtimes for .ipynb files
        _assert_mtimes(copied_ipy_0, copied_ipy_1, different=use_different)


@pytest.mark.parametrize('name, want', [
    ('future/plot_future_imports_broken',
     '.*RuntimeError.*Forcing this example to fail on Python 3.*'),
    ('plot_scraper_broken',
     '.*ValueError.*zero-size array to reduction.*'),
])
def test_error_messages(sphinx_app, name, want):
    """Test that informative error messages are added."""
    src_dir = sphinx_app.srcdir
    example_rst = op.join(src_dir, 'auto_examples', name + '.rst')
    with codecs.open(example_rst, 'r', 'utf-8') as fid:
        rst = fid.read()
    print(rst)
    rst = rst.replace('\n', ' ')
    assert re.match(want, rst) is not None


def test_alt_text_image(sphinx_app):
    """Test alt text for matplotlib images in html and rst"""
    out_dir = sphinx_app.outdir
    src_dir = sphinx_app.srcdir
    # alt text is fig titles, rst
    example_rst = op.join(src_dir, 'auto_examples', 'plot_matplotlib_alt.rst')
    with codecs.open(example_rst, 'r', 'utf-8') as fid:
        rst = fid.read()
    # suptitle and axes titles
    assert ':alt: This is a sup title, subplot 1, subplot 2' in rst
    # multiple titles
    assert ':alt: Left Title, Center Title, Right Title' in rst

    # no fig title - alt text is file name, rst
    example_rst = op.join(src_dir, 'auto_examples',
                          'plot_numpy_matplotlib.rst')
    with codecs.open(example_rst, 'r', 'utf-8') as fid:
        rst = fid.read()
    assert ':alt: plot numpy matplotlib' in rst
    # html
    example_html = op.join(out_dir, 'auto_examples',
                           'plot_numpy_matplotlib.html')
    with codecs.open(example_html, 'r', 'utf-8') as fid:
        html = fid.read()
    assert 'alt="plot numpy matplotlib"' in html


def test_alt_text_thumbnail(sphinx_app):
    """Test alt text for thumbnail in html and rst."""
    out_dir = sphinx_app.outdir
    src_dir = sphinx_app.srcdir
    # check gallery index thumbnail, html
    generated_examples_index = op.join(out_dir, 'auto_examples', 'index.html')
    with codecs.open(generated_examples_index, 'r', 'utf-8') as fid:
        html = fid.read()
    assert 'alt="&quot;SVG&quot;:-`graphics_`"' in html
    # check backreferences thumbnail, html
    backref_html = op.join(out_dir, 'gen_modules',
                           'sphinx_gallery.backreferences.html')
    with codecs.open(backref_html, 'r', 'utf-8') as fid:
        html = fid.read()
    assert 'alt="Link to other packages"' in html
    # check gallery index thumbnail, rst
    generated_examples_index = op.join(src_dir, 'auto_examples',
                                       'index.rst')
    with codecs.open(generated_examples_index, 'r', 'utf-8') as fid:
        rst = fid.read()
    assert ':alt: Trivial module to provide a value for plot_numpy_matplotlib.py' in rst  # noqa: E501


def test_backreference_labels(sphinx_app):
    """Tests that backreference labels work."""
    src_dir = sphinx_app.srcdir
    out_dir = sphinx_app.outdir
    # Test backreference label
    backref_rst = op.join(src_dir, 'gen_modules',
                          'sphinx_gallery.backreferences.rst')
    with codecs.open(backref_rst, 'r', 'utf-8') as fid:
        rst = fid.read()
    label = '.. _sphx_glr_backref_sphinx_gallery.backreferences.identify_names:'  # noqa: E501
    assert label in rst
    # Test html link
    index_html = op.join(out_dir, 'index.html')
    with codecs.open(index_html, 'r', 'utf-8') as fid:
        html = fid.read()
    link = 'href="gen_modules/sphinx_gallery.backreferences.html#sphx-glr-backref-sphinx-gallery-backreferences-identify-names">'  # noqa: E501
    assert link in html


@pytest.mark.parametrize(
    'test, nlines, filenamesortkey', [
        # first example, no heading
        ('Test 1-N', 6, False),
        # first example, default heading, default level
        ('Test 1-D-D', 8, False),
        # first example, default heading, custom level
        ('Test 1-D-C', 8, False),
        # first example, custom heading, default level
        ('Test 1-C-D', 9, False),
        # both examples, no heading
        ('Test 2-N', 11, True),
        # both examples, default heading, default level
        ('Test 2-D-D', 14, True),
        # both examples, custom heading, custom level
        ('Test 2-C-C', 15, True),
    ]
)
def test_minigallery_directive(sphinx_app, test, nlines, filenamesortkey):
    """Tests the functionality of the minigallery directive."""
    out_dir = sphinx_app.outdir
    minigallery_html = op.join(out_dir, 'minigallery.html')
    with codecs.open(minigallery_html, 'r', 'utf-8') as fid:
        lines = fid.readlines()

    # Regular expressions for matching
    any_heading = re.compile(r'<h([1-6])>.+<\/h\1>')
    explicitorder_example = re.compile(r'(?s)<img .+'
                                       r'(plot_second_future_imports).+'
                                       r'auto_examples\/\1\.html')
    filenamesortkey_example = re.compile(r'(?s)<img .+'
                                         r'(plot_numpy_matplotlib).+'
                                         r'auto_examples\/\1\.html')
    # Heading strings
    heading_str = {
        "Test 1-N": None,
        "Test 1-D-D": r'<h2>Examples using .+ExplicitOrder.+<\/h2>',
        "Test 1-D-C": r'<h3>Examples using .+ExplicitOrder.+<\/h3>',
        "Test 1-C-D": r'<h2>This is a custom heading.*<\/h2>',
        "Test 2-N": None,
        "Test 2-D-D": r'<h2>Examples using one of multiple objects.*<\/h2>',
        "Test 2-C-C": r'<h1>This is a different custom heading.*<\/h1>'
    }

    for i in range(len(lines)):
        if test in lines[i]:
            text = ''.join(lines[i:i+nlines])
            print(f'{test}: {text}')
            # Check headings
            if heading_str[test]:
                heading = re.compile(heading_str[test])
                assert heading.search(text) is not None
            else:
                # Confirm there isn't a heading
                assert any_heading.search(text) is None

            # Check for examples
            assert explicitorder_example.search(text) is not None
            if filenamesortkey:
                assert filenamesortkey_example.search(text) is not None
            else:
                assert filenamesortkey_example.search(text) is None


def test_matplotlib_warning_filter(sphinx_app):
    """Test Matplotlib agg warning is removed."""
    out_dir = sphinx_app.outdir
    example_html = op.join(out_dir, 'auto_examples',
                           'plot_matplotlib_alt.html')
    with codecs.open(example_html, 'r', 'utf-8') as fid:
        html = fid.read()
    warning = ('Matplotlib is currently using agg, which is a'
               ' non-GUI backend, so cannot show the figure.')
    assert warning not in html


def test_jupyter_notebook_pandoc(sphinx_app):
    """Test using pypandoc."""
    src_dir = sphinx_app.srcdir
    fname = op.join(src_dir, 'auto_examples', 'plot_numpy_matplotlib.ipynb')
    with codecs.open(fname, 'r', 'utf-8') as fid:
        md = fid.read()

    md_sg = r"Use :mod:`sphinx_gallery` to link to other packages, like\n:mod:`numpy`, :mod:`matplotlib.colors`, and :mod:`matplotlib.pyplot`."  # noqa
    md_pandoc = r'Use `sphinx_gallery`{.interpreted-text role=\"mod\"} to link to other\npackages, like `numpy`{.interpreted-text role=\"mod\"},\n`matplotlib.colors`{.interpreted-text role=\"mod\"}, and\n`matplotlib.pyplot`{.interpreted-text role=\"mod\"}.'  # noqa

    if any(_has_pypandoc()):
        assert md_pandoc in md
    else:
        assert md_sg in md


def test_md5_hash(sphinx_app):
    """Test MD5 hashing."""
    src_dir = sphinx_app.srcdir
    fname = op.join(src_dir, 'auto_examples', 'plot_log.py.md5')
    expected_md5 = '0edc2de97f96f3b55f8b4a21994931a8'
    with open(fname) as md5_file:
        actual_md5 = md5_file.read()

    assert actual_md5 == expected_md5


def test_binder_logo_exists(sphinx_app):
    """Test that the binder logo path is correct."""
    root = op.join(sphinx_app.outdir, 'auto_examples')
    with codecs.open(op.join(root, 'plot_svg.html'), 'r', 'utf-8') as fid:
        html = fid.read()
    path = re.match(r'.*<img alt="Launch binder" src="(.*)" width=.*\/>.*',
                    html, re.DOTALL)
    assert path is not None
    path = path.groups()[0]
    img_fname = op.abspath(op.join(root, path))
    assert 'binder_badge_logo' in img_fname  # can have numbers appended
    assert op.isfile(img_fname)
    assert 'https://mybinder.org/v2/gh/sphinx-gallery/sphinx-gallery.github.io/master?urlpath=lab/tree/notebooks/auto_examples/plot_svg.ipynb' in html  # noqa: E501


def test_defer_figures(sphinx_app):
    """Test the deferring of figures."""
    root = op.join(sphinx_app.outdir, 'auto_examples')
    fname = op.join(root, 'plot_defer_figures.html')
    with codecs.open(fname, 'r', 'utf-8') as fid:
        html = fid.read()

    # The example has two code blocks with plotting commands, but the first
    # block has the flag ``sphinx_gallery_defer_figures``.  Thus, there should
    # be only one image, not two, in the output.
    assert '../_images/sphx_glr_plot_defer_figures_001.png' in html
    assert '../_images/sphx_glr_plot_defer_figures_002.png' not in html


def test_no_dummy_image(sphinx_app):
    """Test that sphinx_gallery_dummy_images are NOT created (when executable
    is True)."""
    img1 = op.join(sphinx_app.srcdir, 'auto_examples', 'images',
                   'sphx_glr_plot_repr_001.png')
    img2 = op.join(sphinx_app.srcdir, 'auto_examples', 'images',
                   'sphx_glr_plot_repr_002.png')
    assert not op.isfile(img1)
    assert not op.isfile(img2)
