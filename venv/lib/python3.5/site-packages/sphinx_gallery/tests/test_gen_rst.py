# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
Testing the rst files generator
"""
from __future__ import (division, absolute_import, print_function,
                        unicode_literals)
import ast
import codecs
import copy
import tempfile
import re
import os
import shutil
import zipfile
import pytest

import sphinx_gallery.gen_rst as sg
from sphinx_gallery import gen_gallery, downloads
from sphinx_gallery.gen_gallery import generate_dir_rst
from sphinx_gallery.utils import _TempDir

# Need to import gen_rst before matplotlib.pyplot to set backend to 'Agg'
import matplotlib.pyplot as plt

CONTENT = [
    '"""',
    '================',
    'Docstring header',
    '================',
    '',
    'This is the description of the example',
    'which goes on and on, Óscar',
    '',
    '',
    'And this is a second paragraph',
    '"""',
    '',
    '# and now comes the module code',
    'import logging',
    'import sys',
    'x, y = 1, 2',
    'print(u"Óscar output") # need some code output',
    'logger = logging.getLogger()',
    'logger.setLevel(logging.INFO)',
    'lh = logging.StreamHandler(sys.stdout)',
    'lh.setFormatter(logging.Formatter("log:%(message)s"))',
    'logger.addHandler(lh)',
    'logger.info(u"Óscar")',
    'print(r"$\\langle n_\\uparrow n_\\downarrow \\rangle$")',
]


def test_split_code_and_text_blocks():
    """Test if a known example file gets properly split"""

    file_conf, blocks = sg.split_code_and_text_blocks(
        'examples/no_output/just_code.py')

    assert file_conf == {}
    assert blocks[0][0] == 'text'
    assert blocks[1][0] == 'code'


def test_bug_cases_of_notebook_syntax():
    """Test over the known requirements of supported syntax in the
    notebook styled comments"""

    with open('sphinx_gallery/tests/reference_parse.txt') as reference:
        ref_blocks = ast.literal_eval(reference.read())
        file_conf, blocks = sg.split_code_and_text_blocks(
            'tutorials/plot_parse.py')

        assert file_conf == {}
        assert blocks == ref_blocks


def test_direct_comment_after_docstring():
    # For more details see
    # https://github.com/sphinx-gallery/sphinx-gallery/pull/49
    with tempfile.NamedTemporaryFile('w', delete=False) as f:
        f.write('\n'.join(['"Docstring"',
                           '# and now comes the module code',
                           '# with a second line of comment',
                           'x, y = 1, 2',
                           '']))
    try:
        file_conf, result = sg.split_code_and_text_blocks(f.name)
    finally:
        os.remove(f.name)

    assert file_conf == {}
    expected_result = [
        ('text', 'Docstring', 1),
        ('code', '\n'.join(['# and now comes the module code',
                            '# with a second line of comment',
                            'x, y = 1, 2',
                            '']), 2)]
    assert result == expected_result


def test_codestr2rst():
    """Test the correct translation of a code block into rst"""
    output = sg.codestr2rst('print("hello world")')
    reference = """
.. code-block:: python

    print("hello world")"""
    assert reference == output


def test_extract_intro_and_title():
    intro, title = sg.extract_intro_and_title('<string>',
                                              '\n'.join(CONTENT[1:10]))
    assert title == 'Docstring header'
    assert 'Docstring' not in intro
    assert intro == 'This is the description of the example which goes on and on, Óscar'  # noqa
    assert 'second paragraph' not in intro


def test_md5sums():
    """Test md5sum check functions work on know file content"""

    with tempfile.NamedTemporaryFile('wb', delete=False) as f:
        f.write(b'Local test\n')
    try:
        file_md5 = sg.get_md5sum(f.name)
        # verify correct md5sum
        assert 'ea8a570e9f3afc0a7c3f2a17a48b8047' == file_md5
        # False because is a new file
        assert not sg.md5sum_is_current(f.name)
        # Write md5sum to file to check is current
        with open(f.name + '.md5', 'w') as file_checksum:
            file_checksum.write(file_md5)
        try:
            assert sg.md5sum_is_current(f.name)
        finally:
            os.remove(f.name + '.md5')
    finally:
        os.remove(f.name)


@pytest.fixture
def gallery_conf():
    """Sets up a test sphinx-gallery configuration"""

    gallery_conf = copy.deepcopy(gen_gallery.DEFAULT_GALLERY_CONF)
    gallery_conf.update(examples_dir=_TempDir(), gallery_dir=_TempDir())
    gallery_conf['src_dir'] = gallery_conf['gallery_dir']

    return gallery_conf


def test_fail_example(gallery_conf, log_collector):
    """Test that failing examples are only executed until failing block"""

    gallery_conf.update(filename_pattern='raise.py')

    failing_code = CONTENT + ['#' * 79,
                              'First_test_fail', '#' * 79, 'second_fail']

    with codecs.open(os.path.join(gallery_conf['examples_dir'], 'raise.py'),
                     mode='w', encoding='utf-8') as f:
        f.write('\n'.join(failing_code))

    sg.generate_file_rst('raise.py', gallery_conf['gallery_dir'],
                         gallery_conf['examples_dir'], gallery_conf)
    assert len(log_collector.calls['warning']) == 1
    assert 'not defined' in log_collector.calls['warning'][0].args[2]

    # read rst file and check if it contains traceback output

    with codecs.open(os.path.join(gallery_conf['gallery_dir'], 'raise.rst'),
                     mode='r', encoding='utf-8') as f:
        ex_failing_blocks = f.read().count('pytb')
        if ex_failing_blocks == 0:
            raise ValueError('Did not run into errors in bad code')
        elif ex_failing_blocks > 1:
            raise ValueError('Did not stop executing script after error')


def test_gen_dir_rst(gallery_conf, fakesphinxapp):
    """Test gen_dir_rst."""
    print(os.listdir(gallery_conf['examples_dir']))
    fname_readme = os.path.join(gallery_conf['src_dir'], 'README.txt')
    with open(fname_readme, 'wb') as fid:
        fid.write(u"Testing\n=======\n\nÓscar here.".encode('utf-8'))
    args = (gallery_conf['src_dir'], gallery_conf['gallery_dir'],
            gallery_conf, [])
    out = generate_dir_rst(*args)
    assert u"Óscar here" in out[0]


def test_pattern_matching(gallery_conf, log_collector):
    """Test if only examples matching pattern are executed"""

    gallery_conf.update(filename_pattern=re.escape(os.sep) + 'plot_0')

    code_output = ('\n Out::\n'
                   '\n'
                   '    Óscar output\n'
                   '    log:Óscar\n'
                   '    $\\langle n_\\uparrow n_\\downarrow \\rangle$\n\n'
                   )
    # create three files in tempdir (only one matches the pattern)
    fnames = ['plot_0.py', 'plot_1.py', 'plot_2.py']
    for fname in fnames:
        with codecs.open(os.path.join(gallery_conf['examples_dir'], fname),
                         mode='w', encoding='utf-8') as f:
            f.write('\n'.join(CONTENT))
        # generate rst file
        sg.generate_file_rst(fname, gallery_conf['gallery_dir'],
                             gallery_conf['examples_dir'], gallery_conf)
        # read rst file and check if it contains code output
        rst_fname = os.path.splitext(fname)[0] + '.rst'
        with codecs.open(os.path.join(gallery_conf['gallery_dir'], rst_fname),
                         mode='r', encoding='utf-8') as f:
            rst = f.read()
        if re.search(gallery_conf['filename_pattern'],
                     os.path.join(gallery_conf['gallery_dir'], rst_fname)):
            assert code_output in rst
        else:
            assert code_output not in rst


@pytest.mark.parametrize('test_str', [
    '# sphinx_gallery_thumbnail_number= 2',
    '# sphinx_gallery_thumbnail_number=2',
    '#sphinx_gallery_thumbnail_number = 2',
    '    # sphinx_gallery_thumbnail_number=2'])
def test_thumbnail_number(test_str):
    # which plot to show as the thumbnail image
    with tempfile.NamedTemporaryFile('w', delete=False) as f:
        f.write('\n'.join(['"Docstring"',
                           test_str]))
    try:
        file_conf, blocks = sg.split_code_and_text_blocks(f.name)
    finally:
        os.remove(f.name)
    assert file_conf == {'thumbnail_number': 2}


def test_save_figures(gallery_conf):
    """Test file naming when saving figures. Requires mayavi."""
    try:
        from mayavi import mlab
    except ImportError:
        raise pytest.skip('Mayavi not installed')
    mlab.options.offscreen = True

    gallery_conf.update(find_mayavi_figures=True)

    mlab.test_plot3d()
    plt.plot(1, 1)
    fname_template = os.path.join(gallery_conf['gallery_dir'], 'image{0}.png')
    image_rst, fig_num = sg.save_figures(fname_template, 0, gallery_conf)
    assert fig_num == 2
    assert '/image1.png' in image_rst
    assert '/image2.png' in image_rst

    mlab.test_plot3d()
    plt.plot(1, 1)
    image_rst, fig_num = sg.save_figures(fname_template, 2, gallery_conf)
    assert fig_num == 2
    assert '/image2.png' not in image_rst
    assert '/image3.png' in image_rst
    assert '/image4.png' in image_rst

    shutil.rmtree(gallery_conf['gallery_dir'])


def test_zip_notebooks(gallery_conf):
    """Test generated zipfiles are not corrupt"""
    gallery_conf.update(examples_dir='examples')
    examples = downloads.list_downloadable_sources(
        gallery_conf['examples_dir'])
    zipfilepath = downloads.python_zip(examples, gallery_conf['gallery_dir'])
    zipf = zipfile.ZipFile(zipfilepath)
    check = zipf.testzip()
    if check:
        raise OSError("Bad file in zipfile: {0}".format(check))


def test_figure_rst():
    """Testing rst of images"""
    figure_list = ['sphx_glr_plot_1.png']
    image_rst, fig_num = sg.figure_rst(figure_list, '.')
    single_image = """
.. image:: /sphx_glr_plot_1.png
    :align: center
"""
    assert image_rst == single_image
    assert fig_num == 1

    image_rst, fig_num = sg.figure_rst(figure_list + ['second.png'], '.')

    image_list_rst = """
.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /sphx_glr_plot_1.png
            :scale: 47

    *

      .. image:: /second.png
            :scale: 47
"""
    assert image_rst == image_list_rst
    assert fig_num == 2

    # test issue #229
    local_img = [os.path.join(os.getcwd(), 'third.png')]
    image_rst, fig_num = sg.figure_rst(local_img, '.')

    single_image = sg.SINGLE_IMAGE % "third.png"
    assert image_rst == single_image
    assert fig_num == 1


# TODO: test that broken thumbnail does appear when needed
# TODO: test that examples are not executed twice
# TODO: test that examples are executed after a no-plot and produce
#       the correct image in the thumbnail
