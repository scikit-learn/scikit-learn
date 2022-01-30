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
import importlib
import tempfile
import re
import os
import shutil
from unittest import mock
import zipfile
import codeop

import pytest

from sphinx.errors import ExtensionError
import sphinx_gallery.gen_rst as sg
from sphinx_gallery import downloads
from sphinx_gallery.gen_gallery import generate_dir_rst
from sphinx_gallery.scrapers import ImagePathIterator, figure_rst

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
    '# sphinx_gallery_thumbnail_number = 1',
    '# sphinx_gallery_defer_figures',
    '# and now comes the module code',
    'import logging',
    'import sys',
    'from warnings import warn',
    'x, y = 1, 2',
    'print(u"Óscar output") # need some code output',
    'logger = logging.getLogger()',
    'logger.setLevel(logging.INFO)',
    'lh = logging.StreamHandler(sys.stdout)',
    'lh.setFormatter(logging.Formatter("log:%(message)s"))',
    'logger.addHandler(lh)',
    'logger.info(u"Óscar")',
    'print(r"$\\langle n_\\uparrow n_\\downarrow \\rangle$")',
    'warn("WarningsAbound", RuntimeWarning)',
]


def test_split_code_and_text_blocks():
    """Test if a known example file gets properly split"""

    file_conf, blocks = sg.split_code_and_text_blocks(
        'examples/no_output/just_code.py')

    assert file_conf == {}
    assert blocks[0][0] == 'text'
    assert blocks[1][0] == 'code'


def test_bug_cases_of_notebook_syntax():
    """Test the known requirements of supported syntax in the
    notebook styled comments. Use both '#'s' and '#%%' as cell
    separators"""

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


def test_final_rst_last_word(tmpdir):
    """Tests last word in final rst block included as text"""
    filename = str(tmpdir.join('temp.py'))
    with open(filename, 'w') as f:
        f.write('\n'.join(['"Docstring"',
                           '# comment only code block',
                           '#%%',
                           '# Include this whole sentence.']))

    file_conf, result = sg.split_code_and_text_blocks(f.name)

    assert file_conf == {}
    expected_result = [
        ('text', 'Docstring', 1),
        ('code', '# comment only code block\n', 2),
        ('text', 'Include this whole sentence.', 4)]
    assert result == expected_result


def test_rst_block_after_docstring(gallery_conf, tmpdir):
    """Assert there is a blank line between the docstring and rst blocks."""
    filename = str(tmpdir.join('temp.py'))
    with open(filename, 'w') as f:
        f.write('\n'.join(['"Docstring"',
                           '####################',
                           '# Paragraph 1',
                           '# is long.',
                           '',
                           '#%%',
                           '# Paragraph 2',
                           '',
                           '# %%',
                           '# Paragraph 3',
                           '']))
    file_conf, blocks = sg.split_code_and_text_blocks(filename)

    assert file_conf == {}
    assert len(blocks) == 4
    assert blocks[0][0] == 'text'
    assert blocks[1][0] == 'text'
    assert blocks[2][0] == 'text'
    assert blocks[3][0] == 'text'

    script_vars = {'execute_script': ''}

    output_blocks, time_elapsed = sg.execute_script(
        blocks, script_vars, gallery_conf)

    example_rst = sg.rst_blocks(blocks, output_blocks, file_conf, gallery_conf)
    want_rst = """\
Docstring

.. GENERATED FROM PYTHON SOURCE LINES 3-5

Paragraph 1
is long.

.. GENERATED FROM PYTHON SOURCE LINES 7-8

Paragraph 2

.. GENERATED FROM PYTHON SOURCE LINES 10-11

Paragraph 3

"""
    assert example_rst == want_rst


def test_rst_empty_code_block(gallery_conf, tmpdir):
    """Test that we can "execute" a code block containing only comments."""
    gallery_conf.update(image_scrapers=())
    filename = str(tmpdir.join('temp.py'))
    with open(filename, 'w') as f:
        f.write('\n'.join(['"Docstring"',
                           '####################',
                           '# Paragraph 1',
                           '',
                           '# just a comment'
                           '']))
    file_conf, blocks = sg.split_code_and_text_blocks(filename)

    assert file_conf == {}
    assert len(blocks) == 3
    assert blocks[0][0] == 'text'
    assert blocks[1][0] == 'text'
    assert blocks[2][0] == 'code'

    gallery_conf['abort_on_example_error'] = True
    script_vars = dict(execute_script=True, src_file=filename,
                       image_path_iterator=[], target_file=filename)

    output_blocks, time_elapsed = sg.execute_script(
        blocks, script_vars, gallery_conf)

    example_rst = sg.rst_blocks(blocks, output_blocks, file_conf, gallery_conf)
    want_rst = """\
Docstring

.. GENERATED FROM PYTHON SOURCE LINES 3-4

Paragraph 1

.. GENERATED FROM PYTHON SOURCE LINES 4-5

.. code-block:: python


    # just a comment"""
    assert example_rst.rstrip('\n') == want_rst


def test_script_vars_globals(gallery_conf, tmpdir):
    """Assert the global vars get stored."""
    gallery_conf.update(image_scrapers=())
    filename = str(tmpdir.join('temp.py'))
    with open(filename, 'w') as f:
        f.write("""
'''
My example
----------

This is it.
'''
a = 1.
b = 'foo'
""")
    file_conf, blocks = sg.split_code_and_text_blocks(filename)
    assert len(blocks) == 2
    assert blocks[0][0] == 'text'
    assert blocks[1][0] == 'code'
    assert file_conf == {}
    script_vars = {'execute_script': True, 'src_file': filename,
                   'image_path_iterator': [],
                   'target_file': filename}
    output_blocks, time_elapsed = sg.execute_script(
        blocks, script_vars, gallery_conf)
    assert 'example_globals' in script_vars
    assert script_vars['example_globals']['a'] == 1.
    assert script_vars['example_globals']['b'] == 'foo'


def test_codestr2rst():
    """Test the correct translation of a code block into rst."""
    output = sg.codestr2rst('print("hello world")')
    reference = """\
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

    # SG incorrectly grabbing description when a label is defined (gh-232)
    intro_label, title_label = sg.extract_intro_and_title(
        '<string>', '\n'.join(['.. my_label', ''] + CONTENT[1:10]))
    assert intro_label == intro
    assert title_label == title

    intro_whitespace, title_whitespace = sg.extract_intro_and_title(
        '<string>', '\n'.join(CONTENT[1:4] + [''] + CONTENT[5:10]))
    assert intro_whitespace == intro
    assert title_whitespace == title

    # Make example title optional (gh-222)
    intro, title = sg.extract_intro_and_title('<string>', 'Title\n-----')
    assert intro == title == 'Title'

    # Title beginning with a space (gh-356)
    intro, title = sg.extract_intro_and_title('filename',
                                              '^^^^^\n   Title  two  \n^^^^^')
    assert intro == title == 'Title  two'

    # Title with punctuation (gh-517)
    intro, title = sg.extract_intro_and_title('<string>',
                                              '    ------------\n"-`Header"-with:; `punct` mark\'s\n----------------')  # noqa: E501
    assert title == '"-`Header"-with:; `punct` mark\'s'

    # Long intro paragraph gets shortened
    intro_paragraph = '\n'.join(['this is one line' for _ in range(10)])
    intro, _ = sg.extract_intro_and_title(
        'filename',
        'Title\n-----\n\n' + intro_paragraph)
    assert len(intro_paragraph) > 100
    assert len(intro) < 100
    assert intro.endswith('...')
    assert intro_paragraph.replace('\n', ' ')[:95] == intro[:95]

    # Errors
    with pytest.raises(ExtensionError, match='should have a header'):
        sg.extract_intro_and_title('<string>', '')  # no title
    with pytest.raises(ExtensionError, match='Could not find a title'):
        sg.extract_intro_and_title('<string>', '=====')  # no real title


@pytest.mark.parametrize(
    'mode,expected_md5',
    (['b', 'a546be453c8f436e744838a4801bd3a0'],
     ['t', 'ea8a570e9f3afc0a7c3f2a17a48b8047']))
def test_md5sums(mode, expected_md5):
    """Test md5sum check functions work on know file content."""
    file_content = b'Local test\r\n'
    with tempfile.NamedTemporaryFile('wb', delete=False) as f:
        f.write(file_content)
    try:
        file_md5 = sg.get_md5sum(f.name, mode)
        # verify correct md5sum
        assert file_md5 == expected_md5
        # False because is a new file
        assert not sg.md5sum_is_current(f.name)
        # Write md5sum to file to check is current
        with open(f.name + '.md5', 'w') as file_checksum:
            file_checksum.write(file_md5)
        try:
            assert sg.md5sum_is_current(f.name, mode)
        finally:
            os.remove(f.name + '.md5')
    finally:
        os.remove(f.name)


@pytest.mark.parametrize('failing_code, want', [
    (CONTENT + ['#' * 79, 'First_test_fail', '#' * 79, 'second_fail'],
     'not defined'),
    (CONTENT + ['#' * 79, 'input("foo")', '#' * 79, 'second_fail'],
     'Cannot use input'),
    (CONTENT + ['#' * 79, 'bad syntax', '#' * 79, 'second_fail'],
     'invalid syntax'),
])
def test_fail_example(gallery_conf, failing_code, want,
                      log_collector, req_pil):
    """Test that failing examples are only executed until failing block."""
    gallery_conf.update(image_scrapers=(), reset_modules=())
    gallery_conf.update(filename_pattern='raise.py')

    with codecs.open(os.path.join(gallery_conf['examples_dir'], 'raise.py'),
                     mode='w', encoding='utf-8') as f:
        f.write('\n'.join(failing_code))

    sg.generate_file_rst('raise.py', gallery_conf['gallery_dir'],
                         gallery_conf['examples_dir'], gallery_conf)
    assert len(log_collector.calls['warning']) == 1
    msg = log_collector.calls['warning'][0].args[2]
    assert want in msg
    assert 'gen_gallery' not in msg
    # can only check that gen_rst is removed on non-input ones
    if 'Cannot use input' not in msg:
        assert 'gen_rst' not in msg
    assert '_check_input' not in msg

    # read rst file and check if it contains traceback output

    with codecs.open(os.path.join(gallery_conf['gallery_dir'], 'raise.rst'),
                     mode='r', encoding='utf-8') as f:
        ex_failing_blocks = f.read().count('pytb')
        assert ex_failing_blocks != 0, 'Did not run into errors in bad code'
        assert ex_failing_blocks <= 1, \
            'Did not stop executing script after error'


def _generate_rst(gallery_conf, fname, content):
    """Return the rST text of a given example content.

    This writes a file gallery_conf['examples_dir']/fname with *content*,
    creates the corresponding rst file by running generate_file_rst() and
    returns the generated rST code.

    Parameters
    ----------
    gallery_conf
        A gallery_conf as created by the gallery_conf fixture.
    fname : str
        A filename; e.g. 'test.py'. This is relative to
        gallery_conf['examples_dir']
    content : str
        The content of fname.

    Returns
    -------
    rst : str
        The generated rST code.
    """
    with codecs.open(os.path.join(gallery_conf['examples_dir'], fname),
                     mode='w', encoding='utf-8') as f:
        f.write('\n'.join(content))
    # generate rst file
    sg.generate_file_rst(fname, gallery_conf['gallery_dir'],
                         gallery_conf['examples_dir'], gallery_conf)
    # read rst file and check if it contains code output
    rst_fname = os.path.splitext(fname)[0] + '.rst'
    with codecs.open(os.path.join(gallery_conf['gallery_dir'], rst_fname),
                     mode='r', encoding='utf-8') as f:
        rst = f.read()
    return rst


ALPHA_CONTENT = '''
"""
Make a plot
===========

Plot.
"""
import matplotlib.pyplot as plt
plt.plot([0, 1], [0, 1])
'''.split('\n')


def _alpha_mpl_scraper(block, block_vars, gallery_conf):
    import matplotlib.pyplot as plt
    image_path_iterator = block_vars['image_path_iterator']
    image_paths = list()
    for fig_num, image_path in zip(plt.get_fignums(), image_path_iterator):
        fig = plt.figure(fig_num)
        assert image_path.endswith('.png')
        # use format that does not support alpha
        image_path = image_path[:-3] + 'jpg'
        fig.savefig(image_path)
        image_paths.append(image_path)
    plt.close('all')
    return figure_rst(image_paths, gallery_conf['src_dir'])


def test_custom_scraper_thumbnail_alpha(gallery_conf, req_mpl_jpg):
    """Test that thumbnails without an alpha channel work w/custom scraper."""
    gallery_conf['image_scrapers'] = [_alpha_mpl_scraper]
    rst = _generate_rst(gallery_conf, 'plot_test.py', ALPHA_CONTENT)
    assert '.jpg' in rst


def test_remove_config_comments(gallery_conf, req_pil):
    """Test the gallery_conf['remove_config_comments'] setting."""
    rst = _generate_rst(gallery_conf, 'test.py', CONTENT)
    assert '# sphinx_gallery_thumbnail_number = 1' in rst
    assert '# sphinx_gallery_defer_figures' in rst

    gallery_conf['remove_config_comments'] = True
    rst = _generate_rst(gallery_conf, 'test.py', CONTENT)
    assert '# sphinx_gallery_thumbnail_number = 1' not in rst
    assert '# sphinx_gallery_defer_figures' not in rst


def test_dummy_image_error(gallery_conf, req_pil):
    """Test correct error is raised if int not provided to
    sphinx_gallery_dummy_images."""
    content_image = CONTENT + ["# sphinx_gallery_dummy_images=False", ]
    msg = "sphinx_gallery_dummy_images setting is not a number"
    with pytest.raises(ExtensionError, match=msg):
        _generate_rst(gallery_conf, 'test.py', content_image)


def test_final_empty_block(gallery_conf, req_pil):
    """Test empty final block is removed. Empty final block can occur after
    sole config comment is removed from final block."""
    content_block = \
        CONTENT + ['# %%', '', '# sphinx_gallery_line_numbers = True']
    gallery_conf['remove_config_comments'] = True
    rst = _generate_rst(gallery_conf, 'test.py', content_block)
    want = "RuntimeWarning)\n\n\n.. rst-class:: sphx-glr-timing"
    assert want in rst


def test_download_link_note_only_html(gallery_conf, req_pil):
    """Test html only directive for download_link."""
    rst = _generate_rst(gallery_conf, 'test.py', CONTENT)
    download_link_note = (".. only:: html\n\n"
                          "    .. note::\n"
                          "        :class: sphx-glr-download-link-note\n\n"
                          )
    assert download_link_note in rst


def test_download_link_classes(gallery_conf, req_pil):
    """Test classes for download links."""
    rst = _generate_rst(gallery_conf, 'test.py', CONTENT)
    for kind in ('python', 'jupyter'):
        assert 'sphx-glr-download sphx-glr-download-' + kind in rst


@pytest.mark.parametrize('ext', ('.txt', '.rst', '.bad'))
def test_gen_dir_rst(gallery_conf, fakesphinxapp, ext):
    """Test gen_dir_rst."""
    print(os.listdir(gallery_conf['examples_dir']))
    fname_readme = os.path.join(gallery_conf['src_dir'], 'README.txt')
    with open(fname_readme, 'wb') as fid:
        fid.write(u"Testing\n=======\n\nÓscar here.".encode('utf-8'))
    fname_out = os.path.splitext(fname_readme)[0] + ext
    if fname_readme != fname_out:
        shutil.move(fname_readme, fname_out)
    args = (gallery_conf['src_dir'], gallery_conf['gallery_dir'],
            gallery_conf, [])
    if ext == '.bad':  # not found with correct ext
        with pytest.raises(ExtensionError, match='does not have a README'):
            generate_dir_rst(*args)
    else:
        out = generate_dir_rst(*args)
        assert u"Óscar here" in out[0]


def test_pattern_matching(gallery_conf, log_collector, req_pil):
    """Test if only examples matching pattern are executed."""
    gallery_conf.update(image_scrapers=(), reset_modules=())
    gallery_conf.update(filename_pattern=re.escape(os.sep) + 'plot_0')

    code_output = ('\n Out:\n\n .. code-block:: none\n'
                   '\n'
                   '    Óscar output\n'
                   '    log:Óscar\n'
                   '    $\\langle n_\\uparrow n_\\downarrow \\rangle$'
                   )
    warn_output = 'RuntimeWarning: WarningsAbound'
    # create three files in tempdir (only one matches the pattern)
    fnames = ['plot_0.py', 'plot_1.py', 'plot_2.py']
    for fname in fnames:
        rst = _generate_rst(gallery_conf, fname, CONTENT)
        rst_fname = os.path.splitext(fname)[0] + '.rst'
        if re.search(gallery_conf['filename_pattern'],
                     os.path.join(gallery_conf['gallery_dir'], rst_fname)):
            assert code_output in rst
            assert warn_output in rst
        else:
            assert code_output not in rst
            assert warn_output not in rst


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


@pytest.mark.parametrize('test_str', [
    "# sphinx_gallery_thumbnail_path= '_static/demo.png'",
    "# sphinx_gallery_thumbnail_path='_static/demo.png'",
    "#sphinx_gallery_thumbnail_path = '_static/demo.png'",
    "    # sphinx_gallery_thumbnail_path='_static/demo.png'"])
def test_thumbnail_path(test_str):
    # which plot to show as the thumbnail image
    with tempfile.NamedTemporaryFile('w', delete=False) as f:
        f.write('\n'.join(['"Docstring"',
                           test_str]))
    try:
        file_conf, blocks = sg.split_code_and_text_blocks(f.name)
    finally:
        os.remove(f.name)
    assert file_conf == {'thumbnail_path': '_static/demo.png'}


def test_zip_notebooks(gallery_conf):
    """Test generated zipfiles are not corrupt."""
    gallery_conf.update(
        examples_dir=os.path.join(gallery_conf['src_dir'], 'examples'))
    shutil.copytree(
        os.path.join(os.path.dirname(__file__), 'tinybuild', 'examples'),
        gallery_conf['examples_dir'])
    examples = downloads.list_downloadable_sources(
        gallery_conf['examples_dir'])
    zipfilepath = downloads.python_zip(examples, gallery_conf['gallery_dir'])
    zipf = zipfile.ZipFile(zipfilepath)
    check = zipf.testzip()
    assert not check, "Bad file in zipfile: {0}".format(check)


def test_rst_example(gallery_conf):
    """Test generated rst file includes the correct paths for binder."""
    gallery_conf.update(binder={'org': 'sphinx-gallery',
                                'repo': 'sphinx-gallery.github.io',
                                'binderhub_url': 'https://mybinder.org',
                                'branch': 'master',
                                'dependencies': './binder/requirements.txt',
                                'notebooks_dir': 'notebooks',
                                'use_jupyter_lab': True,
                                },
                        gallery_dirs=None)

    example_file = os.path.join(gallery_conf['gallery_dir'], "plot.py")
    sg.save_rst_example("example_rst", example_file, 0, 0, gallery_conf)

    test_file = re.sub(r'\.py$', '.rst', example_file)
    with codecs.open(test_file) as f:
        rst = f.read()

    assert "lab/tree/notebooks/plot.ipy" in rst

    # CSS classes
    assert "rst-class:: sphx-glr-signature" in rst
    assert "rst-class:: sphx-glr-timing" in rst


@pytest.fixture(scope='function')
def script_vars(tmpdir):
    fake_main = importlib.util.module_from_spec(
        importlib.util.spec_from_loader('__main__', None))
    fake_main.__dict__.update({'__doc__': ''})
    script_vars = {
        "execute_script": True,
        "image_path_iterator": ImagePathIterator(str(tmpdir.join("temp.png"))),
        "src_file": __file__,
        "memory_delta": [],
        "fake_main": fake_main,
    }
    return script_vars


def test_output_indentation(gallery_conf, script_vars):
    """Test whether indentation of code output is retained."""
    gallery_conf.update(image_scrapers=())
    compiler = codeop.Compile()

    test_string = r"\n".join([
        "  A B",
        "A 1 2",
        "B 3 4"
    ])
    code = "print('" + test_string + "')"
    code_block = ("code", code, 1)
    output = sg.execute_code_block(
        compiler, code_block, None, script_vars, gallery_conf
    )
    output_test_string = "\n".join(
        [line[4:] for line in output.strip().split("\n")[-3:]]
    )
    assert output_test_string == test_string.replace(r"\n", "\n")


def test_output_no_ansi(gallery_conf, script_vars):
    """Test ANSI characters are removed.

    See: https://en.wikipedia.org/wiki/ANSI_escape_code
    """
    gallery_conf.update(image_scrapers=())
    compiler = codeop.Compile()

    code = 'print("\033[94m0.25")'
    code_block = ("code", code, 1)
    output = sg.execute_code_block(
        compiler, code_block, None, script_vars, gallery_conf
    )
    output_test_string = "\n".join(
        [line[4:] for line in output.strip().split("\n")[-3:]]
    )

    assert output_test_string.split('\n')[-1] == "0.25"


def test_empty_output_box(gallery_conf, script_vars):
    """Tests that `print(__doc__)` doesn't produce an empty output box."""
    gallery_conf.update(image_scrapers=())
    compiler = codeop.Compile()

    code_block = ("code", "print(__doc__)", 1)

    output = sg.execute_code_block(
        compiler, code_block, None, script_vars, gallery_conf
    )
    assert output.isspace()


code_repr_only = """
class repr_only_class():

    def __init__(self):
        pass

    def __repr__(self):
        return "This is the __repr__"

class_inst = repr_only_class()
class_inst
"""

code_repr_and_html = """
class repr_and_html_class():
    def __init__(self):
        pass

    def __repr__(self):
        return "This is the __repr__"

    def _repr_html_(self):
        return "<div> This is the _repr_html_ div </div>"

class_inst = repr_and_html_class()
class_inst
"""

code_print_and_repr_and_html = """
print("print statement")

class repr_and_html_class():
    def __init__(self):
        pass

    def __repr__(self):
        return "This is the __repr__"

    def _repr_html_(self):
        return "<div> This is the _repr_html_ div </div>"

class_inst = repr_and_html_class()
class_inst
"""

code_plt = """
import matplotlib.pyplot as plt
fig = plt.figure()
plt.close('all')
fig
"""

html_out = """.. raw:: html

    <div class="output_subarea output_html rendered_html output_result">
    <div> This is the _repr_html_ div </div>
    </div>
    <br />
    <br />"""

text_above_html = """Out:

 .. code-block:: none

    print statement


"""


def _clean_output(output):
    is_text = '.. rst-class:: sphx-glr-script-out' in output

    is_html = '.. raw:: html' in output

    if output.isspace():
        return ''
    elif is_text and is_html:
        output_test_string = "\n".join(output.strip().split("\n")[2:])
        return output_test_string.strip()
    elif is_text:
        output_test_string = "\n".join(
            [line[4:] for line in output.strip().split("\n")[6:]])
        return output_test_string.strip()
    elif is_html:
        output_test_string = "\n".join(output.strip().split("\n"))
        return output_test_string


@pytest.mark.parametrize('capture_repr, code, expected_out', [
    pytest.param(tuple(), 'a=2\nb=3', '', id='assign,()'),
    pytest.param(tuple(), 'a=2\na', '', id='var,()'),
    pytest.param(tuple(), 'a=2\nprint(a)', '2', id='print(var),()'),
    pytest.param(tuple(), 'print("hello")\na=2\na', 'hello',
                 id='print+var,()'),
    pytest.param(('__repr__',), 'a=2\nb=3', '', id='assign,(repr)'),
    pytest.param(('__repr__',), 'a=2\na', '2', id='var,(repr)'),
    pytest.param(('__repr__',), 'a=2\nprint(a)', '2', id='print(var),(repr)'),
    pytest.param(('__repr__',), 'print("hello")\na=2\na', 'hello\n\n2',
                 id='print+var,(repr)'),
    pytest.param(('__repr__',), code_repr_and_html, 'This is the __repr__',
                 id='repr_and_html,(repr)'),
    pytest.param(('__repr__',), code_print_and_repr_and_html,
                 'print statement\n\nThis is the __repr__',
                 id='print and repr_and_html,(repr)'),
    pytest.param(('_repr_html_',), code_repr_only, '', id='repr_only,(html)'),
    pytest.param(('_repr_html_',), code_repr_and_html, html_out,
                 id='repr_and_html,(html)'),
    pytest.param(('_repr_html_',), code_print_and_repr_and_html,
                 ''.join([text_above_html, html_out]),
                 id='print and repr_and_html,(html)'),
    pytest.param(('_repr_html_', '__repr__'), code_repr_and_html, html_out,
                 id='repr_and_html,(html,repr)'),
    pytest.param(('__repr__', '_repr_html_'), code_repr_and_html,
                 'This is the __repr__', id='repr_and_html,(repr,html)'),
    pytest.param(('_repr_html_', '__repr__'), code_repr_only,
                 'This is the __repr__', id='repr_only,(html,repr)'),
    pytest.param(('_repr_html_',), code_plt, '', id='html_none'),
])
def test_capture_repr(gallery_conf, capture_repr, code, expected_out,
                      req_mpl, req_pil, script_vars):
    """Tests output capturing with various capture_repr settings."""
    compiler = codeop.Compile()
    code_block = ('code', code, 1)
    gallery_conf['capture_repr'] = capture_repr
    output = sg.execute_code_block(
        compiler, code_block, None, script_vars, gallery_conf
    )
    assert _clean_output(output) == expected_out


def test_ignore_repr_types(gallery_conf, req_mpl, req_pil, script_vars):
    """Tests output capturing with various capture_repr settings."""
    compiler = codeop.Compile()
    code_block = ('code', 'a=2\na', 1)
    gallery_conf['ignore_repr_types'] = r'int'
    output = sg.execute_code_block(
        compiler, code_block, None, script_vars, gallery_conf
    )
    assert _clean_output(output) == ''


@pytest.mark.parametrize(
    ('order', 'call_count'), [('before', 1), ('after', 1), ('both', 2)]
)
def test_reset_module_order_2_param(gallery_conf, order, call_count, req_pil):
    """Test that reset module with 2 parameters."""

    def cleanup_2_param(gallery_conf, fname):
        pass

    mock_reset_module = mock.create_autospec(cleanup_2_param)
    gallery_conf['reset_modules'] = (mock_reset_module,)
    gallery_conf['reset_modules_order'] = order
    _generate_rst(gallery_conf, 'plot_test.py', CONTENT)
    assert mock_reset_module.call_count == call_count


@pytest.mark.parametrize(
    ('order', 'call_count', 'expected_call_order'),
    [
        ('before', 1, ('before',)),
        ('after', 1, ('after',)),
        ('both', 2, ('before', 'after'))
    ]
)
def test_reset_module_order_3_param(gallery_conf, order, call_count,
                                    expected_call_order, req_pil):
    """Test reset module with 3 parameters."""

    def cleanup_3_param(gallery_conf, fname, when):
        pass

    mock_reset_module = mock.create_autospec(cleanup_3_param)
    gallery_conf['reset_modules'] = (mock_reset_module,)
    gallery_conf['reset_modules_order'] = order
    _generate_rst(gallery_conf, 'plot_test.py', CONTENT)
    assert mock_reset_module.call_count == call_count

    expected_calls = [
        mock.call(mock.ANY, mock.ANY, order) for order in expected_call_order
    ]
    mock_reset_module.assert_has_calls(expected_calls)


def test_reset_module_order_3_param_invalid_when(gallery_conf):
    """Test reset module with unknown 3rd parameter."""

    def cleanup_3_param(gallery_conf, fname, invalid):
        pass

    mock_reset_module = mock.create_autospec(cleanup_3_param)
    gallery_conf['reset_modules'] = (mock_reset_module,)
    gallery_conf['reset_modules_order'] = 'before'
    with pytest.raises(ValueError,
                       match=("3rd parameter in cleanup_3_param "
                              "function signature must be 'when'")):
        _generate_rst(gallery_conf, 'plot_test.py', CONTENT)
    assert mock_reset_module.call_count == 0


class TestLoggingTee:
    def setup(self):
        self.src_filename = 'source file name'
        self.tee = sg._LoggingTee(self.src_filename)
        self.output_file = self.tee.output

    def test_full_line(self, log_collector):
        # A full line is output immediately.
        self.tee.write('Output\n')
        self.tee.flush()
        assert self.output_file.getvalue() == 'Output\n'
        assert len(log_collector.calls['verbose']) == 2
        assert self.src_filename in log_collector.calls['verbose'][0].args
        assert 'Output' in log_collector.calls['verbose'][1].args

    def test_incomplete_line_with_flush(self, log_collector):
        # An incomplete line ...
        self.tee.write('Output')
        assert self.output_file.getvalue() == 'Output'
        assert len(log_collector.calls['verbose']) == 1
        assert self.src_filename in log_collector.calls['verbose'][0].args

        # ... should appear when flushed.
        self.tee.flush()
        assert len(log_collector.calls['verbose']) == 2
        assert 'Output' in log_collector.calls['verbose'][1].args

    def test_incomplete_line_with_more_output(self, log_collector):
        # An incomplete line ...
        self.tee.write('Output')
        assert self.output_file.getvalue() == 'Output'
        assert len(log_collector.calls['verbose']) == 1
        assert self.src_filename in log_collector.calls['verbose'][0].args

        # ... should appear when more data is written.
        self.tee.write('\nMore output\n')
        assert self.output_file.getvalue() == 'Output\nMore output\n'
        assert len(log_collector.calls['verbose']) == 3
        assert 'Output' in log_collector.calls['verbose'][1].args
        assert 'More output' in log_collector.calls['verbose'][2].args

    def test_multi_line(self, log_collector):
        self.tee.write('first line\rsecond line\nthird line')
        assert (self.output_file.getvalue() ==
                'first line\rsecond line\nthird line')
        verbose_calls = log_collector.calls['verbose']
        assert len(verbose_calls) == 3
        assert self.src_filename in verbose_calls[0].args
        assert 'first line' in verbose_calls[1].args
        assert 'second line' in verbose_calls[2].args
        assert self.tee.logger_buffer == 'third line'

    def test_isatty(self, monkeypatch):
        assert not self.tee.isatty()

        monkeypatch.setattr(self.tee.output, 'isatty', lambda: True)
        assert self.tee.isatty()


# TODO: test that broken thumbnail does appear when needed
# TODO: test that examples are executed after a no-plot and produce
#       the correct image in the thumbnail
