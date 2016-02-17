# -*- coding: utf-8 -*-
# Author: Óscar Nájera
# License: 3-clause BSD
"""
==================
RST file generator
==================

Generate the rst files for the examples by iterating over the python
example files.

Files that generate images should start with 'plot'

"""
from __future__ import division, print_function, absolute_import
from time import time
import ast
import hashlib
import os
import re
import shutil
import subprocess
import sys
import traceback
import warnings


# Try Python 2 first, otherwise load from Python 3
from textwrap import dedent
try:
    # textwrap indent only exists in python 3
    from textwrap import indent
except ImportError:
    def indent(text, prefix, predicate=None):
        """Adds 'prefix' to the beginning of selected lines in 'text'.

        If 'predicate' is provided, 'prefix' will only be added to the lines
        where 'predicate(line)' is True. If 'predicate' is not provided,
        it will default to adding 'prefix' to all non-empty lines that do not
        consist solely of whitespace characters.
        """
        if predicate is None:
            def predicate(line):
                return line.strip()

        def prefixed_lines():
            for line in text.splitlines(True):
                yield (prefix + line if predicate(line) else line)
        return ''.join(prefixed_lines())

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

try:
    # make sure that the Agg backend is set before importing any
    # matplotlib
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    # this script can be imported by nosetest to find tests to run: we should
    # not impose the matplotlib requirement in that case.
    pass

from . import glr_path_static
from .backreferences import write_backreferences, _thumbnail_div
from .notebook import Notebook

try:
    basestring
except NameError:
    basestring = str


###############################################################################


class Tee(object):
    """A tee object to redirect streams to multiple outputs"""

    def __init__(self, file1, file2):
        self.file1 = file1
        self.file2 = file2

    def write(self, data):
        self.file1.write(data)
        self.file2.write(data)

    def flush(self):
        self.file1.flush()
        self.file2.flush()


###############################################################################
CODE_DOWNLOAD = """**Total running time of the script:**
({0:.0f} minutes {1:.3f} seconds)\n\n
\n.. container:: sphx-glr-download

    **Download Python source code:** :download:`{2} <{2}>`\n
\n.. container:: sphx-glr-download

    **Download IPython notebook:** :download:`{3} <{3}>`\n"""

# The following strings are used when we have several pictures: we use
# an html div tag that our CSS uses to turn the lists into horizontal
# lists.
HLIST_HEADER = """
.. rst-class:: sphx-glr-horizontal

"""

HLIST_IMAGE_TEMPLATE = """
    *

      .. image:: /%s
            :scale: 47
"""

SINGLE_IMAGE = """
.. image:: /%s
    :align: center
"""


CODE_OUTPUT = """.. rst-class:: sphx-glr-script-out

 Out::

  {0}\n"""


def get_docstring_and_rest(filename):
    """Separate `filename` content between docstring and the rest

    Strongly inspired from ast.get_docstring.

    Returns
    -------
    docstring: str
        docstring of `filename`
    rest: str
        `filename` content without the docstring
    """
    with open(filename) as f:
        content = f.read()

    node = ast.parse(content)
    if not isinstance(node, ast.Module):
        raise TypeError("This function only supports modules. "
                        "You provided {0}".format(node.__class__.__name__))
    if node.body and isinstance(node.body[0], ast.Expr) and \
       isinstance(node.body[0].value, ast.Str):
        docstring_node = node.body[0]
        docstring = docstring_node.value.s
        # This get the content of the file after the docstring last line
        # Note: 'maxsplit' argument is not a keyword argument in python2
        rest = content.split('\n', docstring_node.lineno)[-1]
        return docstring, rest
    else:
        raise ValueError(('Could not find docstring in file "{0}". '
                          'A docstring is required by sphinx-gallery')
                         .format(filename))


def split_code_and_text_blocks(source_file):
    """Return list with source file separated into code and text blocks.

    Returns
    -------
    blocks : list of (label, content)
        List where each element is a tuple with the label ('text' or 'code'),
        and content string of block.
    """
    docstring, rest_of_content = get_docstring_and_rest(source_file)

    blocks = [('text', docstring)]

    pattern = re.compile(
        r'(?P<header_line>^#{20,}.*)\s(?P<text_content>(?:^#.*\s)*)',
        flags=re.M)

    pos_so_far = 0
    for match in re.finditer(pattern, rest_of_content):
        match_start_pos, match_end_pos = match.span()
        code_block_content = rest_of_content[pos_so_far:match_start_pos]
        text_content = match.group('text_content')
        sub_pat = re.compile('^#', flags=re.M)
        text_block_content = dedent(re.sub(sub_pat, '', text_content))
        if code_block_content.strip():
            blocks.append(('code', code_block_content))
        if text_block_content.strip():
            blocks.append(('text', text_block_content))
        pos_so_far = match_end_pos

    remaining_content = rest_of_content[pos_so_far:]
    if remaining_content.strip():
        blocks.append(('code', remaining_content))

    return blocks


def codestr2rst(codestr, lang='python'):
    """Return reStructuredText code block from code string"""
    code_directive = "\n.. code-block:: {0}\n\n".format(lang)
    indented_block = indent(codestr, ' ' * 4)
    return code_directive + indented_block


def text2string(content):
    """Returns a string without the extra triple quotes"""
    try:
        return ast.literal_eval(content) + '\n'
    except Exception:
        return content


def extract_intro(filename):
    """ Extract the first paragraph of module-level docstring. max:95 char"""

    docstring, _ = get_docstring_and_rest(filename)

    # lstrip is just in case docstring has a '\n\n' at the beginning
    paragraphs = docstring.lstrip().split('\n\n')
    if len(paragraphs) > 1:
        first_paragraph = re.sub('\n', ' ', paragraphs[1])
        first_paragraph = (first_paragraph[:95] + '...'
                           if len(first_paragraph) > 95 else first_paragraph)
    else:
        raise ValueError(
            "Example docstring should have a header for the example title "
            "and at least a paragraph explaining what the example is about. "
            "Please check the example file:\n {}\n".format(filename))

    return first_paragraph


def get_md5sum(src_file):
    """Returns md5sum of file"""

    with open(src_file, 'r') as src_data:
        src_content = src_data.read()

        # data needs to be encoded in python3 before hashing
        if sys.version_info[0] == 3:
            src_content = src_content.encode('utf-8')

        src_md5 = hashlib.md5(src_content).hexdigest()
    return src_md5


def check_md5sum_change(src_file):
    """Returns True if src_file has a different md5sum"""

    src_md5 = get_md5sum(src_file)

    src_md5_file = src_file + '.md5'
    src_file_changed = True
    if os.path.exists(src_md5_file):
        with open(src_md5_file, 'r') as file_checksum:
            ref_md5 = file_checksum.read()
        if src_md5 == ref_md5:
            src_file_changed = False

    if src_file_changed:
        with open(src_md5_file, 'w') as file_checksum:
            file_checksum.write(src_md5)

    return src_file_changed


def _plots_are_current(src_file, image_file):
    """Test existence of image file and no change in md5sum of
    example"""

    first_image_file = image_file.format(1)
    has_image = os.path.exists(first_image_file)
    src_file_changed = check_md5sum_change(src_file)

    return has_image and not src_file_changed


def save_figures(image_path, fig_count, gallery_conf):
    """Save all open matplotlib figures of the example code-block

    Parameters
    ----------
    image_path : str
        Path where plots are saved (format string which accepts figure number)
    fig_count : int
        Previous figure number count. Figure number add from this number

    Returns
    -------
    list of strings containing the full path to each figure
    """
    figure_list = []

    fig_managers = matplotlib._pylab_helpers.Gcf.get_all_fig_managers()
    for fig_mngr in fig_managers:
        # Set the fig_num figure as the current figure as we can't
        # save a figure that's not the current figure.
        fig = plt.figure(fig_mngr.num)
        kwargs = {}
        to_rgba = matplotlib.colors.colorConverter.to_rgba
        for attr in ['facecolor', 'edgecolor']:
            fig_attr = getattr(fig, 'get_' + attr)()
            default_attr = matplotlib.rcParams['figure.' + attr]
            if to_rgba(fig_attr) != to_rgba(default_attr):
                kwargs[attr] = fig_attr

        current_fig = image_path.format(fig_count + fig_mngr.num)
        fig.savefig(current_fig, **kwargs)
        figure_list.append(current_fig)

    if gallery_conf.get('find_mayavi_figures', False):
        from mayavi import mlab
        e = mlab.get_engine()
        last_matplotlib_fig_num = len(figure_list)
        total_fig_num = last_matplotlib_fig_num + len(e.scenes)
        mayavi_fig_nums = range(last_matplotlib_fig_num, total_fig_num)

        for scene, mayavi_fig_num in zip(e.scenes, mayavi_fig_nums):
            current_fig = image_path.format(mayavi_fig_num)
            mlab.savefig(current_fig, figure=scene)
            # make sure the image is not too large
            scale_image(current_fig, current_fig, 850, 999)
            figure_list.append(current_fig)
        mlab.close(all=True)

    return figure_list


def scale_image(in_fname, out_fname, max_width, max_height):
    """Scales an image with the same aspect ratio centered in an
       image with a given max_width and max_height
       if in_fname == out_fname the image can only be scaled down
    """
    # local import to avoid testing dependency on PIL:
    try:
        from PIL import Image
    except ImportError:
        import Image
    img = Image.open(in_fname)
    width_in, height_in = img.size
    scale_w = max_width / float(width_in)
    scale_h = max_height / float(height_in)

    if height_in * scale_w <= max_height:
        scale = scale_w
    else:
        scale = scale_h

    if scale >= 1.0 and in_fname == out_fname:
        return

    width_sc = int(round(scale * width_in))
    height_sc = int(round(scale * height_in))

    # resize the image
    img.thumbnail((width_sc, height_sc), Image.ANTIALIAS)

    # insert centered
    thumb = Image.new('RGB', (max_width, max_height), (255, 255, 255))
    pos_insert = ((max_width - width_sc) // 2, (max_height - height_sc) // 2)
    thumb.paste(img, pos_insert)

    thumb.save(out_fname)
    # Use optipng to perform lossless compression on the resized image if
    # software is installed
    if os.environ.get('SKLEARN_DOC_OPTIPNG', False):
        try:
            subprocess.call(["optipng", "-quiet", "-o", "9", out_fname])
        except Exception:
            warnings.warn('Install optipng to reduce the size of the \
                          generated images')


def save_thumbnail(image_path, base_image_name, gallery_conf):
    """Save the thumbnail image"""
    first_image_file = image_path.format(1)
    thumb_dir = os.path.join(os.path.dirname(first_image_file), 'thumb')
    if not os.path.exists(thumb_dir):
        os.makedirs(thumb_dir)

    thumb_file = os.path.join(thumb_dir,
                              'sphx_glr_%s_thumb.png' % base_image_name)

    if os.path.exists(first_image_file):
        scale_image(first_image_file, thumb_file, 400, 280)
    elif not os.path.exists(thumb_file):
        # create something to replace the thumbnail
        default_thumb_file = os.path.join(glr_path_static(), 'no_image.png')
        default_thumb_file = gallery_conf.get("default_thumb_file",
                                              default_thumb_file)
        scale_image(default_thumb_file, thumb_file, 200, 140)


def generate_dir_rst(src_dir, target_dir, gallery_conf, seen_backrefs):
    """Generate the gallery reStructuredText for an example directory"""
    if not os.path.exists(os.path.join(src_dir, 'README.txt')):
        print(80 * '_')
        print('Example directory %s does not have a README.txt file' %
              src_dir)
        print('Skipping this directory')
        print(80 * '_')
        return ""  # because string is an expected return type

    fhindex = open(os.path.join(src_dir, 'README.txt')).read()
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    sorted_listdir = [fname for fname in sorted(os.listdir(src_dir))
                      if fname.endswith('.py')]
    entries_text = []
    for fname in sorted_listdir:
        amount_of_code = generate_file_rst(fname, target_dir, src_dir,
                                           gallery_conf)
        new_fname = os.path.join(src_dir, fname)
        intro = extract_intro(new_fname)
        write_backreferences(seen_backrefs, gallery_conf,
                             target_dir, fname, intro)
        this_entry = _thumbnail_div(target_dir, fname, intro) + """

.. toctree::
   :hidden:

   /%s/%s\n""" % (target_dir, fname[:-3])
        entries_text.append((amount_of_code, this_entry))

    # sort to have the smallest entries in the beginning
    entries_text.sort()

    for _, entry_text in entries_text:
        fhindex += entry_text

    # clear at the end of the section
    fhindex += """.. raw:: html\n
    <div style='clear:both'></div>\n\n"""

    return fhindex


def execute_script(code_block, example_globals, image_path, fig_count,
                   src_file, gallery_conf):
    """Executes the code block of the example file"""
    time_elapsed = 0
    stdout = ''

    # We need to execute the code
    print('plotting code blocks in %s' % src_file)

    plt.close('all')
    cwd = os.getcwd()
    # Redirect output to stdout and
    orig_stdout = sys.stdout

    try:
        # First cd in the original example dir, so that any file
        # created by the example get created in this directory
        os.chdir(os.path.dirname(src_file))
        my_buffer = StringIO()
        my_stdout = Tee(sys.stdout, my_buffer)
        sys.stdout = my_stdout

        t_start = time()
        exec(code_block, example_globals)
        time_elapsed = time() - t_start

        sys.stdout = orig_stdout

        my_stdout = my_buffer.getvalue().strip().expandtabs()
        if my_stdout:
            stdout = CODE_OUTPUT.format(indent(my_stdout, ' ' * 4))
        os.chdir(cwd)
        figure_list = save_figures(image_path, fig_count, gallery_conf)

        # Depending on whether we have one or more figures, we're using a
        # horizontal list or a single rst call to 'image'.
        image_list = ""
        if len(figure_list) == 1:
            figure_name = figure_list[0]
            image_list = SINGLE_IMAGE % figure_name.lstrip('/')
        elif len(figure_list) > 1:
            image_list = HLIST_HEADER
            for figure_name in figure_list:
                image_list += HLIST_IMAGE_TEMPLATE % figure_name.lstrip('/')

    except Exception:
        formatted_exception = traceback.format_exc()

        print(80 * '_')
        print('%s is not compiling:' % src_file)
        print(formatted_exception)
        print(80 * '_')

        figure_list = []
        image_list = codestr2rst(formatted_exception, lang='pytb')

        # Overrides the output thumbnail in the gallery for easy identification
        broken_img = os.path.join(glr_path_static(), 'broken_example.png')
        shutil.copyfile(broken_img, os.path.join(cwd, image_path.format(1)))
        fig_count += 1  # raise count to avoid overwriting image

        # Breaks build on first example error

        if gallery_conf['abort_on_example_error']:
            raise

    finally:
        os.chdir(cwd)
        sys.stdout = orig_stdout

    print(" - time elapsed : %.2g sec" % time_elapsed)
    code_output = "\n{0}\n\n{1}\n\n".format(image_list, stdout)

    return code_output, time_elapsed, fig_count + len(figure_list)


def generate_file_rst(fname, target_dir, src_dir, gallery_conf):
    """ Generate the rst file for a given example.

        Returns the amout of code (in characters) of the corresponding
        files.
    """

    src_file = os.path.join(src_dir, fname)
    example_file = os.path.join(target_dir, fname)
    shutil.copyfile(src_file, example_file)

    image_dir = os.path.join(target_dir, 'images')
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)

    base_image_name = os.path.splitext(fname)[0]
    image_fname = 'sphx_glr_' + base_image_name + '_{0:03}.png'
    image_path = os.path.join(image_dir, image_fname)

    script_blocks = split_code_and_text_blocks(example_file)

    amount_of_code = sum([len(bcontent)
                          for blabel, bcontent in script_blocks
                          if blabel == 'code'])

    if _plots_are_current(example_file, image_path):
        return amount_of_code

    time_elapsed = 0

    ref_fname = example_file.replace(os.path.sep, '_')
    example_rst = """\n\n.. _sphx_glr_{0}:\n\n""".format(ref_fname)
    example_nb = Notebook(fname, target_dir)

    filename_pattern = gallery_conf.get('filename_pattern')
    if re.search(filename_pattern, src_file) and gallery_conf['plot_gallery']:
        example_globals = {
            # A lot of examples contains 'print(__doc__)' for example in
            # scikit-learn so that running the example prints some useful
            # information. Because the docstring has been separated from
            # the code blocks in sphinx-gallery, __doc__ is actually
            # __builtin__.__doc__ in the execution context and we do not
            # want to print it
            '__doc__': '',
            # Examples may contain if __name__ == '__main__' guards
            # for in example scikit-learn if the example uses multiprocessing
            '__name__': '__main__'}

        fig_count = 0
        # A simple example has two blocks: one for the
        # example introduction/explanation and one for the code
        is_example_notebook_like = len(script_blocks) > 2
        for blabel, bcontent in script_blocks:
            if blabel == 'code':
                code_output, rtime, fig_count = execute_script(bcontent,
                                                               example_globals,
                                                               image_path,
                                                               fig_count,
                                                               src_file,
                                                               gallery_conf)

                time_elapsed += rtime
                example_nb.add_code_cell(bcontent)

                if is_example_notebook_like:
                    example_rst += codestr2rst(bcontent) + '\n'
                    example_rst += code_output
                else:
                    example_rst += code_output
                    if 'sphx-glr-script-out' in code_output:
                        # Add some vertical space after output
                        example_rst += "\n\n|\n\n"
                    example_rst += codestr2rst(bcontent) + '\n'

            else:
                example_rst += text2string(bcontent) + '\n'
                example_nb.add_markdown_cell(text2string(bcontent))
    else:
        for blabel, bcontent in script_blocks:
            if blabel == 'code':
                example_rst += codestr2rst(bcontent) + '\n'
                example_nb.add_code_cell(bcontent)
            else:
                example_rst += bcontent + '\n'
                example_nb.add_markdown_cell(text2string(bcontent))

    save_thumbnail(image_path, base_image_name, gallery_conf)

    time_m, time_s = divmod(time_elapsed, 60)
    example_nb.save_file()
    with open(os.path.join(target_dir, base_image_name + '.rst'), 'w') as f:
        example_rst += CODE_DOWNLOAD.format(time_m, time_s, fname,
                                            example_nb.file_name)
        f.write(example_rst)

    return amount_of_code
